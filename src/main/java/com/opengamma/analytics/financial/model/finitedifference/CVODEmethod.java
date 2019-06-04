/**
 * Copyright (C) 2009 - present by OpenGamma Inc. and the OpenGamma group of companies
 *
 * Please see distribution for license.
 */
package com.opengamma.analytics.financial.model.finitedifference;

import static de.grogra.numeric.cvode.CVODE.CVDense;
import static de.grogra.numeric.cvode.CVODE.CV_BDF;
import static de.grogra.numeric.cvode.CVODE.CV_NEWTON;
import static de.grogra.numeric.cvode.CVODE.CV_NORMAL;
import static de.grogra.numeric.cvode.CVODE.CV_ADAMS;
import static de.grogra.numeric.cvode.CVODE.CV_ONE_STEP;
import static de.grogra.numeric.cvode.CVODE.CV_ROOT_RETURN;
import static de.grogra.numeric.cvode.CVODE.CV_SUCCESS;
import static de.grogra.numeric.cvode.CVODE.CVode;
import static de.grogra.numeric.cvode.CVODE.CVodeGetRootInfo;
import static de.grogra.numeric.cvode.CVODE.CVodeInit;
import static de.grogra.numeric.cvode.CVODE.CVodeRootInit;
import static de.grogra.numeric.cvode.CVODE.CVodeSStolerances;
import static de.grogra.numeric.cvode.CVODE.CVodeSVtolerances;
import static de.grogra.numeric.cvode.N_Vec_Serial.N_VNew_Serial;
import static de.grogra.numeric.cvode.CVODE.CVodeCreate;

import com.opengamma.analytics.math.linearalgebra.Decomposition;
import com.opengamma.analytics.math.linearalgebra.LUDecompositionCommons;
import com.opengamma.analytics.math.linearalgebra.TridiagonalMatrix;
import com.opengamma.analytics.math.matrix.DoubleMatrix1D;
import com.opengamma.analytics.math.matrix.DoubleMatrix2D;
import com.opengamma.analytics.math.matrix.OGMatrixAlgebra;
import com.opengamma.analytics.math.surface.Surface;
import com.opengamma.util.ArgumentChecker;
import com.sun.jna.NativeLong;
import com.sun.jna.Pointer;
import com.sun.jna.Structure;
import com.sun.jna.ptr.DoubleByReference;

import de.grogra.numeric.Monitor;
import de.grogra.numeric.cvode.CVRhsFn;
import de.grogra.numeric.cvode.CVRootFn;
import de.grogra.numeric.cvode.N_Vector;

/**
 * A theta (i.e. weighted between explicit and implicit time stepping) scheme using SOR algorithm to solve the matrix system at each time step
 * This uses the exponentially fitted scheme of duffy
 */
public class CVODEmethod implements ConvectionDiffusionPDESolver {
  private static final Decomposition<?> DCOMP = new LUDecompositionCommons();
  // private static final DEFAULT
  private final double _theta;
  private final boolean _showFullResults;

  /**
   * Sets up a standard Crank-Nicolson scheme
   */
  public CVODEmethod() {
    _theta = 0.5;
    _showFullResults = false;
  }

  /**
   * Sets up a scheme that is the weighted average of an explicit and an implicit scheme
   * @param theta The weight. theta = 0 - fully explicit, theta = 0.5 - Crank-Nicolson, theta = 1.0 - fully implicit
   * @param showFullResults Show the full results
   */
  public CVODEmethod(final double theta, final boolean showFullResults, final int N) {
    ArgumentChecker.isTrue(theta >= 0 && theta <= 1.0, "theta must be in the range 0 to 1");
    _theta = theta;
    _showFullResults = showFullResults;
    
  }
  

  public double getTheta() {
    return _theta;
  }

  @Override
  //TODO This is so ugly
  public PDEResults1D solve(final PDE1DDataBundle<ConvectionDiffusionPDE1DCoefficients> pdeData) {
    ArgumentChecker.notNull(pdeData, "pde data");
    final ConvectionDiffusionPDE1DCoefficients coeff = pdeData.getCoefficients();
    if (coeff instanceof ConvectionDiffusionPDE1DStandardCoefficients) {
      final PDE1DDataBundle<ConvectionDiffusionPDE1DStandardCoefficients> temp = convertPDE1DDataBundle(pdeData);
      final SolverImpl solver = new SolverImpl(temp);
      return solver.res;
      //return solver.solve();
    } /*else if (coeff instanceof ConvectionDiffusionPDE1DFullCoefficients) {
      final ConvectionDiffusionPDE1DFullCoefficients temp = (ConvectionDiffusionPDE1DFullCoefficients) coeff;
      final ExtendedSolverImpl solver = new ExtendedSolverImpl(temp, pdeData.getInitialCondition(), pdeData.getLowerBoundary(), pdeData.getUpperBoundary(),
          pdeData.getFreeBoundary(), pdeData.getGrid());
      return solver.solve();
    }*/
    //throw new IllegalArgumentException(coeff.getClass() + " not handled");
	return null;
  }

  private static PDE1DDataBundle<ConvectionDiffusionPDE1DStandardCoefficients> convertPDE1DDataBundle(final PDE1DDataBundle<ConvectionDiffusionPDE1DCoefficients> pdeData) {
    if (pdeData.getFreeBoundary() == null) {
      return new PDE1DDataBundle<>(
          (ConvectionDiffusionPDE1DStandardCoefficients) pdeData.getCoefficients(), pdeData.getInitialCondition(), pdeData.getLowerBoundary(),
          pdeData.getUpperBoundary(), pdeData.getGrid());
    }
    return new PDE1DDataBundle<>(
        (ConvectionDiffusionPDE1DStandardCoefficients) pdeData.getCoefficients(), pdeData.getInitialCondition(), pdeData.getLowerBoundary(),
        pdeData.getUpperBoundary(), pdeData.getFreeBoundary(), pdeData.getGrid());
  }

  private enum SolverMode {
    tridiagonal,
    luDecomp,
    psor;
  }

  class SolverImpl {

    // grid
    private final int _nNodesX;
    private final int _nNodesT;
    private final double[] _dt;
    private final PDEGrid1D _grid;
    private final double[][] _x1st;
    private final double[][] _x2nd;
    private final double[] _dx;
    private PDEFullResults1D res;
    //initial and boundary conditions
    private final double[] _initial;
    private final BoundaryCondition _lower;
    private final BoundaryCondition _upper;
    //PDE coefficients
    private final ConvectionDiffusionPDE1DStandardCoefficients _coeff;
    //free boundary problems
    private final SolverMode _mode;
    private final Surface<Double, Double, Double> _freeB;

    public SolverImpl(final PDE1DDataBundle<ConvectionDiffusionPDE1DStandardCoefficients> pdeData) {

      //unpack pdeData
    	
      _grid = pdeData.getGrid();
      _coeff = pdeData.getCoefficients();
      _lower = pdeData.getLowerBoundary();
      _upper = pdeData.getUpperBoundary();

      _nNodesX = _grid.getNumSpaceNodes();
      _nNodesT = _grid.getNumTimeNodes();
      final int N =_nNodesX; 
      _x1st = new double[_nNodesX - 2][];
      _x2nd = new double[_nNodesX - 2][];
      for (int ii = 0; ii < _nNodesX - 2; ii++) {
        _x1st[ii] = _grid.getFirstDerivativeCoefficients(ii + 1);
        _x2nd[ii] = _grid.getSecondDerivativeCoefficients(ii + 1);
      }
      _dx = new double[_nNodesX - 1];
      for (int ii = 0; ii < _nNodesX - 1; ii++) {
        _dx[ii] = _grid.getSpaceStep(ii);
      }

      _initial = pdeData.getInitialCondition();
      _dt = new double[_nNodesT - 1];
      for (int jj = 0; jj < _nNodesT - 1; jj++) {
        _dt[jj] = _grid.getTimeStep(jj);
      }

      //free boundary
      _freeB = pdeData.getFreeBoundary();
      if (_freeB == null) {
        _mode = SolverMode.tridiagonal;
      } else {
        _mode = SolverMode.psor;
      }
      int flag;
		N_Vector y = null;
		Pointer cvode_mem = Pointer.NULL;
      y = N_VNew_Serial(new NativeLong(N));
      
      cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
      
      CVRhsFn f = new CVRhsFn() {
  		final double[] state = new double[N];
  		final DoubleMatrix1D rate = new DoubleMatrix1D(N);
  		@Override
  		public int callback(double t, N_Vector y, N_Vector ydot,
  				Pointer user_data) {
  			try {
  				assert N == y.getLength();
  				assert N == ydot.getLength();
  				y.get(state);
  				computeDerivatives(t,state,rate);
  				ydot.set(rate.toArray());
  				return 0;
  			} catch (Throwable throwable) {
  				throwable.printStackTrace();
  			}
  			// return negative value to indicate unrecoverable error
  			return -1;
  		}
      };
  	CVRootFn g = new CVRootFn() {
  			final double[] out = new double[nrtfn];
  			final double[] y = new double[N];
  			@Override
  			public int callback(double t, N_Vector y, Pointer gout,
  					Pointer user_data) {
  				try {
  					assert N == y.getLength();
  					y.get(this.y);
  					monitor.g(out, t, this.y);
  					gout.write(0, out, 0, nrtfn);
  					return 0;
  				} catch (Throwable throwable) {
  					throwable.printStackTrace();
  				}
  				// return non-zero value to indicate error
  				return -1;
  			}
  		};
  

  		
      double[][] full = null;
      full = new double[_nNodesT][_nNodesX];
      if (_showFullResults) {
        
        //full[0] = _initial;
      }
      double[] y0 = _initial;
      y.set(y0);
      double t0 = _grid.getTimeNode(0);
      
  	flag = CVodeInit(cvode_mem, f, t0, y);
  	if (absTol != null) {
		// relative is scalar, absolute is vector 
		assert absTolDefault > 0;
		assert relTolDefault > 0;
		assert absTol.length == N;
		y.set(absTol);
		for (int i = 0; i < N; i++) {
			// replace zero tolerance by default one
			if (absTol[i] == 0)
				y.set(i, absTolDefault);
		}
		flag = CVodeSVtolerances(cvode_mem, relTolDefault, y);
	} else {
		// relative is scalar, absolute is scalar
		flag = CVodeSStolerances(cvode_mem, relTolDefault, absTolDefault);
	}
  	flag = CVDense(cvode_mem, N);
    CVodeRootInit(cvode_mem, nrtfn, g);
  	
  	//assert y0.length == y1.length;
	while (true) {
		// perform actual integration
		DoubleByReference t = new DoubleByReference();
		double t1 = _grid.getTimeNode(_nNodesT-1);
		long startTime = System.currentTimeMillis();
		//int flag; 
		
		
		flag = CVode(cvode_mem, t1, y, t, CV_NORMAL);
		System.out.println(System.currentTimeMillis() - startTime);
		//System.out.println("t = " + t1.getValue()  + " t1 " + t1);
		final double[] resA = new double[N];
		if (flag == CV_SUCCESS) {
			// break loop if target time was reached
			//System.out.println("t = " + flag);
			System.out.println("y = " + y.get(1));
			
			y.get(resA);
			//full[0] = _initial;
			full[_nNodesT-1]= resA;
			res = new PDEFullResults1D(_grid, full)
			;
			
			
			break;
		} else if (flag == CV_ROOT_RETURN) {
			// one of the monitor functions triggered
			// find out which
			int[] rootsfound = new int[nrtfn];
			flag = CVodeGetRootInfo(cvode_mem, rootsfound);
		/* checkFlag(flag, "root was found, but could not determine which"); */
			boolean stop = false;
			for (int i = 0; i < nrtfn; i++) {
				if (rootsfound[i] != 0) {
					// call event handler
					//stop |= monitor.handleEvent(i, t.getValue(), y1);
				}
			}
			if (stop) {
				break;
			}
		} else {
		/* checkFlag(flag, "error during integration"); */
		}
	}
	
			/*
			 * if (_showFullResults) { res = new PDEFullResults1D(_grid, full); } else { res
			 * = new PDETerminalResults1D(_grid, h); }
			 */
    //return res;
	//return y1;
    }
  
  public void computeDerivatives(double t, double[] y,DoubleMatrix1D rate) {
      //for (int jj = 0; jj < _nNodesT - 1; jj++) {
/*        final double dt = _dt[jj];
        //double t = _grid.getTimeNode(0);

        double[] topRow = _lower.getLeftMatrixCondition(_coeff, _grid, t);
        double[] bottomRow = _upper.getLeftMatrixCondition(_coeff, _grid, t);
        final double[] cDag = new double[_nNodesX - 2];
        final double[] lDag = new double[_nNodesX - 2];
        final double[] uDag = new double[_nNodesX - 2];
        for (int ii = 0; ii < _nNodesX - 2; ii++) { //tri-diagonal form
          final double x = _grid.getSpaceNode(ii + 1);
          final double a = _coeff.getA(t, x);
          final double b = _coeff.getB(t, x);
          final double c = _coeff.getC(t, x);
          
          cDag[ii] = _x2nd[ii][1] * a + _x1st[ii][1] * b + c;
          lDag[ii] = _x2nd[ii][0] * a + _x1st[ii][0] * b;
          uDag[ii] = _x2nd[ii][2] * a + _x1st[ii][2] * b;
        }
        //RHS of system
        final double[] y = new double[_nNodesX];
        //main part of RHS
        for (int ii = 1; ii < _nNodesX - 1; ii++) { //tri-diagonal form
          //y[ii] = (1 - (1 - _theta) * dt * cDag[ii - 1]) * h[ii] - (1 - _theta) * dt * (lDag[ii - 1] * h[ii - 1] + +uDag[ii - 1] * h[ii + 1]);
          y[ii]=h[ii]+_coeff.getD(_grid.getTimeNode(jj),_grid.getSpaceNode(ii))*dt;
        }

        double t = _grid.getTimeNode(jj + 1);
        //lower & upper boundaries
        //y[0] = _lower.getConstant(_coeff, t);
        //y[_nNodesX - 1] = _upper.getConstant(_coeff, t);

        //put the LHS of system in tri-diagonal form
        final double[] d = new double[_nNodesX]; //main diag
        final double[] u = new double[_nNodesX - 1]; //upper
        final double[] l = new double[_nNodesX - 1]; //lower
        //lower boundary conditions
        topRow = _lower.getLeftMatrixCondition(_coeff, _grid, t);
        final int p2 = topRow.length;
        d[0] = topRow[0];
        if (p2 > 1) {
          u[0] = topRow[1];
          //Review do we need this?
          ArgumentChecker.isFalse(p2 > 2, "Boundary condition means that system is not tri-diagonal");
        }
        bottomRow = _upper.getLeftMatrixCondition(_coeff, _grid, t);
        final int q2 = bottomRow.length;
        d[_nNodesX - 1] = bottomRow[q2 - 1];
        if (q2 > 1) {
          l[_nNodesX - 2] = bottomRow[q2 - 2];
          ArgumentChecker.isFalse(q2 > 2, "Boundary condition means that system is not tri-diagonal");
        }
*/	
	  final double[] d = new double[_nNodesX]; //main diag
      final double[] u = new double[_nNodesX - 1]; //upper
      final double[] l = new double[_nNodesX - 1]; //lower
      double[] topRow = _lower.getLeftMatrixCondition(_coeff, _grid, t);
      double[] bottomRow = _upper.getLeftMatrixCondition(_coeff, _grid, t);
      final int q2 = bottomRow.length;
      d[_nNodesX - 1] = bottomRow[q2 - 1];
      if (q2 > 1) {
        l[_nNodesX - 2] = bottomRow[q2 - 2];
        ArgumentChecker.isFalse(q2 > 2, "Boundary condition means that system is not tri-diagonal");
      }
      final int p2 = topRow.length;
      d[0] = topRow[0];
      if (p2 > 1) {
        u[0] = topRow[1];
        //Review do we need this?
        ArgumentChecker.isFalse(p2 > 2, "Boundary condition means that system is not tri-diagonal");
      }
	  final double[] cDag = new double[_nNodesX - 2];
      final double[] lDag = new double[_nNodesX - 2];
      final double[] uDag = new double[_nNodesX - 2];  
	  for (int ii = 0; ii < _nNodesX - 2; ii++) { //tri-diagonal form
          final double x = _grid.getSpaceNode(ii + 1);
          final double a = _coeff.getA(t, x);
          final double b = _coeff.getB(t, x);
          final double c = _coeff.getC(t, x);
          //debug - fitting par
          //a = getFittingParameter(a, b, ii);
          d[ii] = _x2nd[ii][1] * a + _x1st[ii][1] * b + c;
          l[ii] = _x2nd[ii][0] * a + _x1st[ii][0] * b;
          u[ii] = _x2nd[ii][2] * a + _x1st[ii][2] * b;
        }

			         
	  OGMatrixAlgebra ma= new OGMatrixAlgebra();
	  DoubleMatrix1D y1 = new DoubleMatrix1D(y);
	  final DoubleMatrix2D lhs = new TridiagonalMatrix(d, u, l).toDoubleMatrix2D();
	  rate = (DoubleMatrix1D) ma.multiply(lhs, y1);
        //solve the system (update h)
     int k=1;
    }

   
}
  double absTolDefault = 1e-3;
	double relTolDefault = 1e-3;
	
	double[] absTol;
//	double[] relTol;
	
	int nrtfn;
	Monitor monitor;
	
	
		  
}
