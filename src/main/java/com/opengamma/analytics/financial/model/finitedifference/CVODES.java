package com.opengamma.analytics.financial.model.finitedifference;

import static de.grogra.numeric.cvode.CVODE.CVDense;
import static de.grogra.numeric.cvode.CVODE.CV_BDF;
import static de.grogra.numeric.cvode.CVODE.CV_NEWTON;
import static de.grogra.numeric.cvode.CVODE.CV_NORMAL;
import static de.grogra.numeric.cvode.CVODE.CV_ROOT_RETURN;
import static de.grogra.numeric.cvode.CVODE.CV_SUCCESS;
import static de.grogra.numeric.cvode.CVODE.CVode;
import static de.grogra.numeric.cvode.CVODE.CVodeCreate;
import static de.grogra.numeric.cvode.CVODE.CVodeGetRootInfo;
import static de.grogra.numeric.cvode.CVODE.CVodeInit;
import static de.grogra.numeric.cvode.CVODE.CVodeRootInit;
import static de.grogra.numeric.cvode.CVODE.CVodeSStolerances;
import static de.grogra.numeric.cvode.CVODE.CVodeSVtolerances;
import static de.grogra.numeric.cvode.N_Vec_Serial.N_VNew_Serial;

import org.apache.commons.math.ode.DerivativeException;
import org.apache.commons.math.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math.ode.IntegratorException;

import com.sun.jna.NativeLong;
import com.sun.jna.Pointer;
import com.sun.jna.ptr.DoubleByReference;

import de.grogra.numeric.Monitor;
import de.grogra.numeric.cvode.CVRhsFn;
import de.grogra.numeric.cvode.CVRootFn;
import de.grogra.numeric.cvode.N_Vector;

public class CVODES {
	double absTolDefault = 1e-10;
	double relTolDefault = 1e-10;
	
	double[] absTol;
//	double[] relTol;
	
	int nrtfn;
	Monitor monitor;
	Pointer cvode_mem = Pointer.NULL;
	
	public CVODES(final int N) {
		    int flag = CVDense(cvode_mem, N);
		    CVRhsFn f = new CVRhsFn() {
				final double[] state = new double[N];
				final double[] rate = new double[N];
				@Override
				public int callback(double t, N_Vector y, N_Vector ydot,
						Pointer user_data) {
					try {
						assert N == y.getLength();
						assert N == ydot.getLength();
						y.get(state);
						//ode.computeDerivatives(t,state,rate);
						ydot.set(rate);
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
			
			
			
			/* checkFlag(flag, "could not set linear solver"); */
				
				// set linear solver optional inputs
				// TODO
				
				// specify rootfinding problem
			CVodeRootInit(cvode_mem, nrtfn, g);
		   // _showFullResults = false;
		
		  }
	private double[] solveCVODES(double t0, double[] y0, double t1, double[] y1)
				throws DerivativeException, IntegratorException {
			assert y0.length == y1.length;
					while (true) {
			// perform actual integration
			DoubleByReference t = new DoubleByReference();
			long startTime = System.currentTimeMillis();
			int flag; 
			
			flag = CVode(cvode_mem, t1, y, t, CV_NORMAL);
			System.out.println(System.currentTimeMillis() - startTime);
			System.out.println("t = " + t.getValue()  + " t1 " + t1);
			System.out.println("y = " + y.get(1));
			if (flag == CV_SUCCESS) {
				// break loop if target time was reached
				System.out.println("t = " + flag);
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
						stop |= monitor.handleEvent(i, t.getValue(), y1);
					}
				}
				if (stop) {
					break;
				}
			} else {
			/* checkFlag(flag, "error during integration"); */
			}
		}
		return y1;
	    }
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
