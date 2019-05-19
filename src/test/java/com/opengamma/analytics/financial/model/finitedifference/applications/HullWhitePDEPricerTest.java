/**
 * Copyright (C) 2012 - present by OpenGamma Inc. and the OpenGamma group of companies
 * 
 * Please see distribution for license.
 */
package com.opengamma.analytics.financial.model.finitedifference.applications;

import org.testng.annotations.Test;

import com.opengamma.analytics.financial.model.finitedifference.ExponentialMeshing;
import com.opengamma.analytics.financial.model.finitedifference.MeshingFunction;
import com.opengamma.analytics.financial.model.finitedifference.PDEGrid1D;
import com.opengamma.analytics.financial.model.finitedifference.SandBox;
import com.opengamma.analytics.financial.model.option.pricing.analytic.BjerksundStenslandModel;
import com.opengamma.util.test.TestGroup;

/**
 * Test.
 */
@Test(groups = TestGroup.UNIT_SLOW)
public class HullWhitePDEPricerTest {
  private static final BjerksundStenslandModel AMERICAN_APPOX_PRCIER = new BjerksundStenslandModel();
  private static final BlackScholesMertonPDEPricer PRICER = new BlackScholesMertonPDEPricer();
  private static final HullWhitePDEPricer PRICERHW = new HullWhitePDEPricer();

  /**
   * This tests that the run time scales as the grid size for moderately sized grids. The grid size is the total number of nodes in the grid (#space nodes x #time nodes)
   * and moderately sized means between 1e5 and 1e9 nodes (e.g. from 100 time nodes and 1000 space nodes to 1e4 time and 1e5 space nodes).
   * We also test that the relative error is inversely proportional to the grid size; the constant of proportionality depends on nu (the ration of space nodes to time nodes) and
   * is lowest for nu = 125 in this case.
   */
  //aaafinal double r0,  final double t,final double rMax,final double rMin, final double sigma, final double thetaHW,final double kappa, final PDEGrid1D[] grid, final double[] theta
  
  @Test
  public void HullWhiteTest() {
	  final double sigma = 0.000212693; 
	  final int T = 111; 
	  final double rMax = .01; 
	  final double rMin = -.01; 
	  final double thetaHW = 0.000023637253287033743; 
	  final double kappa = 0.008622227075117428; 
	  final double r0 = -0.00024031738413887727;
	  //final double r0 = .010191919;
	  
    final double beta = 0.01;
    final double lambda = 0.0;
    final double sd = 6.0;

    final int tSteps = 200;
    final double nu = 10;
    final int sSteps = (int) (nu * tSteps);

    //final double bsPrice = Math.exp(-r * t) * BlackFormulaRepository.price(s0 * Math.exp(b * t), k, t, sigma, isCall);
    final MeshingFunction xMesh = new ExponentialMeshing(rMin, rMax, sSteps,0, new double[] {r0});
    PDEGrid1D[] grid = new PDEGrid1D[1];
    final MeshingFunction tMesh = new ExponentialMeshing(0, T, tSteps, -2.5);
    grid[0] = new PDEGrid1D(tMesh, xMesh);
    double[] theta = new double[] {1};
    final double pdePrice1 = PRICERHW.price(r0, T, rMax,rMin, sigma, thetaHW,kappa, grid, theta);
    //final double pdePrice2 = PRICER.price(s0, k, r, b, t, sigma, isCall, false, sSteps, tSteps, beta, lambda, sd);
    //System.out.println(pdePrice1*SandBox.price(0d, (double) T, r0));
    System.out.println(pdePrice1);
    //System.out.println(100*SandBox.price(0d, (double)T, r0));
    //assertEquals(0, relErr1, 5e-4);
    //assertEquals(0, relErr2, 2e-6); // much better accuracy with non-uniform
  }
}
  
// }