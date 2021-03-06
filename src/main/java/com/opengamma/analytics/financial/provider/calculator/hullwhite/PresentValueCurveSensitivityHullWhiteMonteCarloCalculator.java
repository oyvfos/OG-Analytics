/**
 * Copyright (C) 2011 - present by OpenGamma Inc. and the OpenGamma group of companies
 * 
 * Please see distribution for license.
 */
package com.opengamma.analytics.financial.provider.calculator.hullwhite;

import cern.jet.random.engine.MersenneTwister;

import com.opengamma.analytics.financial.interestrate.InstrumentDerivativeVisitorAdapter;
import com.opengamma.analytics.financial.interestrate.annuity.derivative.AnnuityCouponIborRatchet;
import com.opengamma.analytics.financial.montecarlo.provider.HullWhiteMonteCarloMethod;
import com.opengamma.analytics.financial.provider.description.interestrate.HullWhiteOneFactorProviderInterface;
import com.opengamma.analytics.financial.provider.sensitivity.multicurve.MultipleCurrencyMulticurveSensitivity;
import com.opengamma.analytics.math.random.NormalRandomNumberGenerator;

/**
 * Calculates the present value of an inflation instruments by discounting for a given MarketBundle
 */
public class PresentValueCurveSensitivityHullWhiteMonteCarloCalculator extends InstrumentDerivativeVisitorAdapter<HullWhiteOneFactorProviderInterface, MultipleCurrencyMulticurveSensitivity> {

  /**
   * The default number of path in the Monte Carlo simulation.
   */
  private static final int DEFAULT_NB_PATH = 12500;
  /**
   * The number of paths used in the simulation.
   */
  private final int _nbPath;

  /**
   * Calculator constructor using the default number of paths.
   */
  public PresentValueCurveSensitivityHullWhiteMonteCarloCalculator() {
    _nbPath = DEFAULT_NB_PATH;
  }

  /**
   * Constructor with a given number of simulation paths.
   * @param nbPath The number of paths.
   */
  public PresentValueCurveSensitivityHullWhiteMonteCarloCalculator(final int nbPath) {
    _nbPath = nbPath;
  }

  @Override
  public MultipleCurrencyMulticurveSensitivity visitAnnuityCouponIborRatchet(final AnnuityCouponIborRatchet annuity, final HullWhiteOneFactorProviderInterface hullWhite) {
    HullWhiteMonteCarloMethod methodMC = new HullWhiteMonteCarloMethod(new NormalRandomNumberGenerator(0.0, 1.0, new MersenneTwister()), _nbPath);
    return methodMC.presentValueCurveSensitivity(annuity, annuity.getCurrency(), hullWhite);
  }

}
