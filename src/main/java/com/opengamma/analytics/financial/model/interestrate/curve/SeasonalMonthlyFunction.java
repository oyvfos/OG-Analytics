/**
 * Copyright (C) 2013 - present by OpenGamma Inc. and the OpenGamma group of companies
 *
 * Please see distribution for license.
 */
package com.opengamma.analytics.financial.model.interestrate.curve;

import java.util.Map;
import java.util.Set;

import org.apache.commons.lang.Validate;
import org.joda.beans.Bean;
import org.joda.beans.BeanBuilder;
import org.joda.beans.JodaBeanUtils;
import org.joda.beans.MetaProperty;
import org.joda.beans.Property;
import org.joda.beans.PropertyDefinition;
import org.joda.beans.impl.direct.DirectMetaBean;
import org.joda.beans.impl.direct.DirectMetaProperty;
import org.joda.beans.impl.direct.DirectMetaPropertyMap;

import com.opengamma.analytics.math.function.Function1D;
import com.opengamma.analytics.math.interpolation.CombinedInterpolatorExtrapolatorFactory;
import com.opengamma.analytics.math.interpolation.Interpolator1D;
import com.opengamma.analytics.math.interpolation.Interpolator1DFactory;
import com.opengamma.analytics.math.interpolation.data.Interpolator1DDataBundle;

/**
 * Monthly seasonal adjustments (used in particular for inflation curves).
 */
class SeasonalMonthlyFunction extends Function1D<Double, Double> implements Bean {

  /**
   * The cumulative multiplicative seasonal factors from the reference time to the next. 
   * Array of size multiple of 12.
   */
  @PropertyDefinition(get = "private")
  private final double[] _monthlyCumulativeFactors;

  /**
   * The steps in the monthly seasonal adjustments.
   */
  @PropertyDefinition(get = "private")
  private final double[] _steps;

  /**
   * The number of month in a year.
   */
  private static final int NB_MONTH = 12;
  private static final Interpolator1D INTERPOLATOR_STEP_FLAT =
      CombinedInterpolatorExtrapolatorFactory.getInterpolator(Interpolator1DFactory.STEP,
          Interpolator1DFactory.FLAT_EXTRAPOLATOR, Interpolator1DFactory.FLAT_EXTRAPOLATOR);

  /**
   * Constructor for the monthly seasonal adjustments.
   * @param steps
   * @param monthlyFactors
   * @param isAdditive
   */
  public SeasonalMonthlyFunction(double[] steps, double[] monthlyFactors, boolean isAdditive) {
    Validate.notNull(monthlyFactors, "Monthly factors");
    Validate.notNull(steps, "steps");
    Validate.isTrue(monthlyFactors.length == NB_MONTH - 1, "Monthly factors with incorrect length; should be 11");
    Validate.notNull(isAdditive, "isAdditive");
    _steps = steps;
    double[] cumulativeFactors = new double[NB_MONTH]; // monthly factors
    if (isAdditive) {
      cumulativeFactors[0] = 0.0;
      for (int loopmonth = 1; loopmonth < NB_MONTH; loopmonth++) {
        cumulativeFactors[loopmonth] = cumulativeFactors[loopmonth - 1] + monthlyFactors[loopmonth - 1];
      }
    } else {
      cumulativeFactors[0] = 1.0;
      for (int loopmonth = 1; loopmonth < NB_MONTH; loopmonth++) {
        cumulativeFactors[loopmonth] = cumulativeFactors[loopmonth - 1] * monthlyFactors[loopmonth - 1];
      }
    }
    //Constructing a 12-periodic vector of the same size of the step vector, and using the vector cumulative.
    final int numberOfSteps = steps.length;
    _monthlyCumulativeFactors = new double[numberOfSteps];
    for (int loopmonth = 0; loopmonth < numberOfSteps; loopmonth++) {
      _monthlyCumulativeFactors[loopmonth] = cumulativeFactors[loopmonth % 12];
    }
  }

  @Override
  public Double evaluate(Double x) {
    Interpolator1DDataBundle dataBundle = INTERPOLATOR_STEP_FLAT.getDataBundleFromSortedArrays(_steps, _monthlyCumulativeFactors);
    return INTERPOLATOR_STEP_FLAT.interpolate(dataBundle, x);
  }
  //------------------------- AUTOGENERATED START -------------------------
  ///CLOVER:OFF
  /**
   * The meta-bean for {@code SeasonalFunction}.
   * @return the meta-bean, not null
   */
  public static SeasonalMonthlyFunction.Meta meta() {
    return SeasonalMonthlyFunction.Meta.INSTANCE;
  }

  static {
    JodaBeanUtils.registerMetaBean(SeasonalMonthlyFunction.Meta.INSTANCE);
  }

  @Override
  public SeasonalMonthlyFunction.Meta metaBean() {
    return SeasonalMonthlyFunction.Meta.INSTANCE;
  }

  @Override
  public <R> Property<R> property(String propertyName) {
    return metaBean().<R>metaProperty(propertyName).createProperty(this);
  }

  @Override
  public Set<String> propertyNames() {
    return metaBean().metaPropertyMap().keySet();
  }

  //-----------------------------------------------------------------------
  /**
   * Gets the cumulative multiplicative seasonal factors from the reference time to the next. Array of size 12 (the 1st is 1.0, it is added to simplify the implementation).
   * @return the value of the property
   */
  private double[] getMonthlyCumulativeFactors() {
    return (_monthlyCumulativeFactors != null ? _monthlyCumulativeFactors.clone() : null);
  }

  /**
   * Gets the the {@code monthlyCumulativeFactors} property.
   * @return the property, not null
   */
  public final Property<double[]> monthlyCumulativeFactors() {
    return metaBean().monthlyCumulativeFactors().createProperty(this);
  }

  //-----------------------------------------------------------------------
  /**
   * Gets the cumulative multiplicative seasonal factors from the reference time to the next. Array of size 12 (the 1st is 1.0, it is added to simplify the implementation).
   * @return the value of the property
   */
  private double[] getSteps() {
    return (_steps != null ? _steps.clone() : null);
  }

  /**
   * Gets the the {@code steps} property.
   * @return the property, not null
   */
  public final Property<double[]> steps() {
    return metaBean().steps().createProperty(this);
  }

  //-----------------------------------------------------------------------
  @Override
  public SeasonalMonthlyFunction clone() {
    BeanBuilder<? extends SeasonalMonthlyFunction> builder = metaBean().builder();
    for (MetaProperty<?> mp : metaBean().metaPropertyIterable()) {
      if (mp.style().isBuildable()) {
        Object value = mp.get(this);
        if (value instanceof Bean) {
          value = JodaBeanUtils.clone((Bean) value);
        }
        builder.set(mp.name(), value);
      }
    }
    return builder.build();
  }

  @Override
  public boolean equals(Object obj) {
    if (obj == this) {
      return true;
    }
    if (obj != null && obj.getClass() == this.getClass()) {
      SeasonalMonthlyFunction other = (SeasonalMonthlyFunction) obj;
      return JodaBeanUtils.equal(getMonthlyCumulativeFactors(), other.getMonthlyCumulativeFactors()) &&
          JodaBeanUtils.equal(getSteps(), other.getSteps());
    }
    return false;
  }

  @Override
  public int hashCode() {
    int hash = getClass().hashCode();
    hash += hash * 31 + JodaBeanUtils.hashCode(getMonthlyCumulativeFactors());
    hash += hash * 31 + JodaBeanUtils.hashCode(getSteps());
    return hash;
  }

  @Override
  public String toString() {
    StringBuilder buf = new StringBuilder(96);
    buf.append("SeasonalFunction{");
    int len = buf.length();
    toString(buf);
    if (buf.length() > len) {
      buf.setLength(buf.length() - 2);
    }
    buf.append('}');
    return buf.toString();
  }

  protected void toString(StringBuilder buf) {
    buf.append("monthlyCumulativeFactors").append('=').append(getMonthlyCumulativeFactors()).append(',').append(' ');
    buf.append("steps").append('=').append(getSteps()).append(',').append(' ');
  }

  //-----------------------------------------------------------------------
  /**
   * The meta-bean for {@code SeasonalFunction}.
   */
  public static class Meta extends DirectMetaBean {
    /**
     * The singleton instance of the meta-bean.
     */
    static final Meta INSTANCE = new Meta();

    /**
     * The meta-property for the {@code monthlyCumulativeFactors} property.
     */
    private final MetaProperty<double[]> _monthlyCumulativeFactors = DirectMetaProperty.ofReadOnly(
        this, "monthlyCumulativeFactors", SeasonalMonthlyFunction.class, double[].class);
    /**
     * The meta-property for the {@code steps} property.
     */
    private final MetaProperty<double[]> _steps = DirectMetaProperty.ofReadOnly(
        this, "steps", SeasonalMonthlyFunction.class, double[].class);
    /**
     * The meta-properties.
     */
    private final Map<String, MetaProperty<?>> _metaPropertyMap$ = new DirectMetaPropertyMap(
        this, null,
        "monthlyCumulativeFactors",
        "steps");

    /**
     * Restricted constructor.
     */
    protected Meta() {
    }

    @Override
    protected MetaProperty<?> metaPropertyGet(String propertyName) {
      switch (propertyName.hashCode()) {
        case 457851908:  // monthlyCumulativeFactors
          return _monthlyCumulativeFactors;
        case 109761319:  // steps
          return _steps;
      }
      return super.metaPropertyGet(propertyName);
    }

    @Override
    public BeanBuilder<? extends SeasonalMonthlyFunction> builder() {
      throw new UnsupportedOperationException();
    }

    @Override
    public Class<? extends SeasonalMonthlyFunction> beanType() {
      return SeasonalMonthlyFunction.class;
    }

    @Override
    public Map<String, MetaProperty<?>> metaPropertyMap() {
      return _metaPropertyMap$;
    }

    //-----------------------------------------------------------------------
    /**
     * The meta-property for the {@code monthlyCumulativeFactors} property.
     * @return the meta-property, not null
     */
    public final MetaProperty<double[]> monthlyCumulativeFactors() {
      return _monthlyCumulativeFactors;
    }

    /**
     * The meta-property for the {@code steps} property.
     * @return the meta-property, not null
     */
    public final MetaProperty<double[]> steps() {
      return _steps;
    }

    //-----------------------------------------------------------------------
    @Override
    protected Object propertyGet(Bean bean, String propertyName, boolean quiet) {
      switch (propertyName.hashCode()) {
        case 457851908:  // monthlyCumulativeFactors
          return ((SeasonalMonthlyFunction) bean).getMonthlyCumulativeFactors();
        case 109761319:  // steps
          return ((SeasonalMonthlyFunction) bean).getSteps();
      }
      return super.propertyGet(bean, propertyName, quiet);
    }

    @Override
    protected void propertySet(Bean bean, String propertyName, Object newValue, boolean quiet) {
      switch (propertyName.hashCode()) {
        case 457851908:  // monthlyCumulativeFactors
          if (quiet) {
            return;
          }
          throw new UnsupportedOperationException("Property cannot be written: monthlyCumulativeFactors");
        case 109761319:  // steps
          if (quiet) {
            return;
          }
          throw new UnsupportedOperationException("Property cannot be written: steps");
      }
      super.propertySet(bean, propertyName, newValue, quiet);
    }

  }

  ///CLOVER:ON
  //-------------------------- AUTOGENERATED END --------------------------
}
