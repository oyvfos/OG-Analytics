package com.opengamma.analytics.financial.model.finitedifference;

import com.opengamma.analytics.math.function.Function;

public class SandBox {
	 final static double thetaHW = 0.000023637253287033743; 
	  final static double kappa = 0.008622227075117428; 
	  final static double r0 = -0.00024031738413887727;
	  final static double sigma = 0.000212693; 
	final static Function<Double, Double> vas = new Function<Double, Double>() {
        @Override
        public Double evaluate(final Double... ts) {
          final double t = ts[0];
          final double r0 = ts[1];
          return (r0 - thetaHW/kappa)/Math.exp(kappa*t) + thetaHW/kappa; 
          
        }
      };
      final static Function<Double, Double> B = new Function<Double, Double>() {
          @Override
          public Double evaluate(final Double... ts) { 
            final double t = ts[0];
            final double T = ts[1];
            return (1 - Math.exp((-kappa)*(T - t)))/kappa; 
          }
        };
        final static Function<Double, Double> A = new Function<Double, Double>() {
            @Override
            public Double evaluate(final Double... ts) {
              final double t = ts[0];
              final double T = ts[1];
              return (thetaHW/kappa - Math.pow(sigma,2)/(2*Math.pow(kappa, 2)))*(B.evaluate(t,T) - (T - t)) - (Math.pow(sigma, 2)*Math.pow(B.evaluate(t, T),2)/(4*kappa));               
            }
        };
        final static Function<Double, Double> Z = new Function<Double, Double>() {
            @Override
            public Double evaluate(final Double... ts) {
              final double t = ts[0];
              final double T = ts[1];
              final double r = ts[2];
              return Math.exp(A.evaluate(t, T) - vas.evaluate(t, r)* B.evaluate(t, T));
            }
          };
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		System.out.println(Z.evaluate(new Double[] {0d,110d}));
	}
	public static double price(Double t, Double T,Double r0) {
		// TODO Auto-generated method stub
		
		return Z.evaluate(new Double[] {t,T,r0});
	}
}
