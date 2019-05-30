package com.opengamma.analytics.financial.model.finitedifference.applications;

import com.opengamma.analytics.math.function.Function1D;

public class PolicyFuncTest {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		// TODO Auto-generated method stub
				Function1D<Double, Double> f1 = PolicyFuncs.funcs().get(0);
				 for (int i = 0; i < 112; i++) {
				        System.out.println(f1.evaluate((double)i));
				      }
	}

}
