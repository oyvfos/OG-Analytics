package com.opengamma.analytics.financial.model.finitedifference.applications;

import com.opengamma.analytics.math.function.Function1D;
import com.opengamma.analytics.math.interpolation.DoubleQuadraticInterpolator1D;
import com.opengamma.analytics.math.interpolation.data.ArrayInterpolator1DDataBundle;
import com.opengamma.analytics.math.interpolation.data.Interpolator1DDoubleQuadraticDataBundle;
import com.opengamma.util.csv.CSVDocumentReader;

public class FunctionSandboxTest {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		double[] yValues = new double[] {1., 2., 3. };

	    double[] xValues = new double[] {1., 2., 3. };
		ArrayInterpolator1DDataBundle bundle = new ArrayInterpolator1DDataBundle(xValues, yValues);
		DoubleQuadraticInterpolator1D interp = new DoubleQuadraticInterpolator1D();
		Interpolator1DDoubleQuadraticDataBundle bundle0 = interp.getDataBundleFromSortedArrays(xValues,yValues);
		Function1D<Double, Double> fun = interp.getFunction(bundle0);
		fun.evaluate(1.5);
		CSVDocumentReader _csvDocReader = new CSVDocumentReader(FunctionSandboxTest.class.getResource("/src/test/resources/data/d.csv"));
	}

}
