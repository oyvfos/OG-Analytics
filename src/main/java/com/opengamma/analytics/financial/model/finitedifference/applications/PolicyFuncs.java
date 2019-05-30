package com.opengamma.analytics.financial.model.finitedifference.applications;

import java.util.ArrayList;
import java.util.List;

import com.opengamma.analytics.math.function.Function1D;
import com.opengamma.analytics.math.interpolation.DoubleQuadraticInterpolator1D;
import com.opengamma.analytics.math.interpolation.LinearInterpolator1D;
import com.opengamma.analytics.math.interpolation.data.Interpolator1DDataBundle;
import com.opengamma.analytics.math.interpolation.data.Interpolator1DDoubleQuadraticDataBundle;
import com.opengamma.strata.collect.io.CsvFile;
import com.opengamma.strata.collect.io.CsvRow;
import com.opengamma.strata.collect.io.ResourceLocator;

public class PolicyFuncs {
	public static void main(String[] args) {
		funcs();
	}
	
	public static List<Function1D<Double, Double>> funcs() {
		// TODO Auto-generated method stub
		List<Function1D<Double, Double>> myList = new ArrayList<>();
		//DoubleQuadraticInterpolator1D interp = new DoubleQuadraticInterpolator1D();
		LinearInterpolator1D interp = new LinearInterpolator1D();
		ResourceLocator re = ResourceLocator.of("classpath:test/d.csv");
		CsvFile csv = CsvFile.of(re.getCharSource(), true);
		double[] prem = new double[csv.rows().size()]; 
		double[] qx = new double[csv.rows().size()];;
		double[] seq = new double[csv.rows().size()];;
		for (CsvRow row : csv.rows()) {
		      String str = row.getField("qx");
		      String str1 = row.getField("prem");
		      prem[row.lineNumber()-2] = str.isEmpty() ? 0d : Double.parseDouble(str1);
		      qx[row.lineNumber()-2]=str1.isEmpty() ? 0d : Double.parseDouble(str);
		      seq[row.lineNumber()-2]=row.lineNumber()-2;
		}
		 Interpolator1DDataBundle bundleqx = interp.getDataBundleFromSortedArrays(seq,qx);
		 Interpolator1DDataBundle bundleprem = interp.getDataBundleFromSortedArrays(seq,prem);
		//final ImmutableList<Function1D<Double, Double>> funcs;
		myList.add(interp.getFunction(bundleqx));
		myList.add(interp.getFunction(bundleprem));
		return myList;
		//CSVDocumentReader _csvDocReader = new CSVDocumentReader(FunctionSandboxTest.class.getResource("/test/d.csv"));
		
	}

}
