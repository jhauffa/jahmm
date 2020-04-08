package be.ac.ulg.montefiore.run.jahmm;

/**
 * Log-space arithmetic as described in
 * Mann, 2006, Numerically Stable Hidden Markov Model Implementation.
 * http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf
 */
public class LogSpace {

	public static final Double ZERO = Double.NaN;
	public static final Double ONE = 0.0;

	public static double exp(double x) {
		if (Double.isNaN(x))
			return 0.0;
		return Math.exp(x);
	}

	public static double log(double x) {
		if (x == 0.0)
			return ZERO;
		return Math.log(x);
	}

	public static double sum(double[] x) {
		double max = Double.NEGATIVE_INFINITY;
		for (double v : x)
			if (v > max)
				max = v;

		double sumExp = 0.0;
		for (double v : x)
			sumExp += exp(v - max);
		return max + log(sumExp);
	}

	public static double sum(double[][] x) {
		double max = Double.NEGATIVE_INFINITY;
		for (double[] r : x)
			for (double v : r)
				if (v > max)
					max = v;

		double sumExp = 0.0;
		for (double[] r : x)
			for (double v : r)
				sumExp += exp(v - max);
		return max + log(sumExp);
	}

	public static double product(double x, double y) {
		return x + y;
	}

	public static double quotient(double x, double y) {
		return x - y;
	}

}
