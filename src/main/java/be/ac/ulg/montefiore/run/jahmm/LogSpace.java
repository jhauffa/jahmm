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

	public static double sum(double x, double y) {
		if (Double.isNaN(x))
			return y;
		if (Double.isNaN(y))
			return x;

		if (x > y)
			return x + Math.log(1.0 + Math.exp(y - x));
		return y + Math.log(1.0 + Math.exp(x - y));
	}

	public static double product(double x, double y) {
		return x + y;
	}

	public static double quotient(double x, double y) {
		return x - y;
	}

}
