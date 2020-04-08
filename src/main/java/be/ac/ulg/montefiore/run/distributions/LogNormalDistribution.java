package be.ac.ulg.montefiore.run.distributions;

import java.util.Random;

/* implementation adapted from Apache Commons Statistics and Commons RNG */
public class LogNormalDistribution
implements RandomDistribution {

	private final double mean;
	private final double variance;
	private final double standardDeviation;
	private final double logTerm;
	private final static double SQRT2PI = Math.sqrt(2. * Math.PI);
	private final static Random randomGenerator = new Random();

	public LogNormalDistribution()
	{
		this(0., 1.);
	}

	public LogNormalDistribution(double mean, double variance)
	{
		if (variance <= 0.)
			throw new IllegalArgumentException("Variance must be positive");

		this.mean = mean;
		this.variance = variance;
		this.standardDeviation = Math.sqrt(variance);
		this.logTerm = Math.log(standardDeviation) +
				0.5 * GaussianDistribution.LOG_PI2;
	}

	public double mean()
	{
		return mean;
	}

	public double variance()
	{
		return variance;
	}

	public double arithmeticMean()
	{
		return Math.exp(mean + 0.5 * variance);
	}

	public double arithmeticVariance()
	{
		return (Math.exp(variance) - 1.0) * Math.exp((2.0 * mean) + variance);
	}

	public double generate()
	{
		return Math.exp(mean + standardDeviation *
				randomGenerator.nextGaussian());
	}

	public double logProbability(double n)
	{
		if (n <= 0.)
			return Double.NaN;
		double logN = Math.log(n);
		double n1 = (logN - mean) / standardDeviation;
		return -0.5 * n1 * n1 - (logTerm + logN);
	}

	public double probability(double n)
	{
		if (n <= 0.)
			return 0.;
		double n1 = (Math.log(n) - mean) / standardDeviation;
		return Math.exp(-0.5 * n1 * n1) / (standardDeviation * SQRT2PI * n);
	}


	private static final long serialVersionUID = -6013042126599196058L;
}
