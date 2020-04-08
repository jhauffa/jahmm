package be.ac.ulg.montefiore.run.jahmm;

import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Collection;

import be.ac.ulg.montefiore.run.distributions.LogNormalDistribution;


/**
 * This class represents a (monovariate) log-normal distribution function.
 */
public class OpdfLogNormal
implements Opdf<ObservationReal>
{
	private LogNormalDistribution distribution;


	/**
	 * Builds a new log-normal probability distribution with zero mean and
	 * unit variance.
	 */
	public OpdfLogNormal()
	{
		distribution = new LogNormalDistribution();
	}


	/**
	 * Builds a new log-normal probability distribution with a given mean and
	 * variance.
	 *
	 * @param mean The distribution's mean.
	 * @param variance The distribution's variance.
	 */
	public OpdfLogNormal(double mean, double variance)
	{
		distribution = new LogNormalDistribution(mean, variance);
	}


	/**
	 * Returns this distribution's mean value.
	 *
	 * @return This distribution's mean value.
	 */
	public double mean()
	{
		return distribution.mean();
	}


	/**
	 * Returns this distribution's variance.
	 *
	 * @return This distribution's variance.
	 */
	public double variance()
	{
		return distribution.variance();
	}


	public double arithmeticMean()
	{
		return distribution.arithmeticMean();
	}


	public double arithmeticVariance()
	{
		return distribution.arithmeticVariance();
	}


	public double logProbability(ObservationReal o)
	{
		return distribution.logProbability(o.value);
	}


	public double probability(ObservationReal o)
	{
		return distribution.probability(o.value);
	}


	public ObservationReal generate()
	{
		return new ObservationReal(distribution.generate());
	}


	public void fit(ObservationReal... oa)
	{
		fit(Arrays.asList(oa));
	}


	public void fit(Collection<? extends ObservationReal> co)
	{
		double[] weights = new double[co.size()];
		Arrays.fill(weights, 1. / co.size());

		fit(co, weights);
	}


	public void fit(ObservationReal[] o, double[] weights)
	{
		fit(Arrays.asList(o), weights);
	}


	public void fit(Collection<? extends ObservationReal> co,
			double[] weights)
	{
		if (co.isEmpty() || co.size() != weights.length)
			throw new IllegalArgumentException();

		// Compute mean
		double mean = 0.;
		int i = 0;
		for (ObservationReal o : co) {
			if (o.value <= 0.0) {
				throw new IllegalArgumentException("value " + o.value +
						" out of domain");
			}
			mean += Math.log(o.value) * weights[i++];
		}

		// Compute variance
		double variance = 0.;
		i = 0;
		for (ObservationReal o : co) {
			double d = Math.log(o.value) - mean;

			variance += d * d * weights[i++];
		}

		distribution = new LogNormalDistribution(mean, variance);
	}


	public OpdfLogNormal clone()
	{
		try {
			return (OpdfLogNormal) super.clone();
		} catch(CloneNotSupportedException e) {
            throw new AssertionError(e);
        }
	}


	public String toString()
	{
		return toString(NumberFormat.getInstance());
	}


	public String toString(NumberFormat numberFormat)
	{
		return "Log-normal distribution --- " +
		"Mean: " + numberFormat.format(distribution.mean()) +
		" Variance " + numberFormat.format(distribution.variance());
	}


	private static final long serialVersionUID = 1L;
}
