/*
 * Copyright (c) 2004-2009, Jean-Marc François. All Rights Reserved.
 * Licensed under the New BSD license.  See the LICENSE file.
 */

package be.ac.ulg.montefiore.run.distributions;

import java.util.Random;


/**
 * This class implements a Gaussian distribution.
 */
public class GaussianDistribution
implements RandomDistribution {
	
	private final double mean;
	private final double deviation;
	private final double variance;
	private final double logVarianceTerm;
	private final double varianceTerm;
	private final static Random randomGenerator = new Random();
	
	final static double PI2, LOG_PI2;
	static {
		PI2 = 2. * Math.PI;
		LOG_PI2 = Math.log(PI2);
	}
	
	/**
	 * Creates a new pseudo-random, Gaussian distribution with zero mean
	 * and unitary variance.
	 */
	public GaussianDistribution()
	{
		this(0., 1.);
	}
	
	
	/**
	 * Creates a new pseudo-random, Gaussian distribution.
	 *
	 * @param mean The mean value of the generated numbers.
	 * @param variance The variance of the generated numbers.
	 */
	public GaussianDistribution(double mean, double variance)
	{
		if (variance <= 0.)
			throw new IllegalArgumentException("Variance must be positive");
		
		this.mean = mean;
		this.variance = variance;
		this.logVarianceTerm = LOG_PI2 + Math.log(variance);
		this.varianceTerm = Math.pow(PI2 * variance, -.5);
		this.deviation = Math.sqrt(variance);
	}
	
	
	/**
	 * Returns this distribution's mean value.
	 *
	 * @return This distribution's mean value.
	 */
	public double mean()
	{
		return mean;
	}
	
	
	/**
	 * Returns this distribution's variance.
	 *
	 * @return This distribution's variance.
	 */
	public double variance()
	{
		return variance;
	}
	
	
	public double generate()
	{
		return randomGenerator.nextGaussian() * deviation + mean;
	}
	
	
	public double logProbability(double n)
	{
		double expArg = (n - mean) * (n - mean) / variance;
		return -.5 * (logVarianceTerm + expArg);
	}
	
	
	public double probability(double n)
	{
		double expArg = -.5 * (n - mean) * (n - mean) / variance;
		return varianceTerm * Math.exp(expArg);
	}
	
	
	private static final long serialVersionUID = 8037555496366734729L;
}
