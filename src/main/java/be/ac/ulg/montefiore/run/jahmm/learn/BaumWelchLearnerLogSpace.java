/*
 * Copyright (c) 2004-2009, Jean-Marc François. All Rights Reserved.
 * Licensed under the New BSD license.  See the LICENSE file.
 */

package be.ac.ulg.montefiore.run.jahmm.learn;

import java.util.Arrays;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.ForwardBackwardCalculatorLogSpace;
import be.ac.ulg.montefiore.run.jahmm.HmmLogSpace;
import be.ac.ulg.montefiore.run.jahmm.LogSpace;
import be.ac.ulg.montefiore.run.jahmm.Observation;

/**
 * An implementation of the Baum-Welch learning algorithm in log-space as
 * described in
 * Mann, 2006, Numerically Stable Hidden Markov Model Implementation.
 * http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf
 */
public class BaumWelchLearnerLogSpace
{
	private int nbIterations = 9;

	public <O extends Observation> HmmLogSpace<O>
	iterate(HmmLogSpace<O> hmm, List<? extends List<? extends O>> sequences)
	{
		HmmLogSpace<O> nhmm;
		try {
			nhmm = hmm.clone();
		} catch(CloneNotSupportedException e) {
			throw new InternalError();
		}

		double allGamma[][][] = new double[sequences.size()][][];
		double[] s = new double[allGamma.length];
		double piDen = LogSpace.log(allGamma.length);

		/* update transition probabilities Aij */
		double aijDen[] = new double[hmm.nbStates()];
		double aijNum[][] = new double[hmm.nbStates()][hmm.nbStates()];

		Arrays.fill(aijDen, LogSpace.ZERO);
		for (int i = 0; i < hmm.nbStates(); i++)
			Arrays.fill(aijNum[i], LogSpace.ZERO);

		int g = 0;
		for (List<? extends O> obsSeq : sequences) {
			ForwardBackwardCalculatorLogSpace fbc =
					generateForwardBackwardCalculator(obsSeq, nhmm);

			double xi[][][] = estimateXi(obsSeq, fbc, nhmm);
			double gamma[][] = allGamma[g++] = estimateGamma(xi, fbc);

			double[] os = new double[obsSeq.size()];
			for (int i = 0; i < hmm.nbStates(); i++) {
				os[0] = aijDen[i];
				for (int t = 0; t < os.length - 1; t++)
					os[t + 1] = gamma[t][i];
				aijDen[i] = LogSpace.sum(os);

				for (int j = 0; j < hmm.nbStates(); j++) {
					os[0] = aijNum[i][j];
					for (int t = 0; t < os.length - 1; t++)
						os[t + 1] = xi[t][i][j];
					aijNum[i][j] = LogSpace.sum(os);
				}
			}
		}

		for (int i = 0; i < hmm.nbStates(); i++)
			for (int j = 0; j < hmm.nbStates(); j++)
				nhmm.setAij(i, j, LogSpace.quotient(aijNum[i][j], aijDen[i]));

		/* update initial probabilities Pi */
		for (int i = 0; i < hmm.nbStates(); i++) {
			for (int t = 0; t < allGamma.length; t++)
				s[t] = allGamma[t][0][i];
			nhmm.setPi(i, LogSpace.quotient(LogSpace.sum(s), piDen));
		}

		/* update output PDFs */
		List<O> observations = Observation.flat(sequences);
		double[] weights = new double[observations.size()];
		for (int i = 0; i < hmm.nbStates(); i++) {
			int obsIdx = 0;
			for (double[][] gamma : allGamma) {
				for (double[] gammaT : gamma)
					weights[obsIdx++] = gammaT[i];
			}
			double sum = LogSpace.sum(weights);

			for (int j = 0; j < obsIdx; j++)
				weights[j] = LogSpace.exp(LogSpace.quotient(weights[j], sum));

			nhmm.getOpdf(i).fit(observations, weights);
		}

		return nhmm;
	}


	protected <O extends Observation> ForwardBackwardCalculatorLogSpace
	generateForwardBackwardCalculator(List<? extends O> sequence,
			HmmLogSpace<O> hmm)
	{
		return new ForwardBackwardCalculatorLogSpace(sequence, hmm,
				EnumSet.of(ForwardBackwardCalculatorLogSpace.Computation.ALPHA,
						ForwardBackwardCalculatorLogSpace.Computation.BETA));
	}


	public <O extends Observation> HmmLogSpace<O>
	learn(HmmLogSpace<O> initialHmm,
			List<? extends List<? extends O>> sequences)
	{
		HmmLogSpace<O> hmm = initialHmm;
		for (int i = 0; i < nbIterations; i++)
			hmm = iterate(hmm, sequences);
		return hmm;
	}


	protected <O extends Observation> double[][][]
	estimateXi(List<? extends O> sequence,
			ForwardBackwardCalculatorLogSpace fbc, HmmLogSpace<O> hmm)
	{
		if (sequence.size() <= 1) {
			throw new IllegalArgumentException("Observation sequence too " +
					"short");
		}

		double xi[][][] =
				new double[sequence.size()-1][hmm.nbStates()][hmm.nbStates()];

		Iterator<? extends O> seqIterator = sequence.iterator();
		seqIterator.next();

		for (int t = 0; t < sequence.size() - 1; t++) {
			O o = seqIterator.next();

			for (int i = 0; i < hmm.nbStates(); i++) {
				for (int j = 0; j < hmm.nbStates(); j++) {
					xi[t][i][j] = LogSpace.product(fbc.alphaElement(t, i),
							LogSpace.product(hmm.getAij(i, j),
								LogSpace.product(
									hmm.getOpdf(j).logProbability(o),
									fbc.betaElement(t+1, j))));
				}
			}
			double norm = LogSpace.sum(xi[t]);
			for (int i = 0; i < hmm.nbStates(); i++)
				for (int j = 0; j < hmm.nbStates(); j++)
					xi[t][i][j] = LogSpace.quotient(xi[t][i][j], norm);
		}

		return xi;
	}


	/** gamma is directly computed from alpha and beta, xi is not used */
	protected double[][]
	estimateGamma(double[][][] xi, ForwardBackwardCalculatorLogSpace fbc)
	{
		double[][] gamma = new double[xi.length + 1][xi[0].length];

		for (int t = 0; t < xi.length + 1; t++) {
			for (int i = 0; i < xi[0].length; i++) {
				gamma[t][i] = LogSpace.product(fbc.alphaElement(t, i),
						fbc.betaElement(t, i));
			}
			double norm = LogSpace.sum(gamma[t]);
			for (int i = 0; i < xi[0].length; i++)
				gamma[t][i] = LogSpace.quotient(gamma[t][i], norm);
		}

		return gamma;
	}


	public int getNbIterations()
	{
		return nbIterations;
	}


	public void setNbIterations(int nb)
	{
		if (nb < 0)
			throw new IllegalArgumentException("Positive number expected");
		nbIterations = nb;
	}

}
