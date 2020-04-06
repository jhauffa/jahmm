/*
 * Copyright (c) 2004-2009, Jean-Marc François. All Rights Reserved.
 * Licensed under the New BSD license.  See the LICENSE file.
 */

package be.ac.ulg.montefiore.run.jahmm;

import java.util.EnumSet;
import java.util.List;


/**
 * Performs forward/backward calculations in log-space as described in
 * Mann, 2006, Numerically Stable Hidden Markov Model Implementation.
 * http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf
 */
public class ForwardBackwardCalculatorLogSpace extends ForwardBackwardCalculator
{
	private double lnProbability;

	public <O extends Observation>
	ForwardBackwardCalculatorLogSpace(List<? extends O> oseq,
			HmmLogSpace<O> hmm, EnumSet<Computation> flags)
	{
		super(oseq, hmm, flags);
	}


	public <O extends Observation>
	ForwardBackwardCalculatorLogSpace(List<? extends O> oseq,
			HmmLogSpace<O> hmm)
	{
		super(oseq, hmm);
	}


	@Override
	protected <O extends Observation> void
	computeAlphaInit(Hmm<? super O> hmm, O o, int i)
	{
		alpha[0][i] = LogSpace.product(hmm.getPi(i),
				LogSpace.log(hmm.getOpdf(i).probability(o)));
	}


	@Override
	protected <O extends Observation> void
	computeAlphaStep(Hmm<? super O> hmm, O o, int t, int j)
	{
		double sum = LogSpace.ZERO;

		for (int i = 0; i < hmm.nbStates(); i++) {
			sum = LogSpace.sum(sum,
					LogSpace.product(alpha[t-1][i], hmm.getAij(i, j)));
		}

		alpha[t][j] = LogSpace.product(sum,
				LogSpace.log(hmm.getOpdf(j).probability(o)));
	}


	@Override
	protected <O extends Observation> void
	computeBeta(Hmm<? super O> hmm, List<O> oseq)
	{
		beta = new double[oseq.size()][hmm.nbStates()];

		for (int i = 0; i < hmm.nbStates(); i++)
			beta[oseq.size()-1][i] = LogSpace.ONE;

		for (int t = oseq.size()-2; t >= 0; t--)
			for (int i = 0; i < hmm.nbStates(); i++)
				computeBetaStep(hmm, oseq.get(t+1), t, i);
	}


	@Override
	protected <O extends Observation> void
	computeBetaStep(Hmm<? super O> hmm, O o, int t, int i)
	{
		double sum = LogSpace.ZERO;

		for (int j = 0; j < hmm.nbStates(); j++) {
			sum = LogSpace.sum(sum, LogSpace.product(beta[t+1][j],
					LogSpace.product(hmm.getAij(i, j),
							LogSpace.log(hmm.getOpdf(j).probability(o)))));
		}

		beta[t][i] = sum;
	}


	@Override
	protected <O extends Observation> void
	computeProbability(List<O> oseq, Hmm<? super O> hmm,
			EnumSet<Computation> flags)
	{
		lnProbability = LogSpace.ZERO;

		if (flags.contains(Computation.ALPHA)) {
			for (int i = 0; i < hmm.nbStates(); i++) {
				lnProbability = LogSpace.sum(lnProbability,
						alpha[oseq.size()-1][i]);
			}
		} else {
			for (int i = 0; i < hmm.nbStates(); i++) {
				lnProbability = LogSpace.sum(lnProbability,
						LogSpace.product(hmm.getPi(i), LogSpace.product(
								LogSpace.log(
									hmm.getOpdf(i).probability(oseq.get(0))),
								beta[0][i])));
			}
		}

		probability = LogSpace.exp(lnProbability);
	}


	public double lnProbability()
	{
		return lnProbability;
	}
}
