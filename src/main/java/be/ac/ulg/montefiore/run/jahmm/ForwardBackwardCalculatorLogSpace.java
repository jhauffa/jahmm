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
	private double[] s;

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
	compute(List<? extends O> oseq, Hmm<O> hmm, EnumSet<Computation> flags)
	{
		s = new double[hmm.nbStates()];
		super.compute(oseq, hmm, flags);
	}


	@Override
	protected <O extends Observation> void
	computeAlphaInit(Hmm<? super O> hmm, O o, int i)
	{
		alpha[0][i] = LogSpace.product(hmm.getPi(i),
				hmm.getOpdf(i).logProbability(o));
	}


	@Override
	protected <O extends Observation> void
	computeAlphaStep(Hmm<? super O> hmm, O o, int t, int j)
	{
		for (int i = 0; i < s.length; i++)
			s[i] = LogSpace.product(alpha[t-1][i], hmm.getAij(i, j));

		alpha[t][j] = LogSpace.product(LogSpace.INST.sum(s),
				hmm.getOpdf(j).logProbability(o));
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
		for (int j = 0; j < s.length; j++) {
			s[j] = LogSpace.product(beta[t+1][j],
					LogSpace.product(hmm.getAij(i, j),
							hmm.getOpdf(j).logProbability(o)));
		}

		beta[t][i] = LogSpace.INST.sum(s);
	}


	@Override
	protected <O extends Observation> void
	computeProbability(List<O> oseq, Hmm<? super O> hmm,
			EnumSet<Computation> flags)
	{
		double[] p;
		if (flags.contains(Computation.ALPHA)) {
			p = alpha[oseq.size()-1];
		} else {
			p = s;
			for (int i = 0; i < hmm.nbStates(); i++) {
				p[i] = LogSpace.product(hmm.getPi(i),
						LogSpace.product(
								hmm.getOpdf(i).logProbability(oseq.get(0)),
								beta[0][i]));
			}
		}

		lnProbability = LogSpace.INST.sum(p);
		probability = LogSpace.INST.exp(lnProbability);
	}


	public double lnProbability()
	{
		return lnProbability;
	}
}
