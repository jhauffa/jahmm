package be.ac.ulg.montefiore.run.jahmm;

import java.util.List;

public class ViterbiCalculatorLogSpace extends ViterbiCalculator {

	public <O extends Observation>
	ViterbiCalculatorLogSpace(List<? extends O> oseq, HmmLogSpace<O> hmm)
	{
		super(oseq, hmm);
	}


	@Override
	protected <O extends Observation> void computeFirstStep(Hmm<O> hmm, O o0)
	{
		for (int i = 0; i < hmm.nbStates(); i++) {
			delta[0][i] = -hmm.getPi(i) - hmm.getOpdf(i).logProbability(o0);
			psy[0][i] = 0;
		}
	}


	@Override
	protected <O extends Observation> void
	computeStep(Hmm<O> hmm, O o, int t, int j)
	{
		double minDelta = Double.MAX_VALUE;
		int min_psy = 0;

		for (int i = 0; i < hmm.nbStates(); i++) {
			double thisDelta = delta[t-1][i] - hmm.getAij(i, j);

			if (minDelta > thisDelta) {
				minDelta = thisDelta;
				min_psy = i;
			}
		}

		delta[t][j] = minDelta - hmm.getOpdf(j).logProbability(o);
		psy[t][j] = min_psy;
	}

}
