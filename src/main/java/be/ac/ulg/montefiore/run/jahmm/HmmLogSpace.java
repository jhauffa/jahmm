package be.ac.ulg.montefiore.run.jahmm;

import java.util.Iterator;
import java.util.List;

public class HmmLogSpace<O extends Observation> extends Hmm<O> {

	protected HmmLogSpace(int nbStates)
	{
		super(nbStates);
	}

	public HmmLogSpace(Hmm<O> hmm)
	{
		super(hmm.nbStates());
		for (int i = 0; i < nbStates(); i++) {
			setPi(i, LogSpace.log(hmm.getPi(i)));
			for (int j = 0; j < nbStates(); j++)
				setAij(i, j, LogSpace.log(hmm.getAij(i, j)));
			setOpdf(i, hmm.getOpdf(i).clone());
		}
	}

	public Hmm<O> toHmm()
	{
		Hmm<O> hmm = new Hmm<O>(nbStates());
		for (int i = 0; i < nbStates(); i++) {
			hmm.setPi(i, LogSpace.exp(getPi(i)));
			for (int j = 0; j < nbStates(); j++)
				hmm.setAij(i, j, LogSpace.exp(getAij(i, j)));
			hmm.setOpdf(i, getOpdf(i).clone());
		}
		return hmm;
	}

	@Override
	public HmmLogSpace<O> clone() throws CloneNotSupportedException
	{
		HmmLogSpace<O> hmm = new HmmLogSpace<O>(nbStates());
		for (int i = 0; i < nbStates(); i++) {
			hmm.setPi(i, getPi(i));
			for (int j = 0; j < nbStates(); j++)
				hmm.setAij(i, j, getAij(i, j));
			hmm.setOpdf(i, getOpdf(i).clone());
		}
		return hmm;
	}

	@Override
	public int[] mostLikelyStateSequence(List<? extends O> oseq)
	{
		return (new ViterbiCalculatorLogSpace(oseq, this)).stateSequence();
	}

	@Override
	public double probability(List<? extends O> oseq)
	{
		return (new ForwardBackwardCalculatorLogSpace(oseq, this)).
				probability();
	}

	@Override
	public double lnProbability(List<? extends O> oseq)
	{
		return (new ForwardBackwardCalculatorLogSpace(oseq, this)).
				lnProbability();
	}

	@Override
	public double probability(List<? extends O> oseq, int[] sseq)
	{
		return LogSpace.exp(probability(oseq, sseq));
	}

	public double lnProbability(List<? extends O> oseq, int[] sseq)
	{
		if (oseq.size() != sseq.length || oseq.isEmpty())
			throw new IllegalArgumentException();

		double probability = getPi(sseq[0]);

		Iterator<? extends O> oseqIterator = oseq.iterator();

		for (int i = 0; i < sseq.length-1; i++) {
			probability = LogSpace.product(probability,
				LogSpace.product(
					LogSpace.log(getOpdf(sseq[i]).
							probability(oseqIterator.next())),
					getAij(sseq[i], sseq[i+1])));
		}

		return LogSpace.product(probability,
				LogSpace.log(getOpdf(sseq[sseq.length-1]).
						probability(oseq.get(sseq.length-1))));
	}

	private static final long serialVersionUID = 1L;

}
