package be.ac.ulg.montefiore.run.jahmm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

public class HmmBuilder<O extends Observation> {

	protected final int nbStates;

	protected double pi[];
	protected double a[][];
	protected List<Opdf<O>> opdfs;

	public HmmBuilder(int nbStates)
	{
		this.nbStates = nbStates;
		if (nbStates <= 0) {
			throw new IllegalArgumentException("Number of states must be " +
					"positive");
		}
	}

	public HmmBuilder<O> withPi(double[] pi)
	{
		this.pi = pi;
		return this;
	}

	public HmmBuilder<O> withUniformPi()
	{
		pi = new double[nbStates];
		Arrays.fill(pi, 1.0 / nbStates);
		return this;
	}

	public HmmBuilder<O> withRandomPi()
	{
		pi = new double[nbStates];
		generateRandomUnitVector(pi);
		return this;
	}

	public HmmBuilder<O> withA(double[][] a)
	{
		this.a = a;
		return this;
	}

	public HmmBuilder<O> withUniformA()
	{
		a = new double[nbStates][nbStates];
		for (int i = 0; i < nbStates; i++)
			Arrays.fill(a[i], 1.0 / nbStates);
		return this;
	}

	public HmmBuilder<O> withRandomA()
	{
		a = new double[nbStates][nbStates];
		for (int i = 0; i < nbStates; i++)
			generateRandomUnitVector(a[i]);
		return this;
	}

	public HmmBuilder<O> withOpdfs(List<? extends Opdf<O>> opdfs)
	{
		this.opdfs = new ArrayList<Opdf<O>>(opdfs);
		return this;
	}

	public HmmBuilder<O> withOpdfFactory(
			OpdfFactory<? extends Opdf<O>> opdfFactory)
	{
		opdfs = new ArrayList<Opdf<O>>(nbStates);
		for (int i = 0; i < nbStates; i++)
			opdfs.add(opdfFactory.factor());
		return this;
	}

	@SuppressWarnings("unchecked")	// have to assume that T can be cast to O
	public <T extends Observation & CentroidFactory<? super T>>
	HmmBuilder<O> withOpdfClustering(
			OpdfFactory<? extends Opdf<T>> opdfFactory,
			List<? extends List<? extends T>> sequences)
	{
		opdfs = new ArrayList<Opdf<O>>(nbStates);
		List<? extends T> observations = Observation.flat(sequences);
		KMeansCalculator<T> kmc = new KMeansCalculator<T>(nbStates,
				observations);
		for (int i = 0; i < nbStates; i++) {
			Collection<T> clusterObservations = kmc.cluster(i);

			Opdf<T> opdf = opdfFactory.factor();
			if (!clusterObservations.isEmpty())
				opdf.fit(clusterObservations);
			opdfs.add((Opdf<O>) opdf);
		}
		return this;
	}

	public Hmm<O> done()
	{
		if (pi == null)
			pi = new double[nbStates];
		if (a == null)
			a = new double[nbStates][nbStates];
		if (opdfs == null) {
			opdfs = new ArrayList<Opdf<O>>(nbStates);
			for (int i = 0; i < nbStates; i++)
				opdfs.add(null);
		}
		return new Hmm<O>(pi, a, opdfs);
	}

	private static void generateRandomUnitVector(double[] v)
	{
		double sum = 0.0;
		for (int i = 0; i < v.length; i++) {
			v[i] = Math.random();
			sum += v[i];
		}
		for (int i = 0; i < v.length; i++)
			v[i] /= sum;
	}

}
