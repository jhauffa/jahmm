package be.ac.ulg.montefiore.run.jahmm.test;

import java.util.ArrayList;
import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussian;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussianFactory;
import junit.framework.TestCase;

public class BuilderTest extends TestCase {

	private static double observationData[] = { 1, 3, 478, 2, 379, 2, 448 };
	private static double cluster1Mean = 2.0;
	private static double cluster2Mean = 435.0;

	public void testOpdfClustering()
	{
		List<List<ObservationReal>> sequences =
				new ArrayList<List<ObservationReal>>(1);
		List<ObservationReal> sequence =
				new ArrayList<ObservationReal>(observationData.length);
		for (double v : observationData)
			sequence.add(new ObservationReal(v));
		sequences.add(sequence);

		Hmm<ObservationReal> hmm = Hmm.<ObservationReal>builder(2)
				.withUniformPi().withUniformA()
				.withOpdfClustering(new OpdfGaussianFactory(), sequences)
				.done();

		double m1 = ((OpdfGaussian) hmm.getOpdf(0)).mean();
		double m2 = ((OpdfGaussian) hmm.getOpdf(1)).mean();
		if (m1 > m2) {
			double tmp = m2;
			m2 = m1;
			m1 = tmp;
		}
		assertEquals(cluster1Mean, m1, 0.1);
		assertEquals(cluster2Mean, m2, 0.1);
	}

}
