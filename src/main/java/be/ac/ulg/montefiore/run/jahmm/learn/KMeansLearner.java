/*
 * Copyright (c) 2004-2009, Jean-Marc François. All Rights Reserved.
 * Licensed under the New BSD license.  See the LICENSE file.
 */

package be.ac.ulg.montefiore.run.jahmm.learn;

import java.util.*;

import be.ac.ulg.montefiore.run.jahmm.*;


/**
 * An implementation of the K-Means learning algorithm.
 */
public class KMeansLearner<O extends Observation & CentroidFactory<? super O>>
{	
	private Clusters<O> clusters;
	private int nbStates;
	private List<? extends List<? extends O>> obsSeqs;
	private OpdfFactory<? extends Opdf<O>> opdfFactory;
	private boolean terminated;
	
	
	/**
	 * Initializes a K-Means algorithm implementation.  This algorithm
	 * finds a HMM that models a set of observation sequences.
	 *
	 * @param nbStates  The number of states the resulting HMM will be made of.
	 * @param opdfFactory A class that builds the observation probability
	 *                    distributions associated to the states of the HMM.
	 * @param sequences A vector of observation sequences.  Each observation
	 *                sequences is a vector of
	 *                {@link be.ac.ulg.montefiore.run.jahmm.Observation
	 *                observations} compatible with the 
	 *                {@link be.ac.ulg.montefiore.run.jahmm.CentroidFactory
	 *                k-means algorithm}.
	 */
	public KMeansLearner(int nbStates,
			OpdfFactory<? extends Opdf<O>> opdfFactory,
			List<? extends List<? extends O>> sequences)
	{	
		this.obsSeqs = sequences;
		this.opdfFactory = opdfFactory;
		this.nbStates = nbStates;
		
		List<? extends O> observations = Observation.flat(sequences);
		clusters = new Clusters<O>(nbStates, observations);
		terminated = false;
	}
	
	
	/**
	 * Performs one iteration of the K-Means algorithm.
	 * In one iteration, a new HMM is computed using the current clusters, and
	 * the clusters are re-estimated using this HMM.
	 *
	 * @return A new, updated HMM.
	 */
	public Hmm<O> iterate()
	{	
		Hmm<O> hmm = new Hmm<O>(nbStates, opdfFactory);
		
		learnPi(hmm);
		learnAij(hmm);
		learnOpdf(hmm);
		
		terminated = optimizeCluster(hmm);
		
		return hmm;
	}
	
	
	/**
	 * Returns <code>true</code> if the algorithm has reached a fix point,
	 * else returns <code>false</code>.
	 */
	public boolean isTerminated()
	{
		return terminated;
	}
	
	
	/**
	 * Does iterations of the K-Means algorithm until a fix point is reached.
	 * 
	 * @return The HMM that best matches the set of observation sequences given
	 *         (according to the K-Means algorithm).
	 */
	public Hmm<O> learn()
	{	
		Hmm<O> hmm;
		
		do 
			hmm = iterate();
		while(!isTerminated());
		
		return hmm;
	}
	
	
	private void learnPi(Hmm<?> hmm)
	{	
		double[] pi = new double[nbStates];
		Arrays.fill(pi, 0.);
		
		for (List<? extends O> sequence : obsSeqs)
			pi[clusters.clusterNb(sequence.get(0))]++;
		
		for (int i = 0; i < nbStates; i++)
			hmm.setPi(i, pi[i] / obsSeqs.size());
	}
	
	
	private void learnAij(Hmm<O> hmm)
	{	
		for (int i = 0; i < hmm.nbStates(); i++)
			for (int j = 0; j < hmm.nbStates(); j++)
				hmm.setAij(i, j, 0.);
		
		for (List<? extends O> obsSeq : obsSeqs) {
			if (obsSeq.size() < 2)
				continue;
			
			int first_state;
			int second_state = clusters.clusterNb(obsSeq.get(0));
			for (int i = 1; i < obsSeq.size(); i++) {
				first_state = second_state;
				second_state =
					clusters.clusterNb(obsSeq.get(i));
				
				hmm.setAij(first_state, second_state,
						hmm.getAij(first_state, second_state)+1.);
			}
		}
		
		/* Normalize Aij array */
		for (int i = 0; i < hmm.nbStates(); i++) {
			double sum = 0;
			
			for (int j = 0; j < hmm.nbStates(); j++)
				sum += hmm.getAij(i, j);
			
			if (sum == 0.)
				for (int j = 0; j < hmm.nbStates(); j++) 
					hmm.setAij(i, j, 1. / hmm.nbStates());     // Arbitrarily
			else
				for (int j = 0; j < hmm.nbStates(); j++)
					hmm.setAij(i, j, hmm.getAij(i, j) / sum);
		}
	}
	
	
	private void learnOpdf(Hmm<O> hmm)
	{
		for (int i = 0; i < hmm.nbStates(); i++) {
			Collection<O> clusterObservations = clusters.cluster(i);
			
			if (clusterObservations.isEmpty())
				hmm.setOpdf(i, opdfFactory.factor());
			else
				hmm.getOpdf(i).fit(clusterObservations);
		}
	}
	
	
	/* Return true if no modification */
	private boolean optimizeCluster(Hmm<O> hmm)
	{	
		boolean modif = false;
		
		for (List<? extends O> obsSeq : obsSeqs) {
			ViterbiCalculator vc = new ViterbiCalculator(obsSeq, hmm);
			int states[] = vc.stateSequence();
			
			for (int i = 0; i < states.length; i++) {
				O o = obsSeq.get(i);
				
				int curClusterNb = clusters.clusterNb(o);
				if (curClusterNb != states[i]) {
					modif = true;
					clusters.move(o, curClusterNb, states[i]);
				}
			}
		}
		
		return !modif;
	}
}


/*
 * This class holds the matching between observations and clusters.
 */
class Clusters<O extends CentroidFactory<? super O>>
{	
	private HashMap<O,Integer> clustersHash;
	private ArrayList<Collection<O>> clusters;
	
	
	public Clusters(int k, List<? extends O> observations)
	{
		clustersHash = new HashMap<O,Integer>();
		clusters = new ArrayList<Collection<O>>(k);
		
		KMeansCalculator<O> kmc = new KMeansCalculator<O>(k, observations);
		
		for (int i = 0; i < k; i++) {
			Collection<O> cluster = new HashSet<O>(kmc.cluster(i));
			clusters.add(cluster);
			
			for (O element : cluster) 
				clustersHash.put(element, i);
		}
	}
	
	
	public int clusterNb(O o)
	{
		return clustersHash.get(o);
	}
	
	
	public Collection<O> cluster(int clusterNb)
	{
		return clusters.get(clusterNb);
	}
	
	
	public void move(O o, int fromClusterNb, int toClusterNb)
	{
		clusters.get(fromClusterNb).remove(o);
		clusters.get(toClusterNb).add(o);
		clustersHash.put(o, toClusterNb);
	}
}
