package be.ac.ulg.montefiore.run.jahmm;


/**
 * This class can build <code>OpdfLogNormal</code> observation probability
 * functions.
 */
public class OpdfLogNormalFactory
implements OpdfFactory<OpdfLogNormal>
{
	public OpdfLogNormal factor()
	{
		return new OpdfLogNormal();
	}
}
