package be.ac.ulg.montefiore.run.jahmm.io;

import java.io.IOException;
import java.io.StreamTokenizer;

import be.ac.ulg.montefiore.run.jahmm.OpdfLogNormal;


/**
 * This class implements a {@link OpdfLogNormal} reader.  The syntax of the
 * distribution description is the following.
 * <p>
 * The description always begins with the keyword <tt>LogNormalOPDF</tt>.
 * The next (resp. last) symbol is an opening (resp. closing) bracket.
 * Between the brackets are two numbers separated by a space.  The
 * first is the distribution's mean, the second the variance.
 * <p>
 * For example, reading <tt>LogNormalOPDF [ .2 .3 ]</tt> returns a distribution
 * equivalent to <code>new OpdfLogNormal(.2, .3)</code>.
 */
public class OpdfLogNormalReader
extends OpdfReader<OpdfLogNormal>
{
	String keyword()
	{
		return "LogNormalOPDF";
	}

	public OpdfLogNormal read(StreamTokenizer st)
	throws IOException,	FileFormatException {
		HmmReader.readWords(st, keyword());

		double[] meanVariance = OpdfReader.read(st, 2);

		return new OpdfLogNormal(meanVariance[0], meanVariance[1]);
	}
}
