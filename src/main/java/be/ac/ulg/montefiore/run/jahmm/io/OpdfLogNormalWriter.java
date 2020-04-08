package be.ac.ulg.montefiore.run.jahmm.io;

import java.io.IOException;
import java.io.Writer;

import be.ac.ulg.montefiore.run.jahmm.OpdfLogNormal;


/**
 * This class implements a {@link OpdfLogNormal} writer.  It is compatible
 * with the {@link OpdfLogNormalReader} class.
 */
public class OpdfLogNormalWriter
extends OpdfWriter<OpdfLogNormal>
{
	public void write(Writer writer, OpdfLogNormal opdf)
	throws IOException
	{
		String s = "LogNormalOPDF [";

		s += opdf.mean() + " " + opdf.variance();

		writer.write(s + "]\n");
	}
}
