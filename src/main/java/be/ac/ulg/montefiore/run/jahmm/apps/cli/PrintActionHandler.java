/*
 * Copyright (c) 2004-2009, Jean-Marc François. All Rights Reserved.
 * Licensed under the New BSD license.  See the LICENSE file.
 */

package be.ac.ulg.montefiore.run.jahmm.apps.cli;

import java.io.*;
import java.util.EnumSet;

import be.ac.ulg.montefiore.run.jahmm.*;
import be.ac.ulg.montefiore.run.jahmm.apps.cli.CommandLineArguments.Arguments;
import be.ac.ulg.montefiore.run.jahmm.io.*;


/**
 * Reads a HMM from a file and prints it in a human readable way.
 */
class PrintActionHandler extends ActionHandler
{
	@SuppressWarnings({"unchecked", "rawtypes"}) // We use a generic reader
	public void act()
	throws FileFormatException, IOException, FileNotFoundException,
	AbnormalTerminationException
	{
		EnumSet<Arguments> args = EnumSet.of(Arguments.IN_HMM);
		CommandLineArguments.checkArgs(args);
		
		InputStream in = Arguments.IN_HMM.getAsInputStream();
		OpdfReader opdfReader = new OpdfGenericReader();
		Hmm<?> hmm = HmmReader.read(new InputStreamReader(in), opdfReader);
		
		System.out.println(hmm);
	}
}
