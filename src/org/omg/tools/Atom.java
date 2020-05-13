package org.omg.tools;

import java.util.*;
import java.util.prefs.AbstractPreferences;

public class Atom {
	public static final HashMap<String,  List<Integer>> valenceTable;
	public static final HashMap<Character, Integer> atomWeights;
	public static final Set<Character> atoms;
	public final String symbol;
	public final int maxValence;
	public final List<Integer> valenceList;
	public boolean flag = false;	// can be used for different purposes, e.g., to mark the atom being in a fragment
	
	public Atom (String s, int id) {
		symbol = s;
		valenceList = valenceTable.get(symbol);
		if (valenceList == null) {
			System.err.println("The atom type "+symbol+" is not supported.");
			System.exit(200);
		}
		maxValence = valenceList.get(0);
	}

	static {
		// initialize the table
		valenceTable = new HashMap<String,  List<Integer>>();
		// TODO: read atom symbols from CDK?
		valenceTable.put("H", Arrays.asList(1));
		valenceTable.put("C", Arrays.asList(4));
		valenceTable.put("N", Arrays.asList(3));	// Remember to put the biggest number first!?: Not including 5 for now.
		valenceTable.put("O", Arrays.asList(2));
		valenceTable.put("S", Arrays.asList(6,4,2));
		valenceTable.put("P", Arrays.asList(5,3));
	}

	static {
		atomWeights = new HashMap<Character, Integer>();
		atomWeights.put('C', 12);
		atomWeights.put('H', 1);
		atomWeights.put('N', 14);
		atomWeights.put('O', 16);
		atomWeights.put('P', 31);
		atomWeights.put('S', 32);
	}

	static {
		atoms = new HashSet<Character>();
		atoms.add('C');
		atoms.add('H');
		atoms.add('N');
		atoms.add('O');
		atoms.add('P');
		atoms.add('S');
	}

}
