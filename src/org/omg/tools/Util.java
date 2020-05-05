package org.omg.tools;

import org.omg.MolProcessor;
import org.omg.PMG;
import org.openscience.cdk.tools.SystemOutLoggingTool;

import java.util.*;

public class Util {
	public static ArrayList<String> parseFormula (String formula){
		boolean success = false;
		ArrayList<String> atoms = new ArrayList<String>();
		String sym = "";
		while (formula.length() > 0) {
			success = false;
			sym += formula.charAt(0);
			formula = formula.substring(1); 
			if (!Atom.valenceTable.containsKey(sym)) continue;
			success = true;
			int count=0;
			char ch;
			while (formula.length() > 0) {
				ch = formula.charAt(0);
				if ('0' <= ch && ch <= '9') {
					count = count*10 + (ch-'0');
					formula = formula.substring(1); 
				} else {
					break;
				}
			}
			if (count == 0) count = 1;
			for (int n=0; n<count; n++) {
				atoms.add(sym);
			}
			sym = "";
		}
		if (!success) {
			System.err.println("Could not parse the sub-expression: "+sym);
			return null;
		}
		return atoms;
	}
	
	public static int[] readFragment(final Scanner inFile, int[][] matrix, Atom[] atoms){
		int map[], reverseMapping[]=new int[atoms.length];
		Arrays.fill(reverseMapping, -1);
		if (!inFile.hasNext()) return null;
		try{
			inFile.nextLine();
			inFile.nextLine();
			inFile.nextLine();
			int atomCount = inFile.nextInt();
			int bondCount = inFile.nextInt();
			map = new int[atomCount];
			inFile.nextLine();
			for (int i=0; i<atomCount; i++){
				inFile.next(); inFile.next(); inFile.next();
				String symbol = inFile.next();
				inFile.nextLine();
				int pos;
				for (pos=0;!atoms[pos].symbol.equals(symbol); pos++);
				if (pos==atoms.length) {
					System.err.println("The atom "+symbol+" in the fragment does not exist in the molecule.");
					return null;
				}
				while (pos < atoms.length && reverseMapping[pos] != -1) pos++;
				if (pos == atoms.length || !atoms[pos].symbol.equals(symbol)) {
					System.err.println("The number of "+symbol+" atoms in the fragment is bigger than the original molecule.");
					return null;
				}
				map[i] = pos;
				reverseMapping[pos] = i;
			}
			for (int i=0; i<bondCount; i++){
				int left = inFile.nextInt()-1;
				int right = inFile.nextInt()-1;
				int order = inFile.nextInt();
				matrix[map[left]][map[right]] = order;
				matrix[map[right]][map[left]] = order;
				inFile.nextLine();
			}
			do {
				if (inFile.nextLine().equals("$$$$")) break;				
			} while (true);
		} catch (NoSuchElementException nse){
			System.err.println("Failed to read the SDF file.");
			return null;
		}
		return map;
	}

	public static void print2DIntMatrix(int[][] matrix) {
		System.out.println(Arrays.deepToString(matrix).replace("], ", "] \n"));
	}

	public static double[][] intToDouble2DArray(int[][] arr) {
		int n = arr.length;
		double[][] doubleArr = new double[n][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				doubleArr[i][j] = arr[i][j];
			}
		}
		return doubleArr;
	}

	public static int[][] doubleToInt2DArray(double[][] arr) {
		int n = arr.length;
		int[][] intArr = new int[n][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				intArr[i][j] = (int) Math.round(arr[i][j]);
			}
		}
		return intArr;
	}

	/**
	 * Given a square array A_ij, convert to indicator array for if the entry A_ij is nonzero.
	 * @param arr input array
	 * @return 1-0 indicator array of `arr`
	 */
	public static int[][] toOneZero2DArray(int[][] arr) {
		int n = arr.length;
		int[][] oneZeroArray = new int[n][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (arr[i][j] != 0) {
					oneZeroArray[i][j] = 1;
				}
			}
		}
		return oneZeroArray;
	}

	/**
	 * Convert adjacency matrix to weight array:
	 * in an adjacency matrix, if an edge is present then the weight is 1
	 * if an edge is not present then the weight will be infinity.
	 * @param arr
	 * @return
	 */
	public static int[][] adjacencyToWeight(int[][] arr) {
		int n = arr.length;
		int[][] weightArray = new int[n][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < i; j++) {
				if (arr[i][j] != 0) {
					weightArray[i][j] = 1;
					weightArray[j][i] = 1;
				} else {
					weightArray[i][j] = PMG.INFINITY;
					weightArray[j][i] = PMG.INFINITY;
				}
			}
		}
		return weightArray;
	}

	/**
	 * Split a string by individual capital letters (representing atoms):
	 * 	for a formula this would break the formula into components;
	 * @param atomsString a string with capital letters
	 * @return `atomsString` split according to capital letters
	 */
	public static String[] splitIntoComponents(String atomsString) {
		return atomsString.split("(?=\\p{Upper})");
	}

	/**
	 * Generating all valid molecular formulae under a certain `mass`, and using
	 * atoms from `atoms`.
 	 * @param mass maximum mass limit
	 * @param atoms set of atoms to comprise molecular with
	 * @return Set of formulae, under `mass` using `atoms`.
	 */
	public static Set<String> formulaeUnder(int mass, List<Character> atoms) {
		atoms.sort((a1, a2) -> {
			if (a1 == a2) return 0;
			else return Atom.atomWeights.get(a1) < Atom.atomWeights.get(a2) ? 1 : -1;
		});
		Set<String> formulae = new HashSet<>();
		traverseFormulaeTree(formulae, mass, atoms, "");
		return formulae;
	}

	private static void traverseFormulaeTree(Set<String> formulae, int mass, List<Character> atoms, String formula) {
		// Save current formula if legal.
		if (legal(formula)) {
			formulae.add(formula);
		}

		// No atoms available to use: tree ends
		if (atoms.isEmpty()) {
			return;
		}

		Character atom = atoms.get(0);
		// Pruning for hydrogen: when we are at Hydrogens, and we do not have free valencies.
		// We can terminate the tree here, since Hydrogens are the last atoms being added.
		if (atom == 'H' && !freeValencyCheck(formula))
			return;

		int atomWeight = Atom.atomWeights.get(atom);

		// Two choices: do we want `atom` in our formula, or not?
		// If we don't: (left subtree) then traverse with `atom` removed.
		List<Character> newAtoms = new ArrayList<>(atoms);
		newAtoms.remove(0);
		traverseFormulaeTree(formulae, mass, newAtoms, formula);

		// If we do (and can): then continue traversing the tree, updating the mass used
		if (atomWeight > mass) {
			return; // cannot include new atom: weight too low - tree ends here.
		}
		String childFormula = addAtomToFormula(formula, atom);
		traverseFormulaeTree(formulae, mass - atomWeight, atoms, childFormula);
	}

	/**
	 * Add an atom to a candidate chemical formulae (So adding 'H' to 'CH4' gives 'CH5')
	 * @return `formulae` with `atom` added into it
	 */
	private static String addAtomToFormula(String formula, char atom) {
		if (formula.isEmpty()) return formula + atom;

		// Separate into components by regex
		String[] components = splitIntoComponents(formula);

		// If the atom already exists in the formula: this would find it
		for (int i = 0; i < components.length; i++) {
			if (components[i].charAt(0) == atom) {
				if (components[i].length() == 1) {
					// Single atom: we add 1 to it
					components[i] = atom + "2";
				} else {
					// Not single atom: we increment the counter
					int numberOfAtoms = Integer.parseInt(components[i].substring(1)) + 1;
					components[i] = atom + "" + numberOfAtoms;
				}
				return String.join("", components);
			}
		}
		// If the atom does not exist in the formula:
		return formula + atom;
	}

	private static boolean freeValencyCheck(String formula) {
		if (formula.length() == 0)
			return false;

		int nonHDegreeSum = 0;
		int nH = 0;
		int n = 0;
		String[] components = splitIntoComponents(formula);
		// Read the components of the formula, and compute the degrees
		for (String component : components) {
			char atom = component.charAt(0);
			int nA;
			if (component.length() == 1) {
				nA = 1;
			} else {
				nA = Integer.parseInt(component.substring(1));
			}
			n += nA;

			int valency = Atom.valenceTable.get("" + atom).get(0);
			if (atom == 'H') {
				nH = nA;
			} else {
				nonHDegreeSum += valency * nA;
			}
		}

		// We cannot have more hydrogens than free valencies.
		// Most free valencies we can have is Sum_{A =/= H} n_A - 2*(n - n_H)
		// So need: n_H <= Sum_{A =/= H} val(A)*n_A - 2*(n - n_H)
		int nonH = n - nH;
		int maxFreeValencies = nonHDegreeSum - 2*(nonH - 1);
		if (nH > maxFreeValencies)
			return false;

		// all tests passed
		return true;
	}

	/**
	 * Determine if the formula corresponds a valid molecular graph.
	 * The condition checked is: whether the sum of all degrees (valencies) are even.
	 * @param formula, candidate formula (expanded unary version)
	 * @return boolean, whether the formula corresponds to a valid degree sequence
	 */
	private static boolean legal(String formula) {
		if (formula.length() <= 1)
			return false;

		int nonHDegreeSum = 0;
		int nH = 0;
		int n = 0;
		String[] components = formula.split("(?=\\p{Upper})");
		// Read the components of the formula, and compute the degrees
		for (String component : components) {
			char atom = component.charAt(0);
			int nA;
			if (component.length() == 1) {
				nA = 1;
			} else {
				nA = Integer.parseInt(component.substring(1));
			}
			n += nA;

			int valency = Atom.valenceTable.get("" + atom).get(0);
			if (atom == 'H') {
				nH = nA;
			} else {
				nonHDegreeSum += valency * nA;
			}
		}

		// Check for evenness. Graph cannot be legal if sum of degrees is odd.
		if ((nonHDegreeSum + nH) % 2 != 0)
			return false;

		// We cannot have more hydrogens than free valencies.
		// Most free valencies we can have is Sum_{A =/= H} n_A - 2*(n - n_H)
		// So need: n_H <= Sum_{A =/= H} val(A)*n_A - 2*(n - n_H)
		int nonH = n - nH;
		int maxFreeValencies = nonHDegreeSum - 2*(nonH - 1);
		if (nH > maxFreeValencies)
			return false;

		// all tests passed
		return true;
	}

	/**
	 * Convert a long formula (unary) into short form formula.
	 * e.g. CHHHH -> CH4
	 * @return short form formula String
	 */
	private static String toShortFormula(String formula) {
		// Store in a multiset (as a map)
		Map<Character, Integer> atomsCount = new HashMap<>();
		int n = formula.length();
		for (int i = 0; i < n; i++) {
			char atom = formula.charAt(i);
			int count = atomsCount.containsKey(atom) ? atomsCount.get(atom) : 0;
			atomsCount.put(atom, ++count);
		}

		// Convert the map into a string.
		String shortFormula = "";
		for (Character atom : atomsCount.keySet()) {
			shortFormula += (atomsCount.get(atom) == 1) ? "" + atom : "" + atom + atomsCount.get(atom);
		}

		return shortFormula;
	}
}
