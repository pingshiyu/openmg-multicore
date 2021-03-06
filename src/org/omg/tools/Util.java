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
	 * Check, using the Chungphaisan-Erdos-Gallai algorithm for multigraphical degree sequence
	 * Let `d_i` be the degree of the ith vertex, in decreasing order. (d_1 >= d_2 ... >= d_n)
	 * b be the maximal edge multiplicity;
	 * H_i be the sum of d_1 + ... + d_i
	 * w_i be the Erdos-Gallai-like weights for multigraphs: for i \in 1, ..., n-1
	 * 	w_i = 0 if d_1 < b*i
	 * 	w_i = k such that k is the largest with d_k >= b*i
	 *
	 * The degree sequence is multi-graphical with maximum degree b if and only if:
	 * 	1. sum of degrees is even
	 * 	2. H_i <= b * i * (y_i - 1) + H[n] - H[y_i], where y_i = max(i, w_i)
	 *
	 * @param d degree sequence in decreasing order
	 * @param b the maximal edge multiplicity
	 * @return boolean, whether the degree sequence is multi-graphical
	 */
	private static boolean legalMultigraph(int[] d, int b) {
		int n = d.length;
		int[] H = new int[n+1];
		H[0] = 0;
		for (int i = 1; i <= n; i++) {
			H[i] = H[i-1] + d[i-1];
		}

		// check sum of degrees is even
		if (H[n] % 2 == 1) {
			return false;
		}

		int w = n;
		for (int i = 1; i <= n-1; i++) {
			// Compute w_i (taking advantage of w_i being strictly decreasing
			while(w > 1 && d[w-1] < b*i)
				w--;

			int y = Math.max(i, w);

			if (H[i] > (b*i*(y - 1) + H[n] - H[y]))
				return false;
		}
		return true;
	}

	private static boolean legal(String formula) {
		if (formula.length() == 0)
			return false;
		int[] degreeSequence = toDescendingDegreeSequence(formula);
		return legalMultigraph(degreeSequence, degreeSequence[0]);
	}

	/**
	 * Convert the formula into a degree sequence, in decreasing order
	 * @param formula chemical formula
	 * @return degree sequence, in descending order
	 */
	private static int[] toDescendingDegreeSequence(String formula) {
		String[] components = splitIntoComponents(formula);
		List<Integer> degrees = new ArrayList<>();
		for (String component : components) {
			String atom = "" + component.charAt(0);
			int nA = 1;
			if (component.length() > 1) {
				nA = Integer.parseInt(component.substring(1));
			}

			List<Integer> valencies = Atom.valenceTable.get(atom);
			Integer valency = valencies.get(valencies.size() - 1);
			List<Integer> atomDegrees = Collections.nCopies(nA, valency);
			degrees.addAll(atomDegrees);
		}

		Collections.sort(degrees, Collections.reverseOrder());
		return degrees.stream().mapToInt(i->i).toArray();
	}

	/**
	 * Determine if the formula corresponds a valid molecular graph.
	 * The condition checked is: whether the sum of all degrees (valencies) are even.
	 * @param formula, candidate formula (expanded unary version)
	 * @return boolean, whether the formula corresponds to a valid degree sequence
	 */
	/**
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
	 **/

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

	/**
	 * Computes the crowding bounds automatically based on the H-ratio of the formula
	 * @param formula molecular formula
	 * @param d distance to compute bounds up to
	 * @param aggressive use of aggressive bounds
	 * @return neighbourhoods restriction based on h-ratios
	 */
	public static int[] formulaHRatioCrowdingBounds(String formula, int d, boolean aggressive) {
		double hratio = formulaHRatio(formula);
		return hRatioBounds(hratio, d, aggressive);
	}

	/**
	 * Generates the crowding bounds based on the hydrogen ratio. Bounds goes up to
	 * `d`=3 and has the choice of aggression / no aggression
	 * @param hratio hydrogen ratio of the formula
	 * @param aggressive aggression, leading to different set of bounds
	 * @return
	 */
	private static int[] hRatioBounds(double hratio, int d, boolean aggressive) {
		double[] ratios = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0};

		// define the boundaries based on aggression, for d=2 and d=3
		int[][] points;
		if (aggressive) {
			points = new int[][] {{7, 8},  {9, 13},  {12, 17}, {14, 21}, {15, 24}, {15, 24}, {14, 20}, {9, 11},  {9, 11}};
		} else {
			points = new int[][] {{9, 13}, {11, 17}, {13, 20}, {15, 24}, {15, 26}, {15, 26}, {14, 22}, {12, 12}, {12, 12}};
		}

		// crowding bounds
		int[] bounds = new int[d+1];
		if (d < 2) return new int[] {1, 6};
		bounds[0] = 1; bounds[1] = 6; // max. valency 6

		for (int i = 0; i < ratios.length-1; i++) {
			double ratio_lb = ratios[i]; double ratio_ub = ratios[i+1];

			if ((ratio_lb <= hratio) & (hratio <= ratio_ub)) {
				// assigns values according by interpolation
				for (int k = 2; k <= d; k++) {
					int crowding_lb = points[i][k-2];
					int crowding_ub = points[i+1][k-2];
					bounds[k] = (int) Math.floor(crowding_lb + (hratio-ratio_lb)/(ratio_ub-ratio_lb) * (crowding_ub-crowding_lb));
				}
				return bounds;
			}
		}
		return null;
	}

	/**
	 * Compute the H-ratio in a formula
	 * @param formula a string representing a molecular formula
	 * @return double, the h-ratio
	 */
	private static double formulaHRatio(String formula) {
		String[] components = splitIntoComponents(formula);
		int nh = 0;
		int n = 0;
		for (String component : components) {
			char atom = component.charAt(0);

			int na;
			if (component.length() == 1) na = 1;
			else na = Integer.parseInt(component.substring(1));

			if (atom == 'H') nh += na;
			n += na;
		}

		return ((double) nh) / n;
	}
}
