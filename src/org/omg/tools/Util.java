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
	 * Generate all valid molecular formulae under a certain molecular `mass`.
	 * Taking atoms from `atoms`.
	 */
	public static Set<String> formulaeUnder(int mass){
		// Dynamic programming: compute all solutions under mass `mass`:
		HashMap<Integer, Set<String>> allSolutions = new HashMap<>();

		// Base case: there are no formulae of mass <= min_weight => emptySet.
		int minWeight = PMG.INFINITY;
		for (char atom : Atom.atoms) {
			if (Atom.atomWeights.get(atom) < minWeight)
				minWeight = Atom.atomWeights.get(atom);
		}
		for (int m = 0; m < minWeight; m++) {
			HashSet<String> noFormulaSet = new HashSet<String>();
			noFormulaSet.add("");
			allSolutions.put(m, noFormulaSet);
		}

		/*
		Form the larger set by equation:
		L(w) = U_{w_i <= w} ( {S u {a_i} : S \in L(w - w_i)} u L(w - w_i) )
		 */
		for (int x = minWeight; x <= mass; x++) {
			HashSet<String> solutionsOfX = new HashSet<String>();

			// Higher weight solution also includes all lower weight solutions
			// We only need to add in the solutions of the highest weight subproblem.
			solutionsOfX.addAll(allSolutions.get(x - minWeight));

			// Also includes all possibilities of adding new atoms, as long as they are of weight <= x
			for (char atom : Atom.atoms) {
				int atomWeight = Atom.atomWeights.get(atom);

				// add all suitable atoms to the formulae
				if (atomWeight <= x) {
					// And contains all possibilities of adding the extra atom.
					for(String formula : allSolutions.get(x - atomWeight)) {
						solutionsOfX.add(addAtomToFormula(formula, atom));
					}
				}
			}
			allSolutions.put(x, solutionsOfX);
		}

		// The candidate solutions may not be valid graphs: we need to make sure the degree sequence sum is even,
		// hence legal.
		Set<String> candidateSolutions = allSolutions.get(mass);
		Set<String> finalSolutions = new HashSet<String>();
		for (String formula : candidateSolutions) {
			if (legal(formula))
				finalSolutions.add(toShortFormula(formula));
		}

		return finalSolutions;
	}

	/**
	 * Add an atom to a candidate chemical formulae (list out all of the atoms one by one, sorted by
	 * alphabetical order: e.g. CH4 is written as CHHHH.)
	 * @return `formulae` with `atom` added into it
	 */
	private static String addAtomToFormula(String formula, char atom) {
		for (int i = 0; i < formula.length(); i++) {
			if (formula.charAt(i) >= atom) {
				return formula.substring(0, i) + atom + formula.substring(i);
			}
		}
		return formula + atom;
	}

	/**
	 * Determine if the formula corresponds a valid molecular graph.
	 * The condition checked is: whether the sum of all degrees (valencies) are even.
	 * @param formula, candidate formula (expanded unary version)
	 * @return boolean, whether the formula corresponds to a valid degree sequence
	 */
	private static boolean legal(String formula) {
		int n = formula.length();
		// No single atom molecule.
		if (n <= 1)
			return false;

		boolean even = true;
		int nonHDegreeSum = 0;
		int nH = 0;
		for (int i = 0; i < n; i ++) {
			char atom = formula.charAt(i);
			// odd valencies => flip evenness
			if (atom == 'H' || atom == 'N' || atom == 'P')
				even = !even;

			// Count free valencies & leaves (H) counts.
			if (atom == 'H')
				nH++;
			else
				nonHDegreeSum += Atom.valenceTable.get("" + atom).get(0);
		}

		// Check for evenness. Graph cannot be legal if sum of degrees is odd.
		if (!even)
			return false;

		// Hn for n > 2 is illegal.
		if (nH > 2 && nH == n)
			return false;

		// Check for connectedness: in a graph of not just hydrogens, we cannot have more hydrogens than free valencies.
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
