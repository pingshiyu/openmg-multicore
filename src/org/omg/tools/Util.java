package org.omg.tools;

import org.omg.MolProcessor;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.NoSuchElementException;
import java.util.Scanner;

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
					weightArray[i][j] = MolProcessor.INFINITY;
					weightArray[j][i] = MolProcessor.INFINITY;
				}
			}
		}
		return weightArray;
	}
}
