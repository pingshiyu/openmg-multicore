package org.omg;

import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Scanner;
import java.util.Set;
import java.util.concurrent.atomic.AtomicLong;

import org.omg.tools.Atom;
import org.omg.tools.CDKUtil;
import org.omg.tools.Util;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IChemSequence;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import fi.tkk.ics.jbliss.Graph;

// matrix multiplication / linear algebra library
import org.ejml.simple.SimpleMatrix;

/**
 * Meat of the program.
 * Main class for generating all molecules.
 * Given molecular formula passed in from main function.
 */

public class MolProcessor implements Runnable{
	static final int SEM_CAN = 1;
	static final int MIN_CAN = 2;
	static final int CAN_AUG = 4;
	static final int OPTIMAL = SEM_CAN | MIN_CAN;	// currently mix of sem_can + min_can
	public static final int INFINITY = 99999;
	final int method;		// can be one of the above choices
	final boolean checkBad;
	final boolean hashMap;	// used with semiCan (otherwise, minimality check) for finished molecules
	final boolean cdkCheck;
    final private IAtomContainerSet goodlistquery;
    final private IAtomContainerSet badlistquery;

	// all non-hydrogen atoms are stored.
	public final Atom[] atoms;
	final int nH;
	final Graph graph;

	// neighbourhood constraint function: neighbourhood[i] gives the max. number of neighbours within distance i.
	final int[] neighbourhood;
	final boolean restrictNeighbourhoods;

	// These are all n*n matrices.
    // adjacency_{ij}: indicates the multiplicity of the bond connecting atom i to atom j. 0 or more.
    // fragment: the molecule used as the base for generation
    // connectivity_{ij} indicates how many ways atom i is connected to atom j
    // loopPart_{ij}: in how many loops edge ij participates
    // distance_{ij}: the distance from atom i to atom j
    // final: will always refer to the same array and will stay consistent throughout threads.
	final int[][] adjacency, fragment, connectivity, loopPart, distance;
	final static AtomicLong duplicate = new AtomicLong(0);
	private final IAtomContainer acontainer;
//	final boolean compatible[][][];

	static boolean frag = false;
	int []blocks;
	int maxOpenings;
    // startLeft / startRight indicates the indices of an atom
	int startLeft;
	int startRight;
	String canString="";
	
	private static final Set<String> molSet = Collections.synchronizedSet(new HashSet<String>());
	private static final int MAX_LOOP_SIZE = 8;
	

//
//	void initCompatible() {
//		for (int i=0; i<atoms.length; i++){
//			compatible[i][i][0] = true;
//			for (int j=i+1; j<atoms.length; j++){
//				compatible[i][j][0] = compatible[j][i][0] = (atoms[i].symbol.equals(atoms[j].symbol) && atoms[i].flag == atoms[j].flag); // && adjacency[i][j] == adjacency[j][i];
//			}
//		}
//	}

	/**
	 * Used with the semi-canonization method, which does not keep 
	 * a canonical form, but a semi-canonical form of the molecule.
	 * 
	 * @param atoms
	 * @param nH
	 * @param maxOpenings
	 * @param adjacency
	 * @param gr
	 */
	public MolProcessor(final Atom[] atoms, final int nH, final int maxOpenings, final int[] neighbourhood, final boolean restrictNeighbourhoods,
						final int[][] adjacency, final int[][]fragment, final int[][] connectivity, final int[][] loopPart, final int[][] distance,
			            final Graph gr, String canStr, final int[] blocks,
			            final int stL, final int stR, final IAtomContainer acontainer,
			            final int method, final boolean hm, final boolean cdk, final boolean checkBad, 
			            final IAtomContainerSet good, final IAtomContainerSet bad) {
		this.method = method;
		this.hashMap = hm;
		this.cdkCheck = cdk;
		this.checkBad = checkBad;
		this.canString = canStr;
		if (blocks == null) this.blocks = null;
		else this.blocks = blocks.clone();
		this.atoms = atoms;
		this.goodlistquery = good;
		this.badlistquery = bad;
//		this.compatible = new boolean[atoms.length][atoms.length][atoms.length];
//		initCompatible();
		this.nH = nH;
		this.maxOpenings = maxOpenings;
		graph = gr;
		this.fragment = fragment;
		this.restrictNeighbourhoods = restrictNeighbourhoods;
		if (restrictNeighbourhoods) {
			this.neighbourhood = new int[neighbourhood.length];
			System.arraycopy(neighbourhood, 0, this.neighbourhood, 0, neighbourhood.length);
		} else {
			this.neighbourhood = null;
		}
		this.adjacency    = new int [atoms.length][atoms.length];
//		if (fragment != null) 
//			this.fragment = new int [atoms.length][atoms.length];
//		else 
//			this.fragment = null;
		this.connectivity = new int [atoms.length][atoms.length];
		this.loopPart     = new int [atoms.length][atoms.length];
		this.distance     = new int [atoms.length][atoms.length];
		for (int i=0; i<atoms.length; i++)
			for (int j=0; j<atoms.length; j++){
//				if (fragment != null) this.fragment[i][j] = fragment[i][j];
				this.adjacency[i][j]    = adjacency[i][j];
				this.connectivity[i][j] = connectivity[i][j];
				this.loopPart[i][j]     = loopPart[i][j];
				this.distance[i][j]     = distance[i][j];
			}
		this.startLeft = stL;
		this.startRight = stR;
		this.acontainer = acontainer;
	}
	
	public MolProcessor(final ArrayList<String> atomSymbols, String formula, int[] neighbourhood, boolean restrictNeighbourhoods,
			final int method, final boolean hm, final boolean cdk, final boolean checkBad, final boolean frag, String goodlist, String badlist) throws CDKException, FileNotFoundException{
		this.method = method;
		this.hashMap = hm;
		this.cdkCheck = cdk;
		this.checkBad = checkBad;
		int nH=0; // number of hydrogens
		int maxOpenings=0;
		int atomCount=0;
		for (String symbol:atomSymbols){
			if (symbol.equals("H")) {
				nH++;
				maxOpenings--;
			}
		}
		atoms = new Atom[atomSymbols.size()-nH];
		this.nH = nH;
		// populates array with all atoms, except for hydrogens.
		for (String symbol:atomSymbols){
			if (symbol.equals("H")) continue;
			atoms[atomCount] = new Atom(symbol,atomCount); // note: second parameter is not used at all.
			maxOpenings += atoms[atomCount].maxValence;
			atomCount++;
		}

		// initialise neighbourhood constriant
		this.restrictNeighbourhoods = restrictNeighbourhoods;
		if (restrictNeighbourhoods) {
			this.neighbourhood = new int[neighbourhood.length];
			System.arraycopy(neighbourhood, 0, this.neighbourhood, 0, neighbourhood.length);
		} else {
			this.neighbourhood = null;
		}

		// initialises adjacency matrix, initially to all 0.
		adjacency = new int [atomCount][atomCount];
		if (frag) 
			this.fragment = new int [atoms.length][atoms.length];
		else 
			this.fragment = null;
		loopPart = new int [atomCount][atomCount];
		// Connectivity matrix, initialised to all -1
		connectivity = new int [atomCount][atomCount];
		for (int i=0; i<atoms.length; i++) Arrays.fill(connectivity[i], -1);
		// Distance matrix, all off diagonal entries and infty, and all diagonal entries are 0 initially.
		distance = new int [atomCount][atomCount];
		for (int i=0; i<atoms.length; i++) Arrays.fill(distance[i], INFINITY);
        for (int i=0; i<atoms.length; i++) distance[i][i] = 0;
		
		graph = new Graph(atomCount);
		this.maxOpenings = maxOpenings;
		startLeft=0;
		startRight=1;
		
//		this.compatible = new boolean[atoms.length][atoms.length][atoms.length];
//		initCompatible();
		blocks = new int [atoms.length];
		for (int i=2; i<atoms.length; i++) {
			if (atoms[i].symbol.equals(atoms[i-1].symbol))
				blocks[i] = 1;	// continuing the block
		}
		
		if (cdkCheck) {
			IAtomContainer lcontainer;
			do {	// this stupid loop is here only because CDK is not stable with respect to the order of the atoms it returns
				lcontainer = MolecularFormulaManipulator.getAtomContainer(
					MolecularFormulaManipulator.getMolecularFormula(formula, DefaultChemObjectBuilder.getInstance()));
			}while(!inRightOrder(atomSymbols, lcontainer));	// make sure the order is as we want. lcontainer is CDK class for atoms
			acontainer = lcontainer;
				
	        if(goodlist != null){
	            InputStream ins = new BufferedInputStream(new FileInputStream(goodlist));
	            MDLV2000Reader reader = new MDLV2000Reader(ins);
	            ChemFile fileContents = (ChemFile)reader.read(new ChemFile());               
	            IChemSequence sequence = fileContents.getChemSequence(0);

	            goodlistquery = sequence.getChemModel(0).getMoleculeSet();
	            for (int i=1; i<sequence.getChemModelCount(); i++) {
	                goodlistquery.add(sequence.getChemModel(i).getMoleculeSet());   
	            }
	        }  else 
	        	this.goodlistquery = null;
	        if(badlist != null){
	        	 InputStream ins = new BufferedInputStream(new FileInputStream(badlist));
	             MDLV2000Reader reader = new MDLV2000Reader(ins);
	             ChemFile fileContents = (ChemFile)reader.read(new ChemFile());   

	             IChemSequence sequence = fileContents.getChemSequence(0);
	             badlistquery = sequence.getChemModel(0).getMoleculeSet();
	             for (int i=1; i<sequence.getChemModelCount(); i++) {
	                 badlistquery.add(sequence.getChemModel(i).getMoleculeSet());   
	             }
	        } else 
	        	this.badlistquery = null;
		} else {
			acontainer = null;
			goodlistquery = null;
			badlistquery = null;
		}
	}

	private boolean inRightOrder(ArrayList<String> atomSymbols, IAtomContainer lcontainer) {
		List<IAtom> listcont = new ArrayList<IAtom>();
		for(IAtom atom: lcontainer.atoms()){
			if(atom.getSymbol().equals("H")){
				listcont.add(atom);
			}
		}
		for(IAtom atom: listcont){
			lcontainer.removeAtom(atom);
		}

		int atom = 0;
		for (String symbol:atomSymbols){
			if (symbol.equals("H")) continue;	// skip hydrogens
			String symbol2 = lcontainer.getAtom(atom++).getSymbol();
			if (!symbol2.equals(symbol)) return false;	// compare other atoms
		}
		return true;
	}

	private String molString(int[] canPerm){
        /**
         * Obtains the canonical representation of the molecule.
         * The canonical representation is the lower triangle of the matrix, read left-to-right
         * Which gives the molecule's canonical string.
         * Input:
         *  canPerm: canonical permutation (labelling) of the atoms / nodes
         */
		int[][] adjacency = new int [atoms.length][atoms.length];
		for (int i=0; i<atoms.length; i++)
			for (int j=0; j<atoms.length; j++){
				adjacency[canPerm[i]][canPerm[j]] = this.adjacency[i][j];
			}
		StringBuilder buf = new StringBuilder();
		for (int i=0; i<atoms.length; i++)
			for (int j=0; j<i; j++){
				buf.append(adjacency[i][j]);
			}
		return buf.toString();
	}

	public int getBondOrder (final int l, final int r){
		return adjacency[l][r];
	}

	/**
	 * This recursive function takes account of different valencies for the atoms
	 */
	private int openings;
	private boolean addUpOpenings(final int atomNum) {
		if (atoms.length == atomNum) 
			return openings == nH;
		else {
			int bondSum = 0;
			for (int i=0;i<atoms.length;i++)
				bondSum += adjacency[atomNum][i];
			for (Integer valence:atoms[atomNum].valenceList) {
				if (valence >= bondSum) {
					openings += (valence - bondSum);
					if (openings <= nH && addUpOpenings(atomNum+1)) {
						return true;
					}
					openings -= (valence - bondSum);
				}
			}
		}
		return false;
	}
	
	public boolean isComplete() {
		if (maxOpenings == 0) return true;
		openings = 0;
		return addUpOpenings (0);
	}

	private boolean isConnected() {
        // connectivity check done by checking whether connected to first atom or not
		for (int i=1; i<atoms.length; i++) if (connectivity[0][i]==-1) return false;
		return true;
	}
	
	private boolean[] seen;
	private boolean isConnectedDFS(){
		seen = new boolean[atoms.length];
		dfs(0);
		boolean connected = true;
		for (boolean b : seen) connected &= b;
		return connected;
	}
	
	private void dfs(int num) {
		seen[num] = true;
		for (int i=0; i<atoms.length; i++){
			if (adjacency[i][num]>0 && !seen[i]) dfs(i);
		}
	}
	
//
//	private boolean incBondSemiCanBig(final int left, final int right){
//		if (right<atoms.length-1 && atoms[left].flag == atoms[right].flag && 
//			atoms[right].symbol.equals(atoms[right+1].symbol) && 
//			adjacency[left][right]==adjacency[left][right+1]) {
//			int row;
//			for (row=left; row>0; row--)  {
//				if (adjacency[row-1][right]!=adjacency[row-1][right+1]) break;
//			}
//			if (row == 0)
//				return false;
//		}
//		return incBond(left, right);
//	}
//	
//	private boolean incBondSemiCan(final int left, final int right){
//		if (right>left+1 && !atoms[right].flag && 
//			atoms[right].symbol.equals(atoms[right-1].symbol) && 
//			adjacency[left][right]==adjacency[left][right-1]) {
//			int row;
//			for (row=left; row>0; row--)  {
//				if (adjacency[row-1][right]!=adjacency[row-1][right-1]) break;
//			}
//			if (row == 0)
//				return false;
//		}
//		return incBond(left, right);
//	}
//	

	/**
	 * Add / remove bond between `left` and `right`, doing the necessary checks for bad structures as well as
	 * maintaining the matrics along the way (matrices representing current molecule)
	 * @param left
	 * @param right
	 * @return true / false for whether the bond was incremented / deleted
	 */
	private boolean incBond(final int left, final int right) {
		if ((method & SEM_CAN) != 0 && !(blocks[right]==0 || adjacency[left][right]<adjacency[left][right-1])) return false;
		if (adjacency[left][right] > 2) return false;	// no more than Triple bonds
		if (isFull(left))  return false;
		if (isFull(right)) return false;
		if (checkBad && beforeLoop(left, right)) return false;
		if (adjacency[left][right] == 0 && connectIfNotBad(left, right)) return false;
		
		adjacency[left][right]++;
		adjacency[right][left]++;
		return true;
	}

	private boolean decBond(final int left, final int right) {
		if (adjacency[left][right] <= 0) return false;	
		adjacency[left][right]--;
		adjacency[right][left]--;
		if (adjacency[left][right]==0) disconnect(left, right);
		return true;
	}

	/**
	 * connects `left` to `right` atom, and returns whether this new bond closes a loop.
     * Whilst maintaining the global constants, e.g. matrices for the current molecule
	 * @param left
	 * @param right
	 * @return true -> bad connection, not connected. false -> was ble to connect
	 */
	private boolean connectIfNotBad(int left, int right) {
		boolean loop = connectivity[left][right] != -1;
		adjacency[left][right]++;
		adjacency[right][left]++;

		if (loop && checkBad && afterLoop(left, right)) {
			adjacency[left][right]--;
			adjacency[right][left]--;
			return true;
		}
		adjacency[left][right]--;
		adjacency[right][left]--;

		if (restrictNeighbourhoods) {
			// Check for distance matrix adherence & update if satisfies constraints
			int[][] distanceMatrixUpdated = distanceMatrixAddBond(left, right);
			if (distanceMatrixUpdated != null) {
				// Distance matrix has been changed with addition of the new bond. We check if it is still compliant with
				// distance neighbourhood constraints.
				if (isNeighbourhoodCompliant(distanceMatrixUpdated)) {
					// compliant: we go ahead with the updated matrix
					System.arraycopy(distanceMatrixUpdated, 0, this.distance, 0, this.atoms.length);
				} else {
					// not compliant: we do not go ahead and indicate this bond should not be added
					return true;
				}
			}
		}

		if (loop) {
			// if the nodes are already connected, which means we will now make a loop, we don't mark it here, because 
			// it makes it easier to open the loop: then we don't need to do anything! We can however detect that here
			// a loop was made because the two nodes are not directly connected! Note that we remove edges as a result 
			// of backtracking, and therefore the edges are removed in the reverse order of their addition!
			// mark the loop, so we know also later
			int prev_i = right;
			for (int i=left; i!=right; prev_i = i, i = connectivity[i][right]){
				loopPart[prev_i][i]++;
				loopPart[i][prev_i]++;
			}
			loopPart[prev_i][right]++;	// also count the last edge to right
			loopPart[right][prev_i]++;	// also count the last edge to right
		} else {
			// add a direct connection between the two nodes
			connectivity[right][left] = left;	
			connectivity[left][right] = right;	
			// discover other connectivities (there would be none in case of closing a loop)
			for (int i=0; i<atoms.length; i++){
				if (i!=right && connectivity[left][i] != -1){
					connectivity[i][right] = connectivity[i][left]; 
					connectivity[right][i] = left;
				}
			}
			for (int i=0; i<atoms.length; i++){
				if (i!=left && connectivity[right][i] != -1 && connectivity[right][i] != left){
					connectivity[i][left] = connectivity[i][right]; 
					connectivity[left][i] = right;
				}
			}
			for (int i=0; i<atoms.length; i++){
				for (int j=0; j<atoms.length; j++){
					if (i!=right && connectivity[left][i] != -1 && connectivity[left][i] != right){
						if (j!=left && connectivity[right][j] != -1 && connectivity[right][j] != left){
							connectivity[i][j] = connectivity[i][left];
							connectivity[j][i] = connectivity[j][right];
						}
					}
				}
			}
		}
		return false;
	}

	/**
	 * Calculate the distance matrix after adding bond between `left` and `right`. Returning a copy of the matrix.
	 * If the bond will not change neighbourhood (i.e. only increases multiplicity of edge), return `null`.
	 * @param left
	 * @param right
	 * @return `null` if does not change distance matrix, or a copy of the updated distance matrix otherwise.
	 */
	private int[][] distanceMatrixAddBond(int left, int right) {
		// if we are only increasing bond multiplicity: distance matrix won't change.
		if (adjacency[left][right] > 0) {
			return null;
		}
		int[][] updatedMatrix = new int[atoms.length][atoms.length];
		for (int i = 0; i < updatedMatrix.length; i++) {
			for (int j = 0; j < i; j++) {
				// upon adding an edge, if the distance were to change, then it must be because of a shorter path
				// using the newly added edge. The newly added edge is (left, right). Hence the shortest path will be
				// minimum of all paths that either 1) uses (left, right), or 2) does not use (left, right).
				updatedMatrix[i][j] = Math.min(distance[i][right] + 1 + distance[left][j],
						Math.min(distance[i][j], distance[i][left] + 1 + distance[right][j]));
				updatedMatrix[j][i] = updatedMatrix[i][j];
			}
			// distance of i to oneself is 0.
			updatedMatrix[i][i] = 0;
		}
		return updatedMatrix;
	}

	/**
	 * Check if the distance matrix is compliant with the local neighbourhood constraint
	 * Constraint as defined in this.neighbourhood
	 * @param distanceMatrix
	 * @return true if is indeed compliant with neighbourhood constraints, false otherwise.
	 */
	private boolean isNeighbourhoodCompliant(int[][] distanceMatrix) {
		// calculate each of the `i`-neighbourhood sizes, and see if any exceeds the `neighbourhood` limit at any point.
		for (int i = 0; i < this.neighbourhood.length; i++) {
			int[] iNeighbourhood = kNeighbourhood(distanceMatrix, i);
			for (int j = 0; j < iNeighbourhood.length; j++) {
				// check if i-neighbourhood of entry (atom) j exceeds the limitation neighbourhood[i].
				if (iNeighbourhood[j] > neighbourhood[i]) return false;
			}
		}
		return true;
	}

	/**
	 * Calculates the neighbourhood sizes of each molecule within distance k. If k is greater than the number of atoms,
	 * (more strictly the diameter of the graph), then we will just return up to the atom count.
	 * @param distanceMatrix
	 * @param k
	 * @return int[] of the neighbourhood sizes of each molecule.
	 */
	private int[] kNeighbourhood(int[][] distanceMatrix, int k) {
		int[] neighbourhoodSizes = new int[distanceMatrix.length];
		for (int i = 0; i < distanceMatrix.length; i++) {
			// calculate the neighbourhood sizes for each atom.
			for (int j = 0; j < distanceMatrix.length; j++) {
				if (distanceMatrix[i][j] <= k) neighbourhoodSizes[i]++;
			}
		}
		return neighbourhoodSizes;
	}

	private void test() {
		for (int i=0; i<atoms.length; i++)
			for (int j=0; j<atoms.length; j++) {
				int j2 = connectivity[i][j];
				if (j2 != -1 && adjacency[i][j2] == -1)
					i = -2;
			}
	}

	private boolean afterLoop(int left, int right) {
		int prev_i=left;
		int i=right; 
		int singleInLoop = 0, doubleInLoop = 0;
		do {	
			int next_i = connectivity[i][left];
			if (next_i == -1) next_i = right; 	// close the loop (logically only)
			singleInLoop = doubleInLoop = 0;
			if (adjacency[prev_i][i] == 3) return true;		// triple in loop is always bad
			for (int x=0; x<atoms.length; x++) {
				if (x == i) continue;
				if (loopPart[x][i] > 0 || x == prev_i || x == next_i) {	// we haven't marked the current loop, so do consider them
					if (adjacency[x][i] == 1) singleInLoop++;
					else if (adjacency[x][i] == 2) doubleInLoop++;
				}
			}
			// detect: - two doubles, - one double and three singles, in any loop and connected to i
			if (doubleInLoop>1 || doubleInLoop>0 && singleInLoop>2) return true;
			prev_i=i; 
			i=next_i;
		} while (i != right); // "right"  is reached only by i=next_i above after loop is complete

		// Detect a square with two double bonds
		for (i=0; i<atoms.length; i++){
			if (i == right || i == left) continue;
			for (int j=0; j<atoms.length; j++){
				if (j==left || j==right) continue;
				if (adjacency[i][left] * adjacency[j][right] * adjacency[i][j] >= 4) 
					return true;
			}
		}
		return false;
	}

    /**
     * Analogous to connect. Disconnect `left` to `right` atom
     * Whilst maintaining the global constants, e.g. matrices for the current molecule
     * @param left
     * @param right
     */
	private void disconnect(int left, int right) {
		// Update the distance matrix with the `adjacency` matrix, which already the bond (left, right) decremented
		if (restrictNeighbourhoods)
			updateDistanceMatrix();

//		++c37;
		//		System.out.print(" -"+left+right);
		if (connectivity[left][right] != right) {	// at this point a loop was closed; see connect...() method
//			System.out.println(" - loop: "+(--loopCount));
//			if (loopCount<0)
//				loopCount--;
			int prev_i = right;
			for (int i=left; i!=right; prev_i = i, i = connectivity[i][right]){
				loopPart[prev_i][i]--;
				loopPart[i][prev_i]--;
			}
			loopPart[prev_i][right]--;	// also count the last edge to right
			loopPart[right][prev_i]--;	// also count the last edge to right
//			connectivity[right][left] = prev_i;
		} else {
			connectivity[left][right] = connectivity[right][left] = -1;
			for (int i=0; i<atoms.length; i++) {
				if (connectivity[right][i] == left) {
//					if (right==3 && i==7 || i==3 && right==7) 
//						System.out.println("removing 3-7 at "+left+" "+right+" --> "+c37++);
					connectivity[right][i] = connectivity[i][right] = -1;
				}
			}
			for (int i=0; i<atoms.length; i++) {
				if (connectivity[left][i] == right) {
//					if (left==3 && i==7 || i==3 && left==7) 
//						System.out.println("removing 3-7 at "+left+" "+right+" --> "+c37++);
					connectivity[left][i] = connectivity[i][left] = -1;
				}
			}
			for (int i=0; i<atoms.length; i++){
				for (int j=0; j<atoms.length; j++){
					if (i!=right && connectivity[left][i] != -1 && connectivity[left][i] != right){
						if (j!=left && connectivity[right][j] != -1 && connectivity[right][j] != left){
							connectivity[i][j] = connectivity[j][i] = -1;
						}
					}
				}
			}
		}
	}

	/**
	 * Update the `distance` matrix according to the `adjacency` matrix currently.
	 * This is done via Seidel's algorithm: https://en.wikipedia.org/wiki/Seidel%27s_algorithm
	 */
	private void updateDistanceMatrix() {
		int[][] newDistanceMatrix = floydWarshall(Util.adjacencyToWeight(adjacency));
		System.arraycopy(newDistanceMatrix, 0, distance, 0, atoms.length);
	}

	/**
	 * Floyd-Warshall algorithm for all-pairs shortest paths
	 */
	private int[][] floydWarshall(int[][] graphWeightArray) {
		int n = graphWeightArray.length;
		int dist[][] = new int[n][n];
		int i, j, k;

		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				dist[i][j] = graphWeightArray[i][j];
			}
		}

		for (k = 0; k < n; k++) {
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					if (dist[i][k] + dist[k][j] < dist[i][j]) {
						dist[i][j] = dist[i][k] + dist[k][j];
					}
				}
			}
		}
		return dist;
	}

	/**
	 * Implementation of Seidel's algorithm: APD returns the distance matrix of the graph with adjacency matrix A.
	 * This will only work on connected graphs.
	 * @param adjacencyArray adjacency matrix
	 * @param n dimension of matrix
	 * @return
	 */
	private int[][] apd(int[][] adjacencyArray, int n) {
		// base case: check if it is the complete graph, in which case return the adjacency matrix itself
		// (for complete graphs, adjacency matrix is its own distance matrix)
		boolean isCompleteGraph = true;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < i; j++) {
				if (adjacencyArray[i][j] == 0) {
					// if any off-diagonal entry is 0, then the graph is not the complete graph.
					isCompleteGraph = false;
					break;
				}
			}
		}
		if (isCompleteGraph) {
			return adjacencyArray;
		}

		SimpleMatrix A = new SimpleMatrix(Util.intToDouble2DArray(adjacencyArray));
		SimpleMatrix Z = A.mult(A);

		// Initialised to 0. Note that both A and Z are symmetric:
		// A is symmetric as it is the adjacency matrix of the graph
		// Z is symmetric as it is an adjacency matrix raised to a power.
		// So we only need to loop through off diagonal entries to find where there should be 1's.
		int[][] BArray = new int[n][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if ((adjacencyArray[i][j] == 1) || (Z.get(i, j) > 0)) {
					BArray[i][j] = 1;
				}
			}
		}

		int[][] TArr = apd(BArray, n);

		SimpleMatrix T = new SimpleMatrix(Util.intToDouble2DArray(TArr));
		SimpleMatrix X = T.mult(A);

		int[] degree = new int[n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				degree[i] += adjacencyArray[i][j];
			}
		}

		int[][] DArr = new int[n][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (X.get(i, j) >= (TArr[i][j] * degree[j])) {
					DArr[i][j] = 2*TArr[i][j];
				} else {
					DArr[i][j] = 2*TArr[i][j] - 1;
				}
			}
		}

		return DArr;
	}

	private boolean beforeLoop(int left, int right) {
		if (triangleWithDouble(left, right)) return true;
		if (tripleInLoop(left, right)) return true;
		if (twoDoubleInLoop(left, right)) return true;
		if (squareWithTwoDouble(left, right)) return true;
		return false;
	}

	private boolean triangleWithDouble(int left, int right) {
		for (int j=0; j<atoms.length; j++){
			int triangle = adjacency[left][j] * adjacency[right][j];
			if (triangle>=1){				
				if (adjacency[left][right]>=1 || triangle>= 2) return true;	// triangle with double
				for (int jj=0; jj<atoms.length; jj++){		// detect two connected triangles 
					if (jj == left || jj == right || jj == j) continue;
					if (adjacency[left][jj] * adjacency[right][jj]>0) return true;
					if (adjacency[j][jj] * adjacency[right][jj]>0) return true;
					if (adjacency[left][jj] * adjacency[j][jj]>0) return true;
				}
			}
		}
		return false;
	}

	private boolean twoDoubleInLoop(int left, int right) {
		int singleInLoopLeft = 0, singleInLoopRight = 0;
		if (adjacency[left][right]==1){		// it will become a double bond
			if (loopPart[left][right]==0) return false;	
			for (int i=0; i<atoms.length; i++) {
				if (i == left || i == right) continue;
				if (loopPart[left][i]>0){
					if (adjacency[i][left]  > 1) return true;	// if two consecutive doubles in a loop or different loops
					singleInLoopLeft++;
				}
				if (loopPart[right][i]>0) {
					if (adjacency[i][right] > 1) return true;	// if two consecutive doubles in a loop or different loops
					singleInLoopRight++;
				}
			}
		}
		return (singleInLoopLeft > 2 || singleInLoopRight > 2);	// a node with a double and at least three singles in loops
	}

	private boolean tripleInLoop(int left, int right) {
		if (adjacency[left][right]==2){		// it will become a triple bond
			return loopPart[left][right]>0;
		}
		return false;
	}
	
	private boolean squareWithTwoDouble(int left, int right){
		if (adjacency[left][right]==1){		// it will become a double bond
			for (int i=0; i<atoms.length; i++){
				if (i == right || i == left) continue;
				for (int j=0; j<atoms.length; j++){
					if (j==left || j==right) continue;
					if (adjacency[i][left] * adjacency[j][right] * adjacency[i][j] > 1)
						return true;
				}
			}
		}
		return false;
	}

	/**
	 * @param atom atom number, as an index to the atoms array / adjacency matrix row
     * Outputs if the atom with index `atom` has free valencies or not.
	 */
	private boolean isFull(final int atom) {
		int bondSum = 0;
		for (int i=0;i<atoms.length;i++)
			bondSum += adjacency[atom][i];
		return (atoms[atom].maxValence <= bondSum);
	}
	
	//---------------------------------------------------------------
	class SortCompare {
		int maxRow = atoms.length - 1;
		boolean equal;	
		int[] iPerm;
		
		SortCompare(){}
		SortCompare (int max){
			maxRow = max;
		}
		
		boolean isSmallerThan (int currentRow, int base) {
			equal = false;
			int []tPerm = iPerm.clone();
			for (int xCol=base+1; xCol<atoms.length; xCol++){
				int maxIndex=-1;
				int maxVal = -1;
				for (int yCol=atoms.length-1; yCol>=0; yCol--) {
					if (maxVal <= adjacency[currentRow][yCol] && tPerm[yCol] == -1 && 
//							compatible[yCol][xCol][base]) { 
							 !atoms[xCol].flag && !atoms[yCol].flag && atoms[xCol].symbol.equals(atoms[yCol].symbol) && areCompatible(yCol, xCol, base)) {
						maxIndex = yCol;
						maxVal = adjacency[currentRow][yCol];
					}
				}
//				if (maxIndex == -1) return false;	must not happen
				if (adjacency[base][xCol] < maxVal) return true;
				if (adjacency[base][xCol] > maxVal) return false;
				tPerm[maxIndex] = xCol;	// in order to avoid using the same value twice
			}
			equal = true;
			return false;
		}
		
		boolean isMinimal(){
			iPerm = new int[atoms.length];
			Arrays.fill(iPerm, -1);
			return sortCheckMin2(0);
		}
		
		private boolean sortCheckMin2(int base) {
			if (base == maxRow+1 || atoms[base].flag) return true;
			for (int row=0; row<=maxRow; row++){	// even check row == base
				// check what happens if base <- row ?
//				if (iPerm[row] != -1 || !compatible[row][base][base]) continue;
				if (iPerm[row] != -1 || !atoms[base].symbol.equals(atoms[row].symbol) || atoms[row].flag) continue;
				if (!areCompatible(row, base, base)) continue;
				iPerm[row] = base;
//				update(row, base);
				if (isSmallerThan(row, base)) return false;
				if (equal && !sortCheckMin2(base+1)) return false;
				iPerm[row] = -1;
			}
			return true;
		}
		
//		public void update(int row, int base) {
//			for (int src=0; src<atoms.length; src++)
//				for (int dst=0; dst<atoms.length; dst++)
//					compatible[src][dst][base+1] = compatible[src][dst][base] && adjacency[row][src] == adjacency[base][dst];
////					compatible[src][dst][base+1] = areCompatible(src, dst, base+1);
//		}
//		
		private boolean areCompatible(int src, int dst, int base) {
			for (int above=0; above<atoms.length; above++) {
				if (iPerm[above] != -1 && iPerm[above]!=base && adjacency[above][src] != adjacency[iPerm[above]][dst]) return false;
			}
			return true;
		}
	}

	//----------------------------------------------------------------------------------------

	private void outputMatrix(PrintStream out, int[][] fragment) {
		StringWriter writer = new StringWriter();
		writer.write(PMG.formula+"\n");
		writer.write("-----------"+(PMG.molCounter.get()-duplicate.get())+"\n");
		for (int i=0; i<atoms.length; i++) {
			for (int j=0; j<atoms.length; j++){
				writer.write(""+fragment[i][j]);
			}
			writer.write("\n");
		}
//		try {
			out.print(writer.toString());
//		} catch (IOException e) {
//			System.err.println("Could not output the Matrix.");
//			e.printStackTrace();
//		}
	}

    /**
     * Multithreading methods
     */
	@Override
	public void run() {
		PMG.availThreads.incrementAndGet();
		dispatch();
		PMG.startedTasks.incrementAndGet();
		PMG.pendingTasks.decrementAndGet();
	}

	private void dispatch() {
		switch(method){
		case SEM_CAN:
		case MIN_CAN:
		case OPTIMAL:
			generateOrderly(null, null);
			break;
		case CAN_AUG:
			augment();
			break;
		}
	}

//	private void submitNewTask() {
//		if (PMG.availThreads.decrementAndGet()>=0) {
//			PMG.pendingTasks.incrementAndGet();
//			PMG.executor.execute(new MolProcessor(atoms, nH, maxOpenings, adjacency, graph, canString, blocks, startLeft, startRight, acontainer, method, hashMap, cdkCheck));
//			return;
//		}
//		PMG.availThreads.incrementAndGet();
//		dispatch();
//	}

	private void submitNewTask() {
		int avail = PMG.availThreads.get();
		while (avail>1) {
			if (PMG.availThreads.compareAndSet(avail, avail-1)) {
				PMG.pendingTasks.incrementAndGet();
				PMG.executor.execute(new MolProcessor(atoms, nH, maxOpenings, neighbourhood, restrictNeighbourhoods, adjacency, fragment, connectivity, loopPart, distance, graph, canString, blocks, startLeft, startRight, acontainer, method, hashMap, cdkCheck, checkBad, goodlistquery, badlistquery));
				return;
			}
			avail = PMG.availThreads.get();
		}
		dispatch();
	}
	
	/*
	 * This function should be called first with startLeft=0 and startRight=1 as parameters.
	 * Input:
	 *  visited: set of visited graphs nodes (may not necessarily be complete molecules)
	 *  pString: canonical String of current graph, used with the Bliss library.
	 */
	void generateOrderly (Set<String> visited, String pString) {
    	if (visited == null && (method & CAN_AUG)!=0) {
    		visited = new HashSet<String>();
    		pString = canString;
    	}

    	// If the new edge (bond) (startLeft, startRight) cannot be added (due to valency constraints e.g.)
		if (startRight == atoms.length || isFull(startLeft)) {
			if (startLeft == atoms.length-2) return;
			int right = startRight;
			int [] oldBlocks = null;
			if ((method & SEM_CAN) != 0) { // if method is semi-canonicity
				oldBlocks = new int [blocks.length];
				System.arraycopy(blocks, 0, oldBlocks, 0, blocks.length);
				updateBlocks(startLeft);
			}
			startLeft+=1;
			startRight=startLeft+1;
			generateOrderly(visited, pString);
			blocks = oldBlocks;	// is relevant only if SEM_CAN, does not hurt otherwise!
			startLeft-=1;
			startRight=right;
			return;
		}
		else if (incBond(startLeft, startRight)) {
			maxOpenings-=2;
			// molecule final check: output / add to counter if passes
			if (acceptableStructure(visited, pString)) {
				checkMolecule();
				if (maxOpenings>0) submitNewTask();	// if there are still open places for new bonds
			}
			decBond(startLeft, startRight);
			maxOpenings+=2;
		}
		startRight++;
		generateOrderly(visited, pString);
		startRight--;
	}

    /**
     * Determines if the structure is canonical, return whether it will be kept or not.
     * @param visited HashSet of visited molecules. Will only be used when Bliss library is also being used.
     * @param pString a canonical string of the molecule, as calculated by Bliss.
     *                Hence this will be only used when canonical augmentation options has been selected,
     *                as it will make use of the Bliss graph library.
     */
	private boolean acceptableStructure(Set<String> visited, String pString) {
		if ((method&MIN_CAN) != 0) return new SortCompare(startLeft).isMinimal();
		if ((method&CAN_AUG) != 0) {
			// canonize 
			int[] perm1 = graph.canonize(this);
			canString = molString(perm1);
			return (visited.add(canString) && (pString.equals("") || pString.equals(degrade(perm1)))); 
		}
		return true;
	}

	
	private void updateBlocks(int row) {
		blocks[row+2] = 0;	// new block
		for (int i=row+3; i<atoms.length; i++)
			if (adjacency[row][i] != adjacency[row][i-1])
				blocks[i] = 0;	// new block
	}

	private void checkMolecule() {
		if (isComplete() && isConnected()) {
			if (!(hashMap || frag) && !new SortCompare().isMinimal()) {
				duplicate.incrementAndGet();
			} else {
				finalProcess(null);
			}
		}	
	}

	private void finalProcess(String canString) {
		if (cdkCheck && !acceptedByCDK()) {
			PMG.rejectedByCDK.incrementAndGet();
		} else {
		    // accepted molecule: add to set of molecules to keep
			if (hashMap || frag) {
				if (canString == null) canString = molString(graph.canonize(this));
				if (!molSet.add(canString)) {
					duplicate.incrementAndGet();
					return;
				}
			}
			long currentCount = PMG.molCounter.incrementAndGet();
			if(PMG.wFile){
				BufferedWriter outFile = PMG.outFile;
				StringWriter writer = new StringWriter();
				int bondsCount = 0;

				for (Atom atom : atoms) {
					writer.write(String.format("    0.0000    0.0000    0.0000 %-3s 0  0  0  0  0  0%n", atom.symbol));
				}
				for (int i=0; i<atoms.length; i++)
					for (int j=0; j<i; j++) {
						int order = adjacency[i][j];
						if (order > 0){
							bondsCount++;
							writer.write(String.format("%3d%3d%3d%3d%3d%3d%n",i+1, j+1, order, 0, 0, 0));
						}
					}

				writer.write("M  END\n> <Id>\n"+currentCount+"\n\n> <can_string>\n"+canString+"\n\n$$$$\n");
				// For the sake of calculating the bondsCount, we write the first line last!
				String sdfStr = String.format("%n  PMG%n%n%3d%3d  0  0  0  0  0  0  0  0  0%n",atoms.length, bondsCount)
						+ writer.toString();
				//			try {
				//				outFile.write(sdfStr);
				//			} catch (IOException e) {
				//				System.err.println("Could not write molecule to output file.");
				//				e.printStackTrace();
				//			} 
				PMG.fileWriterExecutor.submit(new MoleculeWriter(sdfStr, outFile));
		}
		}
	}


	private boolean acceptedByCDK() {
		try {
			IAtomContainer acprotonate = (IAtomContainer) acontainer.clone();
			for (int r=0; r<atoms.length; r++)
				for (int c=r+1; c<atoms.length; c++) 
					if (adjacency[r][c] > 0) acprotonate.addBond(r, c, bondOrder(adjacency[r][c]));
			return  CDKUtil.acceptedByCDK(acprotonate, nH) && 
					CDKUtil.checkGood(acprotonate, goodlistquery) &&
					CDKUtil.checkBad(acprotonate, badlistquery);	// if things go wrong, we assume it is accepted! :P
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		} catch (CDKException e) {
			e.printStackTrace();
			return false;
		}
		return true;
	}


	private Order bondOrder(int i) {
		switch (i){
			case 1: return IBond.Order.SINGLE;
			case 2: return IBond.Order.DOUBLE;
			case 3: return IBond.Order.TRIPLE;
			case 4: return IBond.Order.QUADRUPLE;
		}
		return null;
	}

	private void augment(){
        /**
         * Perform canonical augmentation on canString.
         * Hence will only be used when the canonical augmentation option has been set to true.
         */
		String pString = canString;
		if (maxOpenings<=0) return;
    	Set<String> visited = new HashSet<String>();
		for (int left = 0; left < atoms.length-1; left++){
			for (int right = left+1; right < atoms.length; right++){
					if (!incBond(left, right)) continue;	
					// canonize 
					int[] perm1 = graph.canonize(this);
					String childString = molString(perm1);
					if (visited.add(childString)) {	
						if (pString.equals("") || pString.equals(degrade(perm1))){ 
							maxOpenings -= 2;
							if (isComplete() && isConnected()) {
								finalProcess(childString);
							}	
							canString = childString;
							submitNewTask();
							maxOpenings += 2;
						}
					}
					decBond(left,right);
			}
		}

	}

	private String degrade(int[] canPerm) {
		// Find the canonical graph
		int[][] adjacency = new int [atoms.length][atoms.length];
		for (int i=0; i<atoms.length; i++)
			for (int j=0; j<atoms.length; j++){
				adjacency[canPerm[i]][canPerm[j]] = this.adjacency[i][j] - (frag?this.fragment[i][j]:0);
			}
		// In canonical graph, find the last edge
		int rowLast=0,colLast=0;
		for (int i = 0; i < atoms.length; i++){
			for (int j = i+1; j < atoms.length; j++){
				if (adjacency[i][j] != 0) {
					rowLast = i;
					colLast = j;
				}
			}
		}
		// map the last edge to the original graph
		for (int i=0; i<atoms.length; i++) if (canPerm[i] == rowLast){rowLast = i; break;}
		for (int i=0; i<atoms.length; i++) if (canPerm[i] == colLast){colLast = i; break;}
		// remove the last edge
		this.adjacency[rowLast][colLast] --;
		this.adjacency[colLast][rowLast] --;
		// canonize and report the canonical string
		int[] perm1 = graph.canonize(this);
		String decString = molString(perm1);
		this.adjacency[rowLast][colLast] ++;
		this.adjacency[colLast][rowLast] ++;
		return decString;
	}
	
	public void useFragment(String fileName) {
		System.out.println("Using the first molecule in "+fileName+" as the starting fragment.");
		try {
			Scanner inFile = new Scanner(new FileInputStream(new File(fileName)));
			int [] map = Util.readFragment(inFile, fragment, atoms);
			if (map == null) {
				System.out.println("Could not initialize the fragment.");
				return;
			}
			frag = true;
			for (int i=0; i<atoms.length; i++)
				for (int j=i+1; j<atoms.length; j++){
					if (fragment[i][j]>0) {
						connectIfNotBad(i,j);
						maxOpenings -= 2*fragment[i][j];
					}
					adjacency[j][i] = adjacency[i][j] = fragment[i][j];
				}
//			outputMatrix(System.out, adjacency);
			for (int i=0; i<map.length; i++) {
				atoms[map[i]].flag = true;
				if (map[i] < blocks.length-1) 
					blocks[map[i]+1] = 0;
			}
		} catch (FileNotFoundException e) {
			System.err.println("The fragments file "+fileName+" could not be found.");
		}
	}
}
