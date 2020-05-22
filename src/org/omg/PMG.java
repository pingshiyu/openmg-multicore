package org.omg;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import org.omg.tools.Atom;
import org.omg.tools.Util;
import org.openscience.cdk.exception.CDKException;

import fi.tkk.ics.jbliss.Graph;

public class PMG{
	public static final int INFINITY = 99999;

	/**Output File containing the list of graph. */
	static BufferedWriter outFile;
	static AtomicInteger availThreads;
	final static AtomicLong molCounter = new AtomicLong(0);
	final static AtomicLong rejectedByCDK = new AtomicLong(0);
	final static AtomicLong pendingTasks = new AtomicLong(0);
	final static AtomicLong startedTasks = new AtomicLong(0);
	final static LinkedBlockingQueue<Runnable> taskQueue = new LinkedBlockingQueue<Runnable>();
	static ThreadPoolExecutor executor;
	static ExecutorService fileWriterExecutor;
	static int executorCount = 0;
	static boolean wFile = false;
	// Generate isomers of formula in the set. Set will be singleton if we are not generating a chemical space.
	static boolean chemicalSpace = false;
	static Set<String> formulae = new HashSet<>();
	static int weightLimit = 50;
	static int formulaeCounter = 0;
	static List<Character> chemSpaceAtoms = new ArrayList<>(Atom.atoms);
	static String currentFormula = null;
	static int[] neighbourhood = null;
	static boolean restrictNeighbourhoods = false;
	static boolean automaticCrowdingBound = false;
	static int method = MolProcessor.OPTIMAL;
	static boolean hashMap = false;
	static boolean cdk = true;
	static boolean checkBad = false;
	private static boolean verbose = false;
	private static String goodlist = null; //prescribed fragments to use as filter after generation process
    private static String badlist = null; //badlist of fragments to use as filter after generation process. Molecules should not contain them 
    private static String fragFile = null;
    private static MolProcessor mp;
    
	public static void main(String[] args) throws IOException, CDKException {
		if (args.length == 0) {
			usage();
			System.exit(0);
		}
		String out = parseArguments(args);

		startupMessages();

		long before = System.currentTimeMillis();

		if (wFile) {
			outFile = new BufferedWriter(new FileWriter(out));
			fileWriterExecutor = Executors.newSingleThreadExecutor();	// use an independent thread to write to file
		}

		for (String formula : formulae) {
			if (executorCount < 1) executorCount = Runtime.getRuntime().availableProcessors();	// default: use all cores available
			availThreads = new AtomicInteger(executorCount-1);	// number of tasks waiting in the queue; or, to generate in parallel

			formulaeCounter++;
			currentFormula = formula;
			if (automaticCrowdingBound) {
				neighbourhood = Util.formulaHRatioCrowdingBounds(formula, 3, true);
				//System.out.println("Bound used: " + Arrays.toString(neighbourhood));
			}

			// Parse atomic formula, exit if failed to parse.
			ArrayList<String> atomSymbols = Util.parseFormula(currentFormula);
			if (atomSymbols == null) System.exit(1);
			mp = new MolProcessor(atomSymbols, formula, neighbourhood, restrictNeighbourhoods, method, (method==MolProcessor.SEM_CAN) && hashMap, cdk, checkBad, fragFile != null, goodlist, badlist);
			if (fragFile != null) mp.useFragment(fragFile);
			executor = new ThreadPoolExecutor(executorCount, executorCount, 0L, TimeUnit.MILLISECONDS, taskQueue);
			pendingTasks.getAndIncrement();
			executor.execute(mp);	// start the generation
			wait2Finish();			// wait for all tasks to finish (pendingTasks == 0)
			executor.shutdown();	// stop the executor --> kill the threads
		}

		try { // close the output file
			if (wFile) {
				fileWriterExecutor.shutdown();
				while (!fileWriterExecutor.awaitTermination(60, TimeUnit.SECONDS));
				outFile.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			System.err.println("Failed while waiting for file outputs to complete.");
			e.printStackTrace();
		}

		long after = System.currentTimeMillis();
		// Report the number of generated molecules
		report(after - before);
	}

	private static String parseArguments(String[] args){
		int i=0;
		String out=null;
		boolean warned = false;
		try {
			// Add molecular formula (First argument)
			formulae.add(args[0]); currentFormula = args[0];

			for(i = 1; i < args.length; i++){
				if(args[i].equals("-p")){
					executorCount = Integer.parseInt(args[++i]);
				}
				else if(args[i].startsWith("-p")){
					executorCount = Integer.parseInt(args[i].substring(2));
				}
				else if(args[i].equals("-m")){
					method = setMethod(Integer.parseInt(args[++i]));
				}
				else if(args[i].startsWith("-m")){
					method = setMethod(Integer.parseInt(args[i].substring(2)));
				}
				else if(args[i].equals("-hashmap")){
					hashMap = true;
				}
				else if(args[i].equals("-chemspaceatoms")) {
					// If this has been set, then -wt / -weight must also be set.
					chemicalSpace = true;
					String atomsString = args[++i];
					Character[] atoms =
							atomsString.chars().mapToObj(c -> (char)c).toArray(Character[]::new);
					chemSpaceAtoms = Arrays.asList(atoms);
				}
				else if(args[i].equals("-wt") || args[i].equals("-weight")) {
					// If weightLimit option has been chosen, then we override the formula with the set of formulas.
					chemicalSpace = true;
					weightLimit = Integer.parseInt(args[++i]);
					formulae = Util.formulaeUnder(weightLimit, chemSpaceAtoms);
					System.out.println("Generating " + chemSpaceAtoms.toString() + " molecules under weight " + weightLimit + ". " +
							"Overriding inputted formula.");
				}
				else if(args[i].equals("-nocdk")){
					cdk = false;
				}
				else if(args[i].equals("-filter")){
					checkBad = true;
				}
				else if(args[i].equals("-o")){
					out = args[++i];
					wFile = true;
				}
				else if(args[i].equals("-neighbourhood")) {
				    // Expected string: comma separated list of integers representing data for the neighbourhood constraint
                    // function. neighbourhood[i] == max. no of neighbours within distance i. neighbourhood[0] = 1.
				    String neighbourhoodStr = args[++i];
				    if (neighbourhoodStr.startsWith("auto")) {
				    	automaticCrowdingBound = true;
					} else {
						String[] neighbourhoodValues = neighbourhoodStr.split(",");
						neighbourhood = new int[neighbourhoodValues.length];
						for (int j = 0; j < neighbourhoodValues.length; j++)
							neighbourhood[j] = Integer.parseInt(neighbourhoodValues[j]);
					}
				    restrictNeighbourhoods = true;
                }
				else if(args[i].equals("-v")){
					verbose = true;
				}
				else if(args[i].equals("-fr")){
					if (!Graph.blissFound) {
						continue;
					}
					if (fragFile == null) {
						fragFile = args[++i];
					} else {
						System.out.println("Disregarding extra fragment: "+args[++i]+". Consider using them with -goodlist.");
					}				
				}
				else if(args[i].equals("-goodlist")){
					if (goodlist == null)
						goodlist = args[++i];
					else 
						System.out.println("Disregarding the extra goodlist: "+args[++i]);
				}
				else if(args[i].equals("-badlist")){
					if (badlist == null)
						badlist = args[++i];
					else 
						System.out.println("Disregarding the extra badlist: "+args[++i]);
				} else {
					if (!warned){
						System.out.println("Warning: Invalid arguments are ignored. Run the program with no arguments to see the usage.");
						warned = true;
					}
				}
			}
		}catch(ArrayIndexOutOfBoundsException ae){
			System.err.println("Insufficient parameters provided. Run the program with no arguments to see the usage.");
			System.exit(1);
		}catch(NumberFormatException ne){
			System.err.println("Invalid number given: "+args[i]);
			System.exit(1);
		}
		return out;
	}

	private static int setMethod(int m) {
		switch (m){
		case 0:	return MolProcessor.SEM_CAN;
		case 1:	return MolProcessor.MIN_CAN;
		case 2:	if (Graph.blissFound) return MolProcessor.CAN_AUG; 
				else {System.out.println("Canonical augmentation is not available without bliss."); System.exit(2);}
		}
		System.err.println("Invalid method is selected.");
		usage();
		return 0;	// Invalid method
	}

	private static void usage() {
		System.out.println("Usage: java -jar PMG.jar <formula> [options]\n");
		System.out.println("Providing a formula for the elemental compositoin is obligatory, e.g., C4H7NO3.");
		System.out.println("You can further specify the following options.");
		System.out.println("\t-filter \tFilter bad substructures in the molecular structure.");
		System.out.println("\t-fr \tA file containing one substructure to use as initial structure for generation");
		System.out.println("\t-goodlist\tA file containing required substructures of the molecule (checked in the end) - only active if cdk is used");
		System.out.println("\t-badlist\tA file containing forbidden substructures (checked in the end) - only active if cdk is used");
		System.out.println("\t-hashmap\tEnables using a hashmap with semi-canonicity instead of the minimizer");
		System.out.println("\t-nocdk \tDisables using CDK for removing unacceptable molecular structures in the end.");
		System.out.println("\t-p  \tThe number of parallel threads to use; by default will use all available cores");
		System.out.println("\t-o  \tThe name of the output file");
		System.out.println("\t-v  \tverbose mode");
		System.out.println("\t-m \tThe method to use: \n\t\t\t0=semi-canonicity; \n\t\t\t1=minimization; \n\t\t\t2=canonical-augmentation;");
		System.out.println("Note that if you don't specify a method (-m) then the mix of semi-canonicity and minimization will be used.");
	}

	private static void startupMessages() {
		if (!chemicalSpace) {
			System.out.print("Processing " + currentFormula + " using");
			switch (method) {
				case MolProcessor.CAN_AUG:
					System.out.println(" canonical augmentation with bliss as canonizer");
					break;
				case MolProcessor.MIN_CAN:
					System.out.println(" only minimization");
					break;
				case MolProcessor.SEM_CAN:
					System.out.println(" semi-canonization and " + (hashMap ? "hash map" : "minimization test on finished molecules"));
					break;
				case MolProcessor.OPTIMAL:
					System.out.println(" semi-canonization and minimization test on every generated structure");
					break;
			}
		}
		System.out.print("Selected options: ");
		if (restrictNeighbourhoods) System.out.print("neighbourhood restriction " + Arrays.toString(neighbourhood) + ", ");
		if (cdk) System.out.print("cdk"+(goodlist==null?"":" with good-list")+(badlist==null?"":(goodlist==null?" with":" and")+" bad-list")+", ");
		if (checkBad) System.out.print("bad-sub filter, ");
		if (wFile) System.out.print("with output to file, ");
		System.out.println((executorCount > 1)? ""+executorCount+" threads." : "sequential execution.");
		if (!cdk && (goodlist!= null || badlist != null)) System.out.println("Warning: The provided good- and/or bad-list will not be used because CDK is disabled.");
	}

	private static void report(long duration) {
		if (chemicalSpace)
			System.out.println("\r" + formulaeCounter + " formulae has been enumerated. Total molecules: " + molCounter.get());
		else
			System.out.println("\rFinal molecule count:  " + molCounter.get());

		if (verbose) {
			if (method == MolProcessor.SEM_CAN || mp.frag) System.out.println("Duplicates removed in the end: "+MolProcessor.duplicate.get());
			if (cdk) System.out.println("Rejecting "+rejectedByCDK.get()+" by CDK.");
			System.out.println("Started Tasks: "+startedTasks.get());
		}
		System.out.println("Duration: " + (duration) + " milliseconds");
	}

	private static void wait2Finish() {
		int waitTime = 5;
		int timer = waitTime;
		while (0 < pendingTasks.get()){
			try {
				Thread.sleep(waitTime);
				if (timer > 0)
					timer--;
				else {
					if (chemicalSpace)
						System.out.print("\rFormula #" + formulaeCounter + ": " + currentFormula + ", " +
							"molecules thus far: " + molCounter.get());
					else
						System.out.print("\rTotal molecules thus far: " + molCounter.get());
					timer = waitTime;
				}
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}

}
