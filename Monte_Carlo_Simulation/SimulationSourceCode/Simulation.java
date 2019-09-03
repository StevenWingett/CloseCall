package uk.ac.babraham.anacondamontecarlosimulation;

import java.util.HashMap;
import java.util.ArrayList;

/**
###################################################################################
###################################################################################
##This file is Copyright (C) 2019, Steven Wingett (steven.wingett@babraham.ac.uk)##
##                                                                               ##
##                                                                               ##
##This file is part of CloseCall.                                                ##
##                                                                               ##
##CloseCall is free software: you can redistribute it and/or modify              ##
##it under the terms of the GNU General Public License as published by           ##
##the Free Software Foundation, either version 3 of the License, or              ##
##(at your option) any later version.                                            ##
##                                                                               ##
##CloseCall is distributed in the hope that it will be useful,                   ##
##but WITHOUT ANY WARRANTY; without even the implied warranty of                 ##
##MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                  ##
##GNU General Public License for more details.                                   ##
##                                                                               ##
##You should have received a copy of the GNU General Public License              ##
##along with CloseCall.  If not, see <http://www.gnu.org/licenses/>.             ##
###################################################################################
###################################################################################
*/


public class Simulation {

	public static void main(String[] args) {
		
		System.out.println("Running Anaconda Monte Carlo Simulation");
			
		// ***Read arguments***
		String inputFilename = "";
		int simsToRun = 3;    //Default number of simulations
		Boolean qc = false;
		Boolean createRandomDataset = false;
		
		if(args.length > 0){
			inputFilename = args[0];
		}else{
			System.err.println("Please specify a file to process");
			System.exit(1);		
		}
		
		if(args.length > 1){
			simsToRun = Integer.parseInt(args[1]);
			System.out.println("Setting number of simulations to " + simsToRun);
		} else {
			System.out.println("Number of simulations defaulting to " + simsToRun);
		}
		
		if(args.length > 2){
			if(args[2].toLowerCase().equals("qc")){
				qc = true;
				System.out.println("Performing diagnostics QC on data");
			} else if(args[2].toLowerCase().equals("random")){
				createRandomDataset = true;
				qc = true;
			} else {
				System.err.println("Third argument not recognised, should be 'QC'.");
				System.exit(1);
			}
		}

		
		//  ***Read in data***	
		System.out.println("Reading in file " + inputFilename);
		
		Distribution features = new Distribution();
		Distribution featuresValGt1 = new Distribution();   //Features with valencies greater than 1
		Distribution valencies = new Distribution();
		
		HashMap<String, Double> observedInteractionsCounter = new HashMap<String, Double>();
		FileIO dataIO = new FileIO();	
		dataIO.inputData(inputFilename, features, valencies, observedInteractionsCounter, featuresValGt1, qc);
		
		int totalNumberFeatures = features.getSize();
		int totalNumberComplexes = valencies.getSize();
		
		System.out.println("Dataset comprised " + totalNumberFeatures + " features distributed over " + totalNumberComplexes + " complexes (including N=1)");

	
		//  ***Perform the simulations***
		String randomString = Utilities.makeRandomString();
		int[] notAddedPoolResults = new int[simsToRun];
		HashMap<String, Integer> compObsSimCounter = new HashMap<String, Integer>();  //Data structure to tally Obs Vs Sim
		HashMap<String, Integer> simCumulativeInteractions = new HashMap<String, Integer>();  //Records cumulated interactions
		for(String interaction : observedInteractionsCounter.keySet()){
			simCumulativeInteractions.put(interaction, 0);
			compObsSimCounter.put(interaction, 0);	
		}
		
		HashMap<String, Integer> simFeaturesCounter = new HashMap<String, Integer>();  //Record number of times feature observed in simulation
		for( String feature : features.getElements()){    //Initialize
			simFeaturesCounter.put(feature, 0);		
		}
		
		HashMap<String, Integer> simFeaturesValGt1Counter = new HashMap<String, Integer>();  //Record number of times multi-valent feature observed in simulation
		for( String feature : features.getElements()){    //Initialize
			simFeaturesValGt1Counter.put(feature, 0);	
		}
		
		
		for (int currentSimNumber = 1; currentSimNumber <= simsToRun; currentSimNumber++){	
			
			System.out.println("Simulation " + currentSimNumber);
			
			//Initialise data structures
			HashMap<String, Double> simInteractionsDecrementer = new HashMap<String, Double>();  //Records simulated interactions		
			for(String interaction : observedInteractionsCounter.keySet()){
				simInteractionsDecrementer.put(interaction, observedInteractionsCounter.get(interaction));
			}
			
			ArrayList<String> notAddedPool = new ArrayList<String>();    //To prevent biases arising from it not being possible to add a single feature multiple times to a given complex
			Complex simComplex = new Complex();
			ArrayList<String> randomDataset = new ArrayList<String>();  //For when createRandomDataset is 'true'			
				if(qc){   //Organisation this way may cause duplication of code, but minimises number of times qc is evaluated							
					for (int currentComplexNumber = 1; currentComplexNumber <= totalNumberComplexes; currentComplexNumber++){
						int valency = Integer.parseInt(valencies.getRandomElement());    //Set sim complex valency
						
						//System.out.println("Complex: " + currentComplexNumber + "\tValency: " + valency);  //Comment out as required
					
						if( (valency == 1)  && !notAddedPool.isEmpty() ){    //No interactions, but select from pool
							String poolFeatureToAdd = notAddedPool.get(0);
							simFeaturesCounter.put(poolFeatureToAdd, simFeaturesCounter.get(poolFeatureToAdd) + 1);
							notAddedPool.remove(0);	
							randomDataset.add(currentComplexNumber + "\t" + poolFeatureToAdd);
							//System.out.println(poolFeatureToAdd);  //Comment out as required
						
						}else if(valency == 1 ){    //No interactions, randomly select feature
							String featureToAdd = features.getRandomElement();
							simFeaturesCounter.put(featureToAdd, simFeaturesCounter.get(featureToAdd) + 1);
							randomDataset.add(currentComplexNumber + "\t" + featureToAdd);
							//System.out.println(featureToAdd);  //Comment out as required
						
						} else {    //Interaction needs 2 or more features
							int added = 0;
							if(!notAddedPool.isEmpty()){    //Add from the not added pool
								String poolFeatureToAdd = notAddedPool.get(0);
								simComplex.addFeature(poolFeatureToAdd);		
								randomDataset.add(currentComplexNumber + "\t" + poolFeatureToAdd);
								//System.out.println(poolFeatureToAdd);  //Comment out as required
								notAddedPool.remove(0);
								added++;		
								simFeaturesCounter.put(poolFeatureToAdd, simFeaturesCounter.get(poolFeatureToAdd) + 1);
								simFeaturesValGt1Counter.put(poolFeatureToAdd, simFeaturesValGt1Counter.get(poolFeatureToAdd) + 1);
							}
									
							do{
								String featureToAdd = features.getRandomElement();
								if(simComplex.addFeature(featureToAdd)){
									randomDataset.add(currentComplexNumber + "\t" + featureToAdd);
									//System.out.println(featureToAdd);  //Comment out as required
									added++;	
									simFeaturesCounter.put(featureToAdd, simFeaturesCounter.get(featureToAdd) + 1);
									simFeaturesValGt1Counter.put(featureToAdd, simFeaturesValGt1Counter.get(featureToAdd) + 1);
								} else {
									notAddedPool.add(featureToAdd);					
								}						
							}while(added < valency);				
						}			
						simComplex.recordSimInterationResults(simInteractionsDecrementer, simCumulativeInteractions, compObsSimCounter);	
						simComplex.empty();

					}
					
				} else {    //Not performing QC diagnostics
					
					for (int currentComplexNumber = 1; currentComplexNumber <= totalNumberComplexes; currentComplexNumber++){
											
						int valency = Integer.parseInt(valencies.getRandomElement());    //Set sim complex valency
						//System.out.println("Complex: " + currentComplexNumber + "\tValency: " + valency);  //Comment out as required
						
						if(valency == 1){
							if(!notAddedPool.isEmpty() ){    //No interactions, remove 1 item from pool
								notAddedPool.remove(0);
							}						
						} else {    //Interaction needs 2 or more features
							int added = 0;
							if(!notAddedPool.isEmpty()){    //Add from the not added pool
								String poolFeatureToAdd = notAddedPool.get(0);
								simComplex.addFeature(poolFeatureToAdd);
								//System.out.println(poolFeatureToAdd);  //Comment out as required
								notAddedPool.remove(0);
								added++;			
							}							
							do{
								String featureToAdd = features.getRandomElement();
								if(simComplex.addFeature(featureToAdd)){
									//System.out.println(featureToAdd);  //Comment out as required
									added++;					
								} else {
									notAddedPool.add(featureToAdd);					
								}						
							}while(added < valency);				
						}			
						simComplex.recordSimInterationResults(simInteractionsDecrementer, simCumulativeInteractions, compObsSimCounter);	
						simComplex.empty();	
					}					
				}
				notAddedPoolResults[currentSimNumber - 1] = notAddedPool.size();
				if(createRandomDataset){
					dataIO.createRandomDataset(inputFilename, currentSimNumber, randomString, randomDataset);
				}
			}	
		
			
		// ***Write out the results***

		System.out.println("Writing out results");
		dataIO.writeResults(inputFilename, simsToRun, observedInteractionsCounter, simCumulativeInteractions, compObsSimCounter, notAddedPoolResults, randomString);
		if(qc){
			dataIO.writeQCResults(inputFilename, simsToRun, features, featuresValGt1, simFeaturesCounter, simFeaturesValGt1Counter, randomString);
		}	
		System.out.println("Simulations completed");
	}
	
	
	

}
