package uk.ac.babraham.anacondamontecarlosimulation;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.util.zip.*;
import java.util.HashMap;
import java.util.ArrayList;

public class FileIO {
	
	public void inputData(String filename, Distribution features, Distribution valencies, HashMap<String, Double> observedInteractionsCounter, Distribution featuresValGt1, boolean qc){
		
		//Decide whether input is zipped	
		FileInputStream fis = null;
		BufferedReader br = null;
	
		try {
			fis = new FileInputStream(filename);
			br = new BufferedReader(new InputStreamReader(new GZIPInputStream(fis)));
		} catch (IOException ioe) {
			try {
				if (fis != null) {
					fis.close();
				}
				br = new BufferedReader(new FileReader(filename));
			} catch (IOException ex) {
				ex.printStackTrace();
				System.exit(1);
			}
		}
		
		try {
			String line = null;
			br.readLine(); // Skip header
		
			int lineNumber = 0;
			int previousBarcodeID = 0;			
			Complex currentComplex = new Complex();

			while ((line = br.readLine()) != null) {
				lineNumber++;
				if(lineNumber % 1_000_000 == 0){
					System.out.println("Read in " + lineNumber + " lines");
				}	
						
				String[] lineElements = line.split("\t");
				int barcodeID = Integer.parseInt(lineElements[1]);
				String featureName = lineElements[7];
							
				features.addElement(featureName);
				
				if(barcodeID < previousBarcodeID){
					System.err.println("Barcodes in file are not in numerical order (i.e. barcode " + previousBarcodeID + " should be less than barcode " + barcodeID + ")");
					System.exit(1);
					
				}
				
				if ( (barcodeID != previousBarcodeID) && (previousBarcodeID != 0) ){    //New complex   			
					
					int valency = currentComplex.getValency();
					valencies.addElement(Integer.toString(valency));
							
					currentComplex.interactionsCounterIncrementer(observedInteractionsCounter);
					if(qc){    //Performing QC diagnostics
						currentComplex.addToDistribution(featuresValGt1, 2);
					}
					currentComplex.empty();
					currentComplex.addFeature(featureName);
					
									
				} else {
					if(!currentComplex.addFeature(featureName)){
						System.err.println("Impossible complex (contains repeat features) in datafile");
						System.exit(1);
					}
				}			
				previousBarcodeID = barcodeID;
			}	
			
			valencies.addElement(Integer.toString(currentComplex.getValency()));  //Add final valency at end of file
			if(qc){    //Performing QC diagnostics
				currentComplex.addToDistribution(featuresValGt1, 2);
			}
			currentComplex.interactionsCounterIncrementer(observedInteractionsCounter);
			currentComplex.empty();
			
			features.activateElements();    //Create proper data structures
			valencies.activateElements();
			if(qc){
				featuresValGt1.activateElements();
			}
						
		} catch (IOException ex) {
				ex.printStackTrace();
				System.exit(1);
		}		
	}
	

	
	public void writeResults (String inputFilename, int simsToRun, 
			HashMap<String, Double> observedInteractionsCounter,  
			HashMap <String, Integer> simCumulativeInteractions, 
			HashMap<String, Integer> compObsSimCounter, int[] notAddedPoolResults, 
			String randomString){
		
		String resOutFilename = inputFilename + ".MonteCarloResults." + randomString + ".txt.gz";
			
		//Write out the simulation results file
		try {
			BufferedOutputStream resultsFileOutStream = new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(resOutFilename), 2048));
			resultsFileOutStream.write("Name_Feature1\tName_Feature2\tObserved_Frequency\tSimulation_Average_Frequency\tObserved/Simulation\tSimulation_Score\tP_Value\n".getBytes());
			
			for(String interaction : observedInteractionsCounter.keySet() ){
				Double observed = observedInteractionsCounter.get(interaction);
				float simAvFreq = (float)simCumulativeInteractions.get(interaction) / simsToRun;
				Double obsSim = observed / simAvFreq;
				int simScore = compObsSimCounter.get(interaction);
				float pVal = (float)simScore / simsToRun;
				String lineToPrint = interaction + "\t" + Double.toString(observed) + "\t" + Double.toString(simAvFreq) + "\t";
				lineToPrint = lineToPrint + Double.toString(obsSim) + "\t" + Double.toString(simScore) + "\t";
				lineToPrint = lineToPrint + Float.toString(pVal) + "\n";
				resultsFileOutStream.write(lineToPrint.getBytes());			
			}		
			resultsFileOutStream.close();		
		} catch (IOException ioe) {
			ioe.printStackTrace();
			System.exit(1);
		}

	
		//Write out the un-emptied "unable to allocate features" file tally
		String poolOutFilename = inputFilename + ".MonteCarloUnallocatedPools." + randomString + ".txt.gz";
		try {
			BufferedOutputStream poolFileOutStream = new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(poolOutFilename), 2048));
			for(int unalloacted : notAddedPoolResults){
				String lineToPrint = Integer.toString(unalloacted) + "\n";
				poolFileOutStream.write(lineToPrint.getBytes());	
			}
			poolFileOutStream.close();
		} catch (IOException ioe) {
			ioe.printStackTrace();
			System.exit(1);
		}
	}
	
	
	
	public void writeQCResults(String inputFilename, int simsToRun, Distribution features, 
			Distribution featuresValGt1, HashMap<String, Integer> simFeaturesCounter, 
			HashMap<String, Integer> simFeaturesValGt1Counter, String randomString){
		
		//Write out the features count results
		String featuresCountOutFilename = inputFilename + ".MonteCarloFeaturesCount." + randomString + ".txt.gz";	
		try {
			BufferedOutputStream FeaturesCountFileOutStream = new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(featuresCountOutFilename), 2048));
			FeaturesCountFileOutStream.write("Feature\tObserved_Frequency\tSimulated_Average_Frequency\n".getBytes());	 
				
			HashMap <String, Integer> observedFeaturesCounter = features.getElementsCounter();					
			for( String feature : observedFeaturesCounter.keySet()){
				int obsCount = observedFeaturesCounter.get(feature);
				float simCount = (float)simFeaturesCounter.get(feature) / simsToRun;
				String lineToPrint = feature + "\t" + obsCount + "\t" + simCount + "\n";
				FeaturesCountFileOutStream.write(lineToPrint.getBytes());	
			}	
			FeaturesCountFileOutStream.close();
			
		} catch (IOException ioe) {
			ioe.printStackTrace();
			System.exit(1);
		}	
		
		//Write out the features count results (valency > 1)
		String featuresValGt1CountOutFilename = inputFilename + ".MonteCarloFeaturesCountValGt1." + randomString + ".txt.gz";	
		try {
			BufferedOutputStream FeaturesValGt1CountFileOutStream = new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(featuresValGt1CountOutFilename), 2048));
			FeaturesValGt1CountFileOutStream.write("Feature\tObserved_Frequency\tSimulated_Average_Frequency\n".getBytes());	 
				
			HashMap <String, Integer> observedFeaturesValGt1Counter = featuresValGt1.getElementsCounter();					
			for( String feature : observedFeaturesValGt1Counter.keySet()){
				int obsValGt1Count = observedFeaturesValGt1Counter.get(feature);
				float simValGt1Count = (float)simFeaturesValGt1Counter.get(feature) / simsToRun;
				String lineToPrint = feature + "\t" + obsValGt1Count + "\t" + simValGt1Count + "\n";
				FeaturesValGt1CountFileOutStream.write(lineToPrint.getBytes());	
			}	
			FeaturesValGt1CountFileOutStream.close();
			
		} catch (IOException ioe) {
			ioe.printStackTrace();
			System.exit(1);
		}	
	}
	
	
	
	public void createRandomDataset(String inputFilename, int currentSimNumber, String randomString, ArrayList<String> randomDataset){
		String randomDatasetFileName = inputFilename + ".RandomDataset." + currentSimNumber + "." + randomString + ".txt.gz"; 	
		
		try {
			BufferedOutputStream RandomDatasetOutStream = new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(randomDatasetFileName), 2048));
			RandomDatasetOutStream.write("Read_ID\tBarcode_ID\tChromosome	Feature_Start\tFeature_End\tFeature_Strand\tFeature_ID\tFeature_Name\n".getBytes());  //Header	 
			
			System.out.println("Writing random dataset " + currentSimNumber);
			for(int i = 0; i < randomDataset.size(); i++){
				String[] dataElements = randomDataset.get(i).split("\t");
				String barcodeID = dataElements[0];   //Will be an integer, but no point parsing to int
				String feature = dataElements[1];
				int readID = i + 1;
				String lineToPrint = readID + "\t" + barcodeID + "\t0\t0\t0\t0\t0\t" + feature + "\n";    //Allocate zeros to data that is not necessary nor available
				RandomDatasetOutStream.write(lineToPrint.getBytes());			
			}
			RandomDatasetOutStream.close();
			
		} catch (IOException ioe) {
			ioe.printStackTrace();
			System.exit(1);
		}			
	}
	
}
