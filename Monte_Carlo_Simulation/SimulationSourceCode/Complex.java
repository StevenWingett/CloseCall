package uk.ac.babraham.anacondamontecarlosimulation;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.HashMap;

/**
 * Describes the complex comprising multiple features
 * @author wingetts
 * 
 */

public class Complex {
	
	HashSet<String> featuresList = new HashSet<String>();    //For checking whether a features has been added
	ArrayList<String> featuresArray = new ArrayList<String>();    //For generating interactions
	
	
	//Add feature (evaluates whether possible to add feature) 
	public boolean addFeature(String feature){	
		if(featuresList.contains(feature)){
			return false;
		} else {
			featuresList.add(feature);	
			featuresArray.add(feature);
			return true;
		}
	}
	
	
	//Get valency
	public int getValency(){
		return featuresList.size();
	}
	
	
	//Empty the complex
	public void empty(){
		featuresList.clear();
		featuresArray.clear();
	}
	
	
	//Get pairwise interactions and increment an interactions counter as desired
	public void interactionsCounterIncrementer(HashMap <String, Double> interactionsCounter){
		int numberFeatures = featuresList.size();
		//int numberInteractions = ( (numberFeatures * numberFeatures) - numberFeatures ) / 2;
		//String[] interactions = new String[numberInteractions];
		//int k = 0;
		
	//	System.out.println("Number features:" + numberFeatures + "\tNumber interactions:" + numberInteractions);
		
		
		for(int i = 0; i < (numberFeatures - 1); i++){    //Proceed to penultimate feature
			for(int j = i + 1; j < numberFeatures; j++){    //Proceed to last feature	
				
				
			//	System.out.println(i + " " + j);
			//	System.out.println(featuresArray.get(i) + "\t"  + featuresArray.get(j) );
				String interaction;
				if( featuresArray.get(i).compareTo(featuresArray.get(j)) > 0 ) {
					interaction = featuresArray.get(j) + "\t" + featuresArray.get(i);
				} else {
					interaction = featuresArray.get(i) + "\t" + featuresArray.get(j);
				}
				
				if(interactionsCounter.containsKey(interaction)){    //Check this interaction already included in the counter
					interactionsCounter.put(interaction, interactionsCounter.get(interaction) + 1);
				} else {
					interactionsCounter.put(interaction, 1d);				
				}
				//interactions[k] = interaction;
				//k++;
			}		
		}
	
	}
	
	
	
	//Get pairwise interactions and increment an interactions counter as desired
	public void recordSimInterationResults(HashMap <String, Double> simInteractionsDecrementer, HashMap <String, Integer> simCumulativeInteractions, HashMap <String, Integer>compObsSimCounter){
		int numberFeatures = featuresList.size();
			
		for(int i = 0; i < (numberFeatures - 1); i++){    //Proceed to penultimate feature
			for(int j = i + 1; j < numberFeatures; j++){    //Proceed to last feature	
					
				String interaction;
				if( featuresArray.get(i).compareTo(featuresArray.get(j)) > 0 ) {
					interaction = featuresArray.get(j) + "\t" + featuresArray.get(i);
				} else {
					interaction = featuresArray.get(i) + "\t" + featuresArray.get(j);
				}
					
				if(simInteractionsDecrementer.containsKey(interaction)){    //Check this interaction already included in the counter
					simInteractionsDecrementer.put(interaction, simInteractionsDecrementer.get(interaction) - 1);
					simCumulativeInteractions.put(interaction, simCumulativeInteractions.get(interaction) + 1);
					if(simInteractionsDecrementer.get(interaction) == 0){
						compObsSimCounter.put(interaction, compObsSimCounter.get(interaction) + 1);    //Record result interaction occurs as much in simulation as observed
						//simInteractionsDecrementer.remove(interaction);  //Don't do this - it should optimises code but produced biased graphs where the observed can never be more than the simulated
					}
				}
			} 
 
		}		
	}
	
	public void printFeatures(){    //Prints the complex features
		
		int numberFeatures = featuresArray.size();
		for(int i = 0; i < numberFeatures ; i++){
			System.out.println(featuresArray.get(i));			
		}	
	}
	
	
	public void addToDistribution(Distribution features, int min){    //Takes a distribution, and adds features if the complex valency meets threshold minimum
		if (this.getValency() >= min){
			for( int i = 0; i < featuresArray.size(); i++ ){
				features.addElement(featuresArray.get(i));			
			}	
		}	
	}
	
	
}
	
	


