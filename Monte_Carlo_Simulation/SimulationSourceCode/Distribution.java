package uk.ac.babraham.anacondamontecarlosimulation;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import java.util.HashSet;


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


/**
 * 
 * @author wingetts
 * 
 * This class describes a distribution. Data may be added sequentially, which may then be stored in an array.
 * Elements may then be randomly selected from the array (the item is not removed after selecting
 */

public class Distribution {
	
	ArrayList<String> tempElements = new ArrayList<String>();
	String[] elements;
	int totalElements;
	Random rand = new Random();
	
	
	public void addElement(String element){  //Adds an element to the temporary array
		tempElements.add(element);
	}
	
	
	public void activateElements() {    //Create the array from which random features may be selected 
		totalElements = tempElements.size();
		elements = new String[totalElements];   //Set array size			
		for(int i = 0; i < totalElements; i++){
			elements[i] = tempElements.get(i);
		}
		tempElements.clear();    //Empty arraylist
	}
	
	
	public String getRandomElement(){
		return elements[rand.nextInt(totalElements)];
	}
	
	
	public String[] getElements(){    //Return all the elements as an array
		return elements;
	}
	
	public String[] getUniqueElements(){    //Return a unique list of elements as an array
		HashSet<String> uniqueElements = new HashSet<String>();
		for(String element: elements){
			uniqueElements.add(element);
		}
		
		String[] uniqueElementsArray = uniqueElements.toArray( new String[uniqueElements.size()] );
		return uniqueElementsArray;
	}
	
	
	public HashMap <String, Integer> getElementsCounter(){    //Return an HashMap tally of the elements->count
		HashMap <String, Integer> elementsCounter = new HashMap <String, Integer>();
		
		for(String element : elements){
			if( elementsCounter.containsKey(element) ){
				elementsCounter.put(element, elementsCounter.get(element) + 1);
			} else {
				elementsCounter.put(element, 1);
			}
		}		
		return elementsCounter;	
	}
	
	
	
	public void listsElements(){
		for(String element : elements){
			System.out.println(element);
		}	
	}
	
	
	public int getSize(){
		return totalElements;
	}
	
}
