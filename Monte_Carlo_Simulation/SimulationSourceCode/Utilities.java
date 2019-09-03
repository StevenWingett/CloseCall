package uk.ac.babraham.anacondamontecarlosimulation;

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


import java.util.Random;

public class Utilities {
		
	public static String makeRandomString(){	
		String characters = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
		String randomString = "";		
		Random randomNumberGenerator = new Random();

		for (int i = 1; i <= 20; i++) {
			int randomInt = randomNumberGenerator.nextInt(characters.length()); // Gives range 0..(numbFeatures - 1)
			String randomCharacter = characters.substring(randomInt, randomInt + 1); //extends to the end of this string or up to endIndex - 1
			randomString = randomString + randomCharacter;
		}			
		return randomString;	
	}

}
