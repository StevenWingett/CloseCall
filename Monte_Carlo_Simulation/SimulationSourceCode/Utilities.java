package uk.ac.babraham.anacondamontecarlosimulation;

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
