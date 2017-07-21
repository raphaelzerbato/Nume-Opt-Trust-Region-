# HMM
Hidden Markov Model
This program uses the Baum and Welsh algorithm to detect hidden regime through a vector of data. 

The files’ names starting with “hmm” are the central models. In each are mentioned the name of the densities used for each states (ie: hmm GBII/SkewNormal is the model using one generalised beta of type II density for one regime and a SkewNormal law of type II for another regime). For each model using either a SkewNormal law of type II or a Generalised Beta of type II (GBII), there are a few annexe files to download for each density. Each annexe file is a matlab function and the density it relates to is mentioned at the beginning of the name. Hence if you need to use the GBII/SkewNormal model, all the annexe files mentioning either GBII or SkewNormal needs to be downloaded. 
	
  Two annexes files needs to be downloading in addition of these necessary for each density mentioned earlier. These files are “mydoc” and “swtest” which respectively create a word document that sums up all the key estimations found for the model and is the Shapiro-Wilk Parametric test used to test if the estimated parameters match the data. 
