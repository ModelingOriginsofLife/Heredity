Heredity Metric
===============

Scripts for computing the heredity detection metric on a data set.

Requirements: numpy, scikit-learn

Usage: python ./computemetric.py datafile.txt --<parameters>

Output: The detected number of heritable states 

Parameters:
-----------

--convergence FILENAME: Output convergence curve to a file

--convstep STEPSIZE: Resolution of the convergence curve - how many
samples to add each cycle
			  
--trials NUMTRIALS: When computing convergence curve, how many trials should 
be used to estimate variance

--eigvals FILENAME: Output eigenvalues to a file
			    
--eigvecs FILENAME: Output eigenvectors corresponding to heritable
states to a file
			  
--regularization COPIES: Sets number of regularization duplicates to
stabilize the algorithm in cases where there's no heredity

--alpha ALPHA: Gap sensitivity parameter. Higher values take more data to 
converge but are more robust.

--zscore: Rescale features to zero mean and unit standard deviation
			
Data file format:
-----------------

The data file should be a space-separated dense matrix of data with no header row. Each row should be an independent run of the system under a different selection pressure, and each column should  be a numerical or binary feature.          
