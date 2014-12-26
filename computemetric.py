import numpy as np
import sys
import argparse
import math
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cross_validation import ShuffleSplit

#########################################################################
###
### Compute Heredity Metric
###
### Use: python ./computemetric.py datafile.txt
### Output: The detected number of heritable states 
### 
### Data file format:
###
### The data file should be a space-separated dense matrix of data
### with no header row. Each row should be an independent run of the    
### system under a different selection pressure, and each column should 
### be a numerical or binary feature.          
###
#########################################################################

globalargs = {}

def PCAmetric(population, finalpass = False):
	count = 0
	stdev = 0
	fraction = 1
	
	if globalargs.zscore:
		subdata = StandardScaler().fit_transform(population)
	else:
		subdata = population

	pca = PCA().fit(subdata)

	latent = pca.explained_variance_ratio_
	
	if (globalargs.eigvals != None) and finalpass:
		np.savetxt(globalargs.eigvals, latent, delimiter=' ')
	
	threshold = latent[latent.shape[0]-1]
	
	# Iteratively find the eigenvalue gap
	for iter in range(1000):
		lcount = 0
		for i in range(latent.shape[0]):
			if (latent[i] > threshold):
				lcount += 1
			
		count = lcount
		postmean = 0
		norm = 0
		for i in range(lcount, 2*lcount):
			if i<latent.shape[0]:
				postmean += latent[i]
				norm += 1
		
		postmean /= norm
		
		# This is the sensitivity parameter alpha. If this is set to be very
		# small, the algorithm will be more sensitive but will be more likely
		# to blow up and give an anomalously large number of heritable states as its 
		# output.
		
		threshold = (latent[0] - postmean)*globalargs.alpha + postmean
	
	if (globalargs.eigvecs != None) and (count>0) and finalpass:
		np.savetxt(globalargs.eigvecs, pca.components_[0:count,:], delimiter=' ')
		
	return count

# This function adds 10 copies of the first data row to the data set in
# order to regularize the algorithm against the case where there is only
# one heritable state

def doRegularization(data):
	regulator = data[0]
	newdata = np.copy(data)
	for i in range(globalargs.regularization):
		newdata = np.vstack( (newdata, regulator) )
	
	return newdata

def myargs():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description = 
""" 
Compute heredity metric and related quantities on a given set of data.
\nUsage:
\n
\n\tpython computemetric.py dataset.txt --<parameters>
\n
""")
	parser.add_argument('datafile', default = None, type=str,
                        help='input data file')
	parser.add_argument('--convergence', default = None, type = str, help = "Output convergence curve to a file (note, this will take significantly longer to compute!)")
	parser.add_argument('--convstep', default = 1, type = int, help = "Resolution of the convergence curve - how many samples to add each cycle")
	parser.add_argument('--trials', default = 1, type = int, help = "When computing convergence curve, how many trials should be used to estimate variance")
	parser.add_argument('--eigvals', default = None, type = str, help = "Output eigenvalues to a file")
	parser.add_argument('--eigvecs', default = None, type = str, help = "Output eigenvectors corresponding to heritable states to a file")
	parser.add_argument('--regularization', default = 10, type = int, help = "Sets number of regularization duplicates to stabilize no-heredity case")
	parser.add_argument('--alpha', default = 0.1, type = float, help = "Gap sensitivity parameter. Higher values take longer to converge but are more robust.")
	parser.add_argument('--zscore', action = "store_true", help = "Rescale features to zero mean and unit standard deviation")
  
	args = parser.parse_args()
        
	return args

globalargs = myargs()

if globalargs.datafile != None:
	data = np.loadtxt(globalargs.datafile)

	if globalargs.convergence != None:
		for count in range(1,data.shape[0], globalargs.convstep):
			detectedStates = 0.0
			detectedStatesStdev = 0.0
			shuffler = ShuffleSplit(data.shape[0], n_iter = globalargs.trials, test_size = count)
			
			for train_index, test_index in shuffler:
				result = PCAmetric(doRegularization(data[test_index]), False)
				detectedStates += result
				detectedStatesStdev += result*result
				
			outfile = open(globalargs.convergence,"a")
			if globalargs.trials == 1:
				outfile.write("{} {}\n".format(count, detectedStates))
			else:
				detectedStates /= globalargs.trials
				detectedStatesStdev /= globalargs.trials
				detectedStatesStdev = math.sqrt(detectedStatesStdev - detectedStates*detectedStates)
				outfile.write("{} {} {}\n".format(count, detectedStates, detectedStatesStdev))
			outfile.close()
			
	print "{}".format(PCAmetric(doRegularization(data), True))
