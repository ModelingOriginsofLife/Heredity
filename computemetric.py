import numpy as np
import sys
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

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

# This is the sensitivity parameter alpha. If this is set to be very
# small, the algorithm will be more sensitive but will be more likely
# to blow up and give an anomalously large number of heritable states as its 
# output.
ALPHA = 0.1

def PCAmetric(population):
	count = 0
	stdev = 0
	fraction = 1
	
	subdata = population
	# Uncomment this to normalize each feature by its standard deviation
	# subdata = StandardScaler().fit_transform(population)

	pca = PCA().fit(subdata)

	latent = pca.explained_variance_ratio_
	
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
		
		threshold = (latent[0] - postmean)*ALPHA + postmean
	return count

# This function adds 10 copies of the first data row to the data set in
# order to regularize the algorithm against the case where there is only
# one heritable state

def doRegularization(data):
	regulator = data[0]
	newdata = np.copy(data)
	for i in range(10):
		newdata = np.vstack( (newdata, regulator) )
	
	return newdata

data = np.loadtxt(sys.argv[1])

print "{}".format(PCAmetric(doRegularization(data)))