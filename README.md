Heredity
========

Scripts for computing the heredity detection metric on a data set.

Requirements: numpy, scikit-learn

Usage: python ./computemetric.py datafile.txt
Output: The detected number of heritable states 

Data file format:
-----------------

The data file should be a space-separated dense matrix of data with no header row. Each row should be an independent run of the system under a different selection pressure, and each column should  be a numerical or binary feature.          
