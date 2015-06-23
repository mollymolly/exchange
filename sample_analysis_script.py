##################################################################
# This script uses ulfexchange to analyze a set of sample data.  #
##################################################################

from ulfexchange import  runBatch, runBatchRedAndGreen
import os

# Set up inputs
# The directory where the data is stored
directory = os.path.join(os.path.dirname(os.path.realpath(__file__)), \
            'sample data')
#directory = '/Users/mollyrossow/Desktop/10 sec photoconversion to analyse'
# File numbers of the files to be analyzed as strings
fileNumbers = ['002']
# The root name of each file. For example, one file in this set is control002red.tif. 
#rootName = '10 sec photoconversion'
rootName = 'control'

# Set the location and the radius of the converted region. This is determined 
# experimentally. 
(covertedCenter, convertedRadius) = ((147.0, 116.0), 50)

# Call the analysis routine
runBatch(fileNumbers, directory, rootName, covertedCenter, convertedRadius, \
         background = True)
        
#runBatchRedAndGreen(fileNumbers, directory, rootName, covertedCenter, convertedRadius)
