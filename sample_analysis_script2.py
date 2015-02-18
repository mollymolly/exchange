##################################################################
# This script uses ulfexchange to analyze a set of sample data.  #
##################################################################

from ulfexchange import  runBatchRedAndGreen
import os

# Set up inputs
# The directory where the data is stored
directory = os.path.join(os.path.dirname(os.path.realpath(__file__)), \
            'sample data')
# File numbers of the files to be analyzed as strings
fileNumbers = ['003','004','006']
# The root name of each file. For example, one file in this set is control002red.tif. 
rootName = 'control'

# Set the location and the radius of the converted region. This is determined 
# experimentally. 
(covertedCenter, convertedRadius) = ((147.0, 116.0), 50)

# Call the analysis routine
runBatchRedAndGreen(fileNumbers, directory, rootName, covertedCenter, convertedRadius)