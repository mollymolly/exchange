from ulfexchange import  runBatch
import os


##########################
# This script runs the ULF intensity analysis routine on a set of data. 
##########################

# CHANGE HERE
# The directory where the data is stored
directory = os.path.dirname(os.path.realpath(__file__)) # directory of this file. 
# File numbers of the files to be analyzed. 
fileNumbers = ['002','004']
# The root name of each file. For example, a file one file in this set is control002red.tif. 
rootName = 'control'

# Set the location and the radius of the converted region. This is determined experimentally. 
(covertedCenter, convertedRadius) = ((147.0, 116.0), 50)
# Call the analysis routine
runBatch(fileNumbers, directory, rootName, covertedCenter, convertedRadius, background = True)
#runBatch(['01','02','03'], '/data', 'control_', (150,110), 50, background = True)