# -*- coding: utf-8 -*-
""" Introduction
---------------------------------------
This module analyzes intensity changes in confocal images of cells expressing 
mEos-vimentin :superscript:`Y117L` as described in "Vimentin filament 
precursors exchange subunits in an ATP-dependent manner" Amélie Robert, Molly 
Rossow, Caroline Hookway, Stephen A. Adams and Vladimir I. Gelfand (Accepted to 
PNAS June 2015.)


The experiments described by Robert *et al.* use a photoconvertible fluorophore
(mEos3.2) to track the behavior of a mutant of the vimentin protein 
(vimentin :superscript:`Y117L`) that forms oligomers known as unit length 
filaments (ULFs). Fluorescently labeled ULFs are visible as bright particles 
in cells. The fluorophore mEos converts from green to red when exposed to UV 
light. To observe how mEos-vimentin :superscript:`Y117L` is incorporated into 
ULFs we photoconverted one region of the cell and observed fluorescence 
accumulate in the ULFs elsewhere in the cell. This analysis software measures 
the change in intensity in ULFs. 

Work Flow
----------------------------------------
Robert *et al.* used ulfexchange.py to analyze confocal images of confocal images
of cells expressing mEos-vimentin :superscript:`Y117L` with the following steps:

1. **Bleach correction**: Red channel images were bleach corrected by scaling 
   each image so that its mean was the same as the first image.
2. **ULF intensity calculation**: Coordinates of the center of each ULF in each 
   frame (obtained separately using the green channel) were used to measure the 
   intensity of each ULF in each red frame by calculating the average fluorescence 
   in a circle of fixed size centered at these coordinates.
3. **Background subtraction**: The intensity of each ULF before conversion was 
   subtracted from all its post-conversion intensities.
4. **Normalization between data sets**: ULF intensities were normalized by 
   dividing by the background subtracted intensity in the photoconverted region.
5. **Numerical and graphical results** Slope calculations, plotting, values 
   written file, etc.

The function analyzeDataSet implements the workflow described above. Further 
details of the analysis are included in Robert *et. al.* 

Required inputs
------------------------------------------
This analysis takes as input four files: 

1. A TIFF file containing a time series of red channel images of the 
   mEos-vimentin :superscript:`Y117L` cells after photo conversion.
2. A TIFF file with a times series green channel images corresponding to the 
   first tiff file. 
3. A TIFF file with a single red channel image of the same cell before photo 
   conversion for determining the background fluorescence. 
4. A text file containing the coordinates for the ULFs of interest in each 
   frame. For the work described in Robert *et al.* 
   `Diatrack <http://diatrack.org/>`_\ was used to obtain these coordinates. 


Naming convention for batch processing
--------------------------------------------
The function runBatch performs a complete 5-step analysis on multiple sets of 
data. Each data file described above needs to be named with the following 
convention: 

1. rootname + red + number + .tif
2. rootname + green + number + .tif
3. rootname + before + number + .tif
4. rootname + .txt 

Where rootname is a name that is the same for all files in the batch and number
is unique. 

For example if rootname is control the following files are an allowed data set: 

===== =============== ====================== 
#     Contents        Name
===== =============== ======================
1     Red channel     control_red01.tif
2     Green channel   control_green01.tif
3     Red background  control_before01.tif
4     Coordinates     control_01.tif
===== =============== ======================

Outputs
----------------------------------------------
The function runBatch generates a number of out put files. Each individual data 
set (comprised of the four files described above) generates the following five
files: 

1. An excel file containing the intensities of each ULF and it's distance from 
   the photoconverted region. 
2. A png file with the red and green channel of the first frame. The red channel 
   shows the region of photo conversion as a red circle. The green channel shows 
   the tracks of the ULFs. 
3. A png file with a plot of the slope the intensity in the converted region 
   versus time.
4. A png file with all intensities for all ULFs plotted versus time.
5. A png file with ULF intensities plotted versus time grouped by distance from 
   the converted region. 

In addition, two batch files are generated summarizing the results of all data
sets analyzed in the batch. 

1. An excel file allslopes.xlsx contain the slopes for all ULFs in the batch 
   groups by distance from the converted region by sheet of the spread sheet.
2. A png file all_intensities_by_distance.png that contains plots of all ULF 
   intensities grouped by distance from the converted region. 

The file exchangeExample.py contains and example of using runBatch.

Example
-----------------------------------------------
Below is an example script for calling runBatch on three sets of data::
        

    from ulfexchange import  runBatch
    import os
    
    # Set up inputs
    # The directory where the data is stored
    directory = os.path.join(os.path.dirname(os.path.realpath(__file__)), 
    'sample data')
    # File numbers of the files to be analyzed as strings
    fileNumbers = ['002','004','006']
    # The root name of each file. For example, one file in this set is 
    # control002red.tif. 
    rootName = 'control'
    
    # Set the location and the radius of the converted region. This is determined 
    # experimentally. 
    (covertedCenter, convertedRadius) = ((147.0, 116.0), 50)
    
    # Call the analysis routine
    runBatch(fileNumbers, directory, rootName, covertedCenter, convertedRadius, 
    background = True)

Requirements
-----------------------------------------------
* Python 2.7
* csv 1.0
* xlsxWriter 0.3.7-1
* tifffile 2013.11.03 (http://www.lfd.uci.edu/~gohlke)
* Matplotlib 1.3.1
* Numpy 1.8.1
* scipy 0.14.0



Functions
------------------------------------------------
"""


from subprocess import check_output
from os.path import dirname, realpath, basename, splitext, join
from heapq import nsmallest
from copy import deepcopy
import csv
import xlsxwriter
from matplotlib.pylab import  figure, cm, close, Circle, tight_layout, annotate
from tifffile import TiffFile, imsave
from numpy import shape, floor, ceil, swapaxes, sqrt, mean, zeros, amax
from scipy.stats import linregress


#############
# Data I/O  #
############
def readInTif(fileName):
    """Reads in a single frame or multi-frame tiff file into a numpy array. 

       Uses TiffFile written by Christoph Gohlke.
    
    Arguments:
       fileName: Full path to tiff file.
    Returns:
       Numpy array.
    """
    tif = TiffFile(fileName)
    images = tif.asarray()
    images = images.astype('float')
    tif.close()
    s = shape(images)
    # If the tif file contains a series of images,
    # make time the third axis of the matrix.
    if len(s) == 3:
        images = swapaxes(images,0,2)
    return images
    
def saveTif(fileName, data):
    """Saves a numpy array as a TIFF file.

       This function was used for trouble shooting experiments.

    Arguments:
        fileName: Full path to tiff file.
        data: Numpy arrays with the format [x,y,t].
    """
    data[data >= 2**16] =  2**16 - 1  # deal with 16 bit overflow
    data[data <= 0] = 0
    data = data.astype('uint16')
    s = shape(data)
    # tifffile assumes [t,x,y] so swaps axes.
    if len(s) > 2:
        data = swapaxes(data,0,2)
    imsave(fileName, data)
      
def readInCoordinatesFile(fileName, nFrames):
    """Reads in a text file of coordinates. 

       Reads in a text file of coordinates created by Diatrack particle tracking
       software. Excludes particles that don't exist in every frame. 
       
    Arguments:
        fileName: Full path to coordinate file.
        nFrame: Number of frames.
    Return: Returns a list of tracks with the format: 
        [([x11,x12,x13],[y11,y12,y13]),([x21,x22,x23],[y21,y22,y23])]
        where particle  1 has position (X11, y11)  in frame 1 and position 
        (x12, y12) in frame 2.
    """
    # open file
    with open(fileName, 'rU') as coordfile:
        tracksReader = csv.reader(coordfile, delimiter=' ', quotechar='|')
        fileContents =[]
        for row in tracksReader:
            fileContents.append(row)
    
    tracks = []
    titleRow  = fileContents[0]

    for r in fileContents[1:]: #skip first row which contains header information.
        # remove empty strings
        numberRow = []
        for c in r:
            if c != '' and c != '.':
                n = float(c)
                if n != 0:
                    numberRow.append(float(c))
        
        # Check if data is rowwise or columnwise and divide sort into x and y.
        # y is first. Really. It's in the note at the top of each .txt 
        if titleRow[1] == '(columnwise):':
            raise _fileFormatError(fileName)
            return tracks
        else:
            y = numberRow[1:-1:3] 
            x = numberRow[2:-1:3]

        # Exclude tracks that don't exist in every frame. 
        if (len(x) == nFrames) and (len(y) == nFrames):
            tracks.append((x,y))
        
    return tracks


################
# Calculations #
################
def distance(p1, p2):
    """Calculates the distance between two points.

       Arguments: 
           p1: Point (x1, y1).
           p2: Point (x2, y2).
    """
    x1 = p1[0]
    y1 = p1[1]
    x2 = p2[0]
    y2 = p2[1]
    d = sqrt((x1 - x2) ** 2 + (y1 - y2) **2) 
    return d
    
def getIntensities(mat, tracks, radius = 6):
    """Finds the average intensities in circles in each frame of mat. 

    Arguments:
        mat: 3D numpy array (x,y,time).
        tracks: Coordinates of ULF locations, created by readInCoordinatesFile 
             and has the format: 
             [([x11,x12,x13],[y11,y12,y13]),([x21,x22,x23],[y21,y22,y23])]
             where particle 1 has position (X11, y11)  in frame 1 and position 
             (x12, y12) in frame 2. 
        radius: Radius of the circle of the circle around each location. 
    Returns:
         List of lists of intensities. Each sublist corresponds to a location 
         in tracks. Sublists are in the same order as the coordinates in tracks. 
     """
    s = shape(mat)
    if len(s) == 3:
        npoints = s[2]
    else:
        npoints = 1
    intensities = []
    for x,y in tracks:
        trace = []
        for t in range(npoints):
            center = (x[t], y[t])
            trace.append(getAverageIntensityInCircle(mat[:,:,t], center,\
            radius))
        intensities.append(trace)
    return intensities
    
def getAverageIntensityInCircleOverTime(imMat, center, radius):
    """Finds the average intensity in a circle over time. 

       Finds the average intensities in a single circle defined by 
       center and radius in each frame of imMat. This function uses the same 
       center for all frames.

    Arguments:
        inMat: 3D numpy array (x,y,time).
        center: (x, y) coordinate pair.
        radius: Circle radius.
    Returns:
        List of intensities in circle over time. 
    """
    
    s = shape(imMat)
    averageIntensities = []
    # For each frame, find the average intensity in the circle.
    if len(s) > 2:
        for t in range(s[2]):
            averageIntensities.append(getAverageIntensityInCircle(imMat[:,:,t],
            center,radius))
    else:
        averageIntensities.append(getAverageIntensityInCircle(imMat[:,:],
        center,radius))
    return averageIntensities
    
def getAverageIntensityInCircle(frame, center, radius):
    """Finds the average intensity in a circle. 

       Returns the average intensity of pixels in a circle defined  by center 
       and radius in the 2D array frame.
       
      Arguments:
         frame: 2D numpy array (x,y).
         center: (x,y) coordinate pair.
         radius: Circle radius.
     Returns: 
         Average intensity in the circle. 
    """
    
    s = shape(frame)
    includedIntensities = [] # the points in the circle
    centerX = center[0]
    centerY = center[1]
    
    # Define a box around the circle. 
    left = int(floor(centerX - radius))
    right = int(ceil(centerX + radius))
    bottom = int(floor(centerY- radius))
    top = int(ceil(centerY + radius))

    # Check all points in this box to see if they are in the circle. 
    for r in range(bottom, top + 1):
        for c in range(left, right + 1):
            d = distance((c,r), center) # Find distance between point and center
                                        # coordinates are in (x,y) pairs 
                                        # x->c, y->r
            # If the point is in the circle add it to included intensities
            if d <= radius: 
                # Check if the coordinates are within the frame. 
                if ((r < s[0]) and (r >= 0)) and ((c < s[1]) and (c >= 0)):
                    #Access matrix with rowsXcolums
                    includedIntensities.append(frame[r,c]) 
    #Return the mean of the points in the circle  
    return mean(includedIntensities)

def getDistancesToFirstPointInTrack(tracks, point):
    """Finds the distance between a point and the first point in tracks.

    Arguments:
        tracks:  A list returned by readInCoordinatesFile. tracks has the format
                 [([x11,x12,x13],[y11,y12,y13]),([x21,x22,x23],[y21,y22,y23])]
                  where particle 1 has position (X11, y11)  in frame 1 and 
                  position (x12, y12) in frame 2.
       point: (x,y) coordinate pair.
    Returns:
       The distance between the first point in tracks and point. 
    """
    text_file = open("distances002.txt", "w")
    distances = []
    for t in tracks:
        xs = t[0]
        ys = t[1] 
        d = distance((xs[0],ys[0]), point)
        distances.append(d)
        text_file.write("y = {0:.2f} x = {1:.2f} d = {2:.2f}\n".format(xs[0], ys[0], d))
    text_file.close()   
    return distances

def subtractULFBackground(ULFBackgrounds, intensities):
    """ Subtracts the background intensity from ULF intensity. 

    Subtracts the background intensity of a ULF from the same ULF in each
    frame after photoconversion. ULFBackgrounds is a list of background 
    intensities. intensities is a list of ULFs lists of ULF intensity trace over 
    time. Backgrounds and intensity traces at the same index correspond to the
    same ULF.
    
    Arguments: 
        ULFBackgrounds: List of background intensities [ba, bb, bc,....] where
                        ba is the background intensity of ulf a, etc. 
        intensities : List of lists of intensities [[ia1, ia2, ia3, ...],
                      [ib1, ib2, ib3...],...] where ia1 is the intensity of 
                      ulf a in frame 1, etc.
    """
    
    nInts = []
    for intList, ib in zip(intensities, ULFBackgrounds):
        nIntList = []
        for i in intList:
            newInt = i -ib
            if newInt < 0:
                newInt = 0
            nIntList.append(newInt)
        nInts.append(nIntList)      
    return nInts
        
def getBackgroundULFInensities(backgroundMat, tracks):
    """ Finds the intensity of each ULF in the background image. 

        Find the intensity of ULFs in the background image using the first 
        location measured in tracks.

        Arguments: 
             backgroundMat: Numpy array [x,y] of the background image. 
             tracks: list of track coordinates. 
       Returns:
             List of background intensities
    """
   
    # Make the backgroundMat image an [x,y,1] shaped array instead of a [x,y] 
    # shaped array.
    backgroundMat = [backgroundMat]
    backgroundMat = swapaxes(backgroundMat,0,1)
    backgroundMat = swapaxes(backgroundMat,1,2)
    
    # Get first location in each track
    firstPoints = []
    for t in tracks:
        x = t[0][0]
        y = t[1][0]
        firstPoints.append(([x],[y]))
    backgroundInts = getIntensities(backgroundMat, firstPoints)
    
    # Convert from a list of lists to a list of numbers. 
    outInts = []
    for b in backgroundInts:
        outInts.append(b[0])
    return outInts

def bleachCorrectInts(ints, bleachFrac):
    """ Uses the bleaching fraction calculated from a bleaching movie to correct
        intensity of each ULF. An alternative to ratio bleach correction. Not 
        used by Robert *et al.*
    """
    correctedInts = []
    for trace in ints:
        newTrace = []
        for i,b in zip(trace, bleachFrac):
            newTrace.append(float(i) / b);
        correctedInts.append(newTrace)
    return correctedInts  
        
def findSlopes(intensities, times = None, nPoints = 5):
    """ Finds the slopes of a list of lists  of intensities. 

       Calculates a linear regression using intensities  as the y-values and time
       as x-values (or indices if no time is provided) and returns the slope.
       Only the first nPoints are included in this calculation. Default is to 
       use the first 5 points. This corresponds to the first minute in a 
       time-lapse image series acquired at one frame every 15 seconds.
    
       Arguments: 
           intensities: list of list of intensities.
           time: list of time points.
           nPoints: only the first nPoints are used. 
       Returns: 
           List of slopes. 
     """ 
    if times == None:
        times = range(1, len(intensities[0]) + 1)
    slopes = []
    intercepts = []
    r_values = []
    for intTrace in intensities:
        slope, intercept, r_value, p_value, std_err = linregress(times[:nPoints], 
        intTrace[:nPoints])
        slopes.append(slope)
        intercepts.append(intercept)
        r_values.append(r_value)
    return(slopes, intercepts, r_values)

def findBleachFraction(mat):
    """ Calculates the fraction of the movie that has bleached in each frame.

        Not used in published results. 

        Arguments: 
            mat: numpy array [x,y,t]
        Returns: 
           List of fraction bleached in each frame. 
            
    """ 
    bleachMatInts = measureBleachingOverEntireMovie(mat)
    bleachFrac = []
    for i in bleachMatInts:
        p = float(i)/bleachMatInts[0]
        bleachFrac.append(p)
    return bleachFrac
    
def ratioBleachCorrect(movie):
    """ Performs ratio bleach correction on an image stack. 

        Corrects for bleaching by normalizing the intensity of each frame to that
        of the first frame.

       Arguments: 
           movie: Numpy array (x,y,t)
       Returns:
           Ratio bleach corrected movie. 
    """
    frame0 = movie[:,:,0]
    minVal = amax(nsmallest(100000,movie.flatten()))
    m0 = float(mean(frame0) - minVal)
    s = shape(movie)
    newMovie = zeros(s)
    newMovie[:,:,0] = movie[:,:,0] - minVal
    for t in range(1,s[2]):
        frame = movie[:,:,t] 
        frame = frame - minVal
        m = mean(frame)
        frame = m0 * (frame / m)
        newMovie[:,:,t] = frame  
    movie[movie <= 0] = 0  
    return newMovie

def removeIntensitiesInConvertedRegion(convertedRadius, intensities, distances):
    """Remove ULF intensities that are located in the converted region. 

       Removes values from the list intensities and removes corresponding values 
       from the list distance when distances are less than the converted radius. 
       This is used to eliminate ULFs in the photoconverted region from further
       calculations. 

       Arguments: 
           convertedRadius: Radius of the converted region. 
           intensities: List of list of intensities. 
           distances: List of distances. 
       Returns: 
           List of lists of intensities with intensities in the converted region
           removed and a list of corresponding distances. 
     """
    
    newDistances = []
    newIntensities = []
    for i in xrange(len(intensities)):
        if distances[i] > convertedRadius:
           newDistances.append(distances[i])
           newIntensities.append(intensities[i])
    return (newIntensities, newDistances)

def measureBleachingOverEntireMovie(mat):
    """Returns the average intensity of the each frame.
    """
    
    s = shape(mat)
    frameInts = []
    for t in range(s[2]):
        frameInts.append(mean(mat[:,:,t]))
    return frameInts

def normalizeIntensities(intensities, convertedIntensity):
    """Normalized intensities relative to the converted region. 
      
       Normalize intensities by dividing each intensity by the intensity in the
       converted region in the first frame. 
  
       Arguments: 
           intensities: List of intensities.
           convertedIntensity: Intensity in the converted region. 
      Returns: 
           List of normalized intensities. 
    """
    
    normedIntensities = deepcopy(intensities)
    for intList in normedIntensities:
        for i in range(len(intList)):    
            intList[i] = intList[i] / float(convertedIntensity)
    return normedIntensities
 

################
# Plot Results #
################
def plotTracks(tracks, redMat, greenMat, fileNumber, centerConversion, \
             radiusConversion, directory, number):
    """Plots the red and green channels with tracks and conversed region. 

       Creates and save a figure that displays the first frame of the red and 
       green channels with a circle identifying the converted region displayed 
       over the red channel and the ULF tracks displayed over the green channel.
    """ 
        
    fig1 = figure()
    
    # Display the red channel in the left subplot
    axLeft = fig1.add_subplot(1,2,1)
    axLeft.imshow(redMat[:,:,0], cmap=cm.get_cmap('gray'))
    axLeft.hold(True)
    axLeft.set_title('Red')
    
    # dispaly the green channel in the right subplot
    axRight = fig1.add_subplot(1,2,2)
    axRight.imshow(greenMat[:,:,0], cmap=cm.get_cmap('gray'))
    axRight.hold(True)
    axRight.set_title('Green')
    
    # Plot a circle over the red channel
    circ = Circle(centerConversion, radius = radiusConversion, color = 'r')
    circ.set_fill(None)
    axLeft.add_patch(circ)
    
    # Plot tracks over the green channel.
    for (x,y) in tracks:
        axRight.plot(x[:5],y[:5],'r')
       
    # Set axis limits to size of images. 
    s = shape(redMat) 
    axRight.set_xlim([0, s[1]])
    axRight.set_ylim([s[0], 0])
    axLeft.set_xlim([0, s[1]])
    axLeft.set_ylim([s[0], 0])
    
    # Save figure
    fig1.savefig(join(directory, ('tracks' + str(number) + '.png')))
    close(fig1)
    
def plotIntensities(directory, number, normedIntensities, ymax = 1.0):
    """ Plots intensity traces together. 

       Creates and saves a figure with a plot of all the intensities (y-axis is 
       labeled normalized intensities) on the same axis. ymax in the maximum
       value of the y-axis. 
    """
    
    fig2 = figure()
    ax = fig2.add_subplot(111)
    ax.set_ylim([0,ymax])
    ax.hold(True)
    for trace in normedIntensities:
        ax.plot(trace)
        
    ax.set_xlabel('frame')
    ax.set_ylabel('normalized intensity')
        
    #save figures
    fig2.savefig(join(directory, ('intensities' + str(number) + '.png')))
    close(fig2)

def plotByDistances(intensities, distances, directory, fileNumber = None):
    """Plots intensities grouped by distance. 

       Generates and saves a plot of intensities versus time (y-axis is 
       normalized intensities) grouped by distance. The distance ranges 50 to 
       100 pixels, 100 to 150 pixels, 150 to 200 pixels and 200 plus pixels are
       each displayed in a separate subplot. 
    """

    # Group intensities by distance.
    (range50minus, range50to100, range100to150, range150to200, range200plus)\
     = groupByDistance(intensities, distances)         
            
    # Create figure. 
    fig = figure()
    ax = fig.add_subplot(221)
    for trace in range50to100:
        ax.plot(trace)
        ax.set_xlabel('frame')
        ax.set_ylabel('normalized intensity')
        ax.set_title('50 to 100 pixels')
        ax.set_ylim([0,1.2])
    
    ax = fig.add_subplot(222)
    for trace in range100to150:
        ax.plot(trace)
        ax.set_xlabel('frame')
        ax.set_ylabel('normalized intensity')
        ax.set_title('100 to 150 pixels')
        ax.set_ylim([0,1.2])
    
    ax = fig.add_subplot(223)
    for trace in range150to200:
        ax.plot(trace)
        ax.set_xlabel('frame')
        ax.set_ylabel('normalized intensity')
        ax.set_title('150 to 200 pixels')
        ax.set_ylim([0,1.2])
    
    ax = fig.add_subplot(224)
    for trace in range200plus:
        ax.plot(trace)
        ax.set_xlabel('frame')
        ax.set_ylabel('normalized intensity')
        ax.set_title('200 plus pixels')
        ax.set_ylim([0,1.2])
      
    tight_layout()
    
    # Save figure
    if fileNumber == None:
        fileName = join(directory, 'all_intensities_by_distances.png')
    else:
        fileName = join(directory, ('intensities_by_distances' + \
        str(fileNumber) + '.png'))
    
    fig.savefig(fileName)
    close(fig)

def makeBleachingPlot(bleachingFile, nFrames = 5):
    """ Plots the average intensity in each frame. 

        Creates plots of the average intensity in each frame. Saves the plots 
        to the same directory as the bleaching movie. 
        Plots are annotated  with the slope of the line and the percent bleaching.
        Not used in the work published by Roberts *et al.* but used for trouble 
        shooting experiments. 
    """
    
    mat = readInTif(bleachingFile)
    b = measureBleachingOverEntireMovie(mat)
    # calculate percent bleached over nFrames
    p = 1 - b[nFrames] / float(b[0])
    p = p * 100
    # calculate slope over nFrames
    s = findSlopes([b],nPoints=nFrames)[0]
    
    directoryName =dirname(bleachingFile)
    filename = splitext(basename(bleachingFile))[0]
    
    fig = figure()
    ax = fig.add_subplot(111)
    ax.plot(b,'-o')
    ax.set_xlabel('Frame')
    ax.set_ylabel('Average Intensity')
    ax.set_title(filename)
    text1 = 'slope over first %d points: %.2f' % (nFrames, s[0])
    text2 = '%% bleach over first %d points: %.2f' % (nFrames, p)
    annotate(text1, xy=(0.3, 0.95), xycoords='axes fraction')
    annotate(text2, xy=(0.3, 0.9),  xycoords='axes fraction')
    fig.savefig(join(directoryName, filename + '.png'))
    close(fig)
    

    outStr =filename + ','
    for i in b:
        outStr = outStr + str(i) + ',' 
    print outStr
    
def plotConvertedIntensities(ints, redMovieName, nFrames = 5):
    """ Plots the intensities in the converted regions over time. 
    """

    # calculate percent bleached over nFrames
    p = 1 - ints[nFrames] / float(ints[0])
    p = p * 100
    # calcuate slope over nFrames
    s = findSlopes([ints],nPoints=nFrames)[0]
    
    directoryName =dirname(redMovieName)
    filename = splitext(basename(redMovieName))[0]
    
    fig = figure()
    ax = fig.add_subplot(111)
    ax.plot(ints,'r-o')
    ax.set_xlabel('Frame')
    ax.set_ylabel('Average Intensity in Converted Region')
    ax.set_title(filename)
    text1 = 'slope over first %d points: %.2f' % (nFrames, s[0])
    text2 = '%% change over first %d points: %.2f' % (nFrames, p)
    annotate(text1, xy=(0.3, 0.95), xycoords='axes fraction')
    annotate(text2, xy=(0.3, 0.9),  xycoords='axes fraction')
    fig.savefig(join(directoryName, filename + 'bleaching region.png'))
    close(fig)


################################
#  Write results to excel file #
################################
def writeDataToRow(worksheet, wbRow, wbCol, data):
    """Writes values to a row in an Excel file.

       Arguments: 
           worksheet: an Excel worksheet created by xlswriter
           wbRow: Index of the row.
           wbCol: Index of the column.
           data: List of numbers.
     """
    
    for n in data:
        worksheet.write(wbRow, wbCol, n)
        wbCol +=1
        
def writeVersionInfoToFile(worksheet,wbRow):
    """Writes git sha to an Excel file if available. 

       Arguments: 
           worksheet: Excel worksheet
           wbRow: Row index.
           wbCol: Column index.
    """
    
    try:
        homeDir = dirname(realpath(__file__))
        label = check_output(["git", "--git-dir=" + homeDir + '/.git', "log"])
        if len(label) > 150:
            label = label[:150]
    except:
        label = 'No version information available.'
        
    worksheet.write(wbRow,0,label)
    
def writeDataToFile(worksheet, fileNames, data):
    """Writes calculated parameters in data to the Excel worksheet. 

       Arguments:
           worksheet: Excel worksheet.
           fileName: Dictionary of files names. 
           data: Dictionary of data.
    """
    
    # Write version information 
    wbRow = 0
    writeVersionInfoToFile(worksheet,wbRow)
    wbRow += 1
   
    # Generate column labels that include frame numbers.    
    frameLabels = ['frame ' + str(x) for x in \
    range(len(data['convertedIntensities']))]

    # Write file names to the worksheet
    worksheet.write(wbRow, 0, 'Red Channel')
    worksheet.write(wbRow, 1, fileNames['redImageName'])
    wbRow += 1
    
    worksheet.write(wbRow, 0, 'Green Channel')
    worksheet.write(wbRow, 1, fileNames['greenImageName'])
    wbRow += 1

    worksheet.write(wbRow, 0, 'Coordinates File')
    worksheet.write(wbRow, 1, fileNames['CoordFileName'])
    wbRow += 1
    
    worksheet.write(wbRow, 0, 'Background File Name')
    worksheet.write(wbRow, 1, fileNames['backgroundName'])
    wbRow += 1
    
    worksheet.write(wbRow, 0, 'Bleaching File Name')
    worksheet.write(wbRow, 1, fileNames['bleachingFileName'])
    wbRow += 2
    
    # Write intensities in bleaching movie to file
    worksheet.write(wbRow, 0, 'Fraction unbleached in bleaching movie')
    writeDataToRow(worksheet, wbRow, 1, frameLabels)
    wbRow += 1
    writeDataToRow(worksheet, wbRow, 1, data['bleachedFraction'])
    wbRow += 2
    
    # Write convertered region center and radius to file.
    worksheet.write(wbRow, 0, 'Center of converted region: ')  
    worksheet.write(wbRow, 1, str(data['convertedCenter']))
    wbRow += 1
    worksheet.write(wbRow, 0, 'Radius of converted region: ')
    worksheet.write(wbRow, 1, data['convertedRadius'])
    wbRow += 1
    
    # Write converted region average background intensity to file.
    worksheet.write(wbRow, 0, 'Average intensity in converted region in \
                               background image')
    worksheet.write(wbRow,1, data['backgroundInConvertedRegion'])
    wbRow += 1
        
    # Write average intensities in converted region in each frame to file.
    worksheet.write(wbRow, 0, 'Average red intensity in the converted area.')
    writeDataToRow(worksheet, wbRow, 1, frameLabels)
    wbRow += 1
    writeDataToRow(worksheet, wbRow, 1, data['convertedIntensities'])
    wbRow += 2
    
    # Write average of entire frame red channel.
    worksheet.write(wbRow, 0, 'Average red intensity in entire frame.')
    wbRow += 1
    writeDataToRow(worksheet, wbRow, 1, frameLabels)
    wbRow += 1
    writeDataToRow(worksheet, wbRow, 1, data['redFrameInts'])
    wbRow += 2
    
    # Write average of entire frame green channel.
    worksheet.write(wbRow, 0, 'Average green intensity in entire frame.')
    wbRow += 1
    writeDataToRow(worksheet, wbRow, 1, frameLabels)
    wbRow += 1
    writeDataToRow(worksheet, wbRow, 1, data['greenFrameInts'])
    wbRow += 2
    
    # Write ULF background intensities to file.
    worksheet.write(wbRow, 0, 'Background Intensiteis of ULFs')
    wbRow += 1
    worksheet.write(wbRow, 0, 'Distance')
    worksheet.write(wbRow, 1, 'Intensity')
    wbRow += 1
    distances = data['distances']
    backgroundULFs = data['backgroundULFs']
    for n in range(len(backgroundULFs)):
        worksheet.write(wbRow, 0, distances[n])

        worksheet.write(wbRow, 1, backgroundULFs[n])
        wbRow += 1
        
    wbRow += 2
   
    # Write slopes and distances with the corresponding intensities (all in the
    # same row) to file. 
    
    # Write column titles.
    worksheet.write(wbRow, 0, 'Average red intensity in the ULF regions:')
    wbRow += 1
    worksheet.write(wbRow, 0, 'Distance (pixels)')
    worksheet.write(wbRow, 1, 'Slope  (intensity / frame)')
    writeDataToRow(worksheet, wbRow, 2, frameLabels)
    wbRow +=1
    
    # Write data.
    distances = data['distances']
    slopes = data['slopes']
    intensities = data['intensities']
    for n in range(len(data['intensities'])):
        worksheet.write(wbRow, 0, distances[n])
        worksheet.write(wbRow, 1, slopes[n])
        writeDataToRow(worksheet, wbRow, 2, intensities[n])
        wbRow += 1
    wbRow += 2
    
    # Write normalized slopes and distances with the corresponding normalized 
    # intensities (all in the same row) to file.
    
    # Write column titles.
    worksheet.write(wbRow, 0, 'Normalized Average red intensity in the ULF \
                    regions:')
    wbRow += 1
    worksheet.write(wbRow, 0, 'Distance (pixels)')
    worksheet.write(wbRow, 1, 'Slope (intensity / frame)') 
    writeDataToRow(worksheet, wbRow, 2, frameLabels)
    wbRow +=1
   
    # Write data.
    normedIntensities = data['normedIntensities']
    normSlopes = data['normSlopes']
    for n in range(len(normedIntensities)):
        worksheet.write(wbRow, 0, distances[n])
        worksheet.write(wbRow, 1, normSlopes[n])
        writeDataToRow(worksheet, wbRow, 2, normedIntensities[n])

        wbRow += 1

def writeAllDataToFile(workbook, allDistances, allSlopes):
    """Records all data from batch processing to an excel file. 
       
       Arguments: 
          workbook: Excel workbook.
          allDistances: list of distances.
          allSlopes: List of corresponding slopes. 
    """
     
    # Group slopes and distances by distance.  
    (range50minus, range50to100, range100to150, range150to200, range200plus) \
    = groupByDistance(allSlopes, allDistances)                             
    
    # Write to excel file. 
    worksheet0 = workbook.add_worksheet()
    worksheet0.write(0,0,'slopes in the range 0 to 50 pixels')
    rn = 1
    for s in range50minus:
        worksheet0.write(rn,0,s)
        rn += 1
    
    worksheet1 = workbook.add_worksheet() 
    worksheet1.write(0,0,'slopes in the range 50 to 100 pixels')
    rn = 1
    for s in range50to100:
        worksheet1.write(rn,0,s)
        rn += 1
        
    worksheet2 = workbook.add_worksheet() 
    worksheet2.write(0,0,'slopes in the range 100 to 150 pixels')
    rn = 1
    for s in range100to150:
        worksheet2.write(rn,0,s)
        rn += 1
    
    worksheet3 = workbook.add_worksheet() 
    worksheet3.write(0,0,'slopes in the range 150 to 200 pixels')
    rn = 1
    for s in range150to200:
        worksheet3.write(rn,0,s)
        rn += 1
        
    worksheet4 = workbook.add_worksheet() 
    worksheet4.write(0,0,'slopes in the range 200 plus pixels')
    rn = 1
    for s in range200plus:
        worksheet4.write(rn,0,s)
        rn += 1
        
    workbook.close() 
    

##############
# Exceptions #
##############
class _fileFormatError(Exception):
    """Exception raised when a file does not have the correct format.
    """
    def __init__(self, fileName):
        self.fileName = fileName
    def __str__(self):
        return repr(self.fileName)
    
class  _missingFileError(Exception):
    """Exception raised when a file does not exist. Returns the file name.
    """
     
    def __init__(self, fileName):
        self.fileName = fileName
    def __str__(self):
        return repr(self.fileName)
        
class _imageShapeDiscrepencyError(Exception):
    """Exception raised when two images are not the same shape.
    """
    def __init__(self, fileOne, fileTwo):
        self.files = fileOne + ' and ' + fileTwo
    def __str__(self):
        return repr(self.files)


##################
# Batch Analysis #
##################
def openFilesInSet(redImageName, greenImageName):
    """ Opens the red and green files in a data set.
     
        Arguments: 
            redImageName: Name of the red image.
            greenImageName: Name of the green image.
        Returns: 
            Tuple of numpy arrays (red, green).
    """
    data = {}
    try:
        redMat =  readInTif(redImageName)
    except IOError:
        raise _missingFileError(redImageName) 
        return data
    
    try:       
        greenMat = readInTif(greenImageName)
    except IOError:
        raise _missingFileError(greenImageName) 
        return data
    return (redMat, greenMat)

def analyzeDataSet(CoordFileName, redImageName, greenImageName, resultsName, 
                convertedCenter, convertedRadius, fileNumber, dataDirectory,
                backgroundName = None, bleachingFileName = None):
    """Analyze one data set.  

       Analyzes one dataset including a red image series, a green image series, 
       a coordinates file, and a background image file. This is where the 5 step 
       workflow for ULF intensity analysis is implemented.
    """
        
    # Open files. Return if files do not exist.
    (redMat, greenMat) = openFilesInSet(redImageName, greenImageName)
        
    data = {}
    # Check if images are the same size
    sRed = shape(redMat)
    sGreen = shape(greenMat)
    if sRed != sGreen:
        raise _imageShapeDiscrepencyError(redImageName, greenImageName)
        return data
        
    # WORKFLOW STEP 1: Bleach correction
    # ratio bleach correction
    redMat = ratioBleachCorrect(redMat)
    nFrames = sRed[2]  
    
    # WORKFLOW STEP2: ULF intensity calculation
    
    try: 
        tracks = readInCoordinatesFile(CoordFileName, nFrames)
    except IOError:
        raise _missingFileError(CoordFileName) 
        return data     
    
    # Get the red channel intesity of each ULF in each frame. ULF locations are 
    # defined  by tracks (created by readInCoordinatesFile) 
    intensities = getIntensities(redMat, tracks)  
    # print(intensities[0][0:10])      
    
    # Get the distnaces between the center of the converted region and the
    # first point in each track. 
    distances = getDistancesToFirstPointInTrack(tracks, convertedCenter)
    
    # Get the average intensity in the converted region in each frame. 
    convertedIntensities = getAverageIntensityInCircleOverTime(redMat, 
    convertedCenter, convertedRadius)
    
    plotConvertedIntensities(convertedIntensities,redImageName)
    
    # WORKFLOW STEP 3: Background subtraction
    #Background subtraction 
    if backgroundName!= None:
        try: 
            backgroundMat = readInTif(backgroundName)  
        except IOError:
            raise _missingFileError(backgroundName) 
            return data  
        # Check if background image is the same shape as the red and green
        # images.
        
        # Find the intensity of each ULF in the background frame. 
        backgroundULFInts = getBackgroundULFInensities(backgroundMat, tracks) 
        # Subtract the intensity of the ULF in the background frame from the
        # Intensity in subsequent frames.     
        nInts = subtractULFBackground(backgroundULFInts, intensities)
        backgroundInConvertedRegion = getAverageIntensityInCircle(backgroundMat, 
        convertedCenter, convertedRadius)
    else:  
        backgroundInConvertedRegion = 0
        nInts = intensities
        backgroundULFInts = []
    
    # WORKFLOW STEP 4: Normalization
    # Normalize relative to the intensity of the converted region in the first 
    # frame minus the intensity of the converted region in the background frame. 
    normedIntensities = normalizeIntensities(nInts, 
    convertedIntensities[0] - backgroundInConvertedRegion) 
    # normedIntensities = normalizeIntensities(nInts, convertedIntensities[0]) 

    # WORKFLOW STEP 5: Visulaize and report results   
    # Get slopes of raw and normalized intensities. 
    (slopes, intercepts, r_values) = findSlopes(intensities)
    (normSlopes, normIntercepts, normR_values) = findSlopes(normedIntensities)
    
    # Get bleaching rate (intensity in each entire frame in each frame)
    redFrameInts = measureBleachingOverEntireMovie(redMat)
    greenFrameInts = measureBleachingOverEntireMovie(greenMat)
    
    # Write caluluated parameters to excel file. 
    fileNames = {}
    fileNames['redImageName'] = redImageName
    fileNames['greenImageName'] =greenImageName
    fileNames['CoordFileName'] = CoordFileName
    fileNames['backgroundName'] = backgroundName
    
    if bleachingFileName == None:
        fileNames['bleachingFileName'] = 'None'
    else:
        fileNames['bleachingFileName'] =bleachingFileName
    data = {}
    data['convertedIntensities'] = convertedIntensities
    data['intensities'] = intensities
    data['distances'] = distances
    data['normedIntensities'] = normedIntensities
    data['slopes'] = slopes
    data['normSlopes'] = normSlopes
    data['redFrameInts'] = redFrameInts
    data['greenFrameInts'] = greenFrameInts
    data['backgroundULFs'] = backgroundULFInts
    data['convertedCenter'] = convertedCenter
    data['convertedRadius'] = convertedRadius 
    data['backgroundInConvertedRegion'] = backgroundInConvertedRegion
    data['bleachedFraction'] = []
    
        
    workbook = xlsxwriter.Workbook(resultsName)
    worksheet = workbook.add_worksheet()
    writeDataToFile(worksheet, fileNames, data)
    
    
    # Plot raw data and results. 
    plotTracks(tracks, redMat, greenMat, fileNumber, convertedCenter, 
               convertedRadius, dataDirectory, fileNumber)
    plotIntensities(dataDirectory, fileNumber, normedIntensities)
    plotByDistances(normedIntensities, distances, dataDirectory, fileNumber)  
    workbook.close()
    
    return (data)
    
def analyzeRedAndGreenSet(dataSet):
    """ Red and Green channel comparison for multiple data sets. 
    """
    
    CoordFileName = dataSet['cordName']
    redImageName = dataSet['redName']
    greenImageName = dataSet['greenName']
    resultsName = dataSet['resultsName']
    convertedCenter = dataSet['covertedCenter']
    redBackgroundName = dataSet['redBackgroundName']
    greenBackgroundName = dataSet['greenBackgroundName']
    # Open files. Return if files do not exist.
    (redMat, greenMat) = openFilesInSet(redImageName, greenImageName)
        
    data = {}
    # Check if images are the same size
    sRed = shape(redMat)
    sGreen = shape(greenMat)
    if sRed != sGreen:
        raise _imageShapeDiscrepencyError(redImageName, greenImageName)
        return data
        #
    # Skip bleach correction.
    
    nFrames = sRed[2]  
    
    # WORKFLOW STEP2: ULF intensity calculation
    try: 
        tracks = readInCoordinatesFile(CoordFileName, nFrames)
    except IOError:
        raise _missingFileError(CoordFileName) 
        return data     
    
    # Get the red and green channel intesity of each ULF in each frame. ULF locations are 
    # defined  by tracks (created by readInCoordinatesFile) 
    RedIntensities = getIntensities(redMat, tracks) 
    GreenIntensities = getIntensities(greenMat, tracks) 
    
    # Get the distnaces between the center of the converted region and the
    # first point in each track. 
    distances = getDistancesToFirstPointInTrack(tracks, convertedCenter)
    
    # WORKFLOW STEP 3: Background subtraction
    #Background subtaction  
    try: 
        backgroundMatRed = readInTif(redBackgroundName)  
        backgroundMatGreen = readInTif(greenBackgroundName)
    except IOError:
        raise _missingFileError(redBackgroundName) 
        return data  
    
    # Find the intensity of each ULF in the background frame. 
    backgroundULFIntsRed = getBackgroundULFInensities(backgroundMatRed, tracks)
    backgroundULFIntsGreen = getBackgroundULFInensities(backgroundMatGreen, tracks)
      
    # Subtract the intensity of the ULF in the background frame from the
    # Intensity in subsequent frames.     
    nIntsRed = subtractULFBackground(backgroundULFIntsRed, RedIntensities)
    nIntsGreen = subtractULFBackground(backgroundULFIntsGreen, GreenIntensities)
    
    # No normalization for now! 2/18/15

    # WORKFLOW STEP 5: Visualize and report results   
    # Get slopes of raw and normalized intensities. 
    (slopesRed, interceptsRed, r_valuesRed) = findSlopes(nIntsRed)
    (slopesGreen, interceptsGree, r_valuesGreen) = findSlopes(nIntsGreen)

    workbook = xlsxwriter.Workbook(resultsName)
    redSheet = workbook.add_worksheet(name='Red')
    greenSheet = workbook.add_worksheet(name='Green')
    row = 0
    col = 0
    redSheet.write(row, col, 'Distances')
    redSheet.write(row, col + 1, 'Slope')
    redSheet.write(row, col + 2, 'Intensities')
    greenSheet.write(row, col, 'Distances')
    greenSheet.write(row, col + 1, 'Slope')
    greenSheet.write(row, col + 2, 'Intentsites')
    row+=1
    for ind, d in enumerate(distances):
        # Write out distances
        redSheet.write(row, col, d)
        greenSheet.write(row, col, d)
        # Write out slopes
        redSheet.write(row, col + 1, slopesRed[ind])
        greenSheet.write(row, col + 1, slopesGreen[ind])
        # Write out intensities
        writeDataToRow(redSheet, row, col + 2, nIntsRed[ind])
        writeDataToRow(greenSheet, row, col + 2, nIntsGreen[ind])
        row+=1
    workbook.close() 
    
    return ({})
         
def groupByDistance(vals, distances):
    """ Group ULF intensity traces by distance from the converted region.

        Arguments:
            vals: List of list of intensities.
            distances: List of distances.
       Returns: 
           Tuple of 5 lists corresponding to intensity traces in the 
           following ranges: 0 to 50, 50 to 100, 100 to 150, 150 to 200,
           greater then 200.  
    """
        
    range50minus = []
    range50to100 = []
    range100to150 = []
    range150to200 = []
    range200plus = []
    for intVal, d in zip(vals, distances):
        if (d < 50):
            range50minus.append(intVal)
        elif (d > 50) and (d < 100):
            range50to100.append(intVal)
        elif (d >= 100) and (d < 150):
            range100to150.append(intVal)
        elif (d >= 150) and (d < 200):
            range150to200.append(intVal)
        elif (d >= 200):
            range200plus.append(intVal)
    return (range50minus, range50to100, range100to150, range150to200, 
    range200plus)
    
def analyzeFilesTogether(arguments, backgroundNames = None):
    """Analyzes multiple data sets together. Calls analyzeDataSet for each 
       data set and also groups results from all data sets. 
    """
    
    coordFileNames = arguments['coordFileNames']
    redImageNames = arguments['redImageNames']
    greenImageNames = arguments['greenImageNames']   
    covertedCenter = arguments['covertedCenter']
    convertedRadius = arguments['convertedRadius'] 
    dataDirectory = arguments['dataDirectory']
    FileNumbers = arguments['FileNumbers']
    ResultsNames = arguments['ResultsNames']
    bleachingFile = arguments['bleachingFileName']
    
    if backgroundNames == None:
        backgroundNames = []
        for x in xrange(len(coordFileNames)):
            backgroundNames.append(None)
           
    allDistances = []
    allSlopes = []
    allNormedIntensities = [] 
    unProcessedFiles = []
    shapeDiscrepencenyFiles = []
    normedIntensities = []
    normSlopes = []
    allDistances = []
     
    for cordName, redName, greenName, bgName, FileNumber, resultsName in \
        zip(coordFileNames, redImageNames, greenImageNames, backgroundNames,\
        FileNumbers, ResultsNames): 
        try:
            data = analyzeDataSet(cordName, redName, greenName, resultsName,\
            covertedCenter, convertedRadius, FileNumber, dataDirectory, bgName,\
            bleachingFile)
        except _missingFileError as fName:
            unProcessedFiles.append(fName)   
        except _imageShapeDiscrepencyError as imageNamges: 
            shapeDiscrepencenyFiles.append(imageNamges)       
        else:
            distances = data['distances']
            normedIntensities = data['normedIntensities']
            normSlopes = data['normSlopes']
            allDistances = allDistances + distances
            allSlopes = allSlopes + normSlopes
            allNormedIntensities = allNormedIntensities + normedIntensities
    
    allResultsName = dataDirectory + '/allSlopes.xlsx'
    
    # Print warning if there were missing files.
    if len(unProcessedFiles) > 0:  
        warnStr = 'The following files do not exists. These data sets were \
        ignored.'
        for f in unProcessedFiles:
            warnStr = warnStr + '\n' + f.fileName
        print(warnStr)
    
    # Print warning if there were data sets with images of different sizes
    if len(shapeDiscrepencenyFiles) > 0:
        warnStrShape = 'The following files contained images with different \
        shapes. These data sets were ignored.'
        for f in shapeDiscrepencenyFiles:
            warnStrShape = warnStrShape + '\n' + f.files
        print(warnStrShape)
        
    if len(allNormedIntensities) > 0: # If at least some of the files where 
                                      # succesfully processed. 
        workbook = xlsxwriter.Workbook(allResultsName)
        writeAllDataToFile(workbook, allDistances, allSlopes)
        plotByDistances(allNormedIntensities, allDistances, dataDirectory, 
        fileNumber = None)
        workbook.close()
        
def analyzeRedAndGreen(arguments):
    """ Red and Green channel comparison for multiple data sets. 
    """
    
    coordFileNames = arguments['coordFileNames']
    redImageNames = arguments['redImageNames']
    greenImageNames = arguments['greenImageNames']   
    covertedCenter = arguments['covertedCenter']
    #convertedRadius = arguments['convertedRadius'] 
    #dataDirectory = arguments['dataDirectory']
    #FileNumbers = arguments['FileNumbers']
    ResultsNames = arguments['ResultsNames']
    redBackgroundsRed = arguments['allBackgroundsRed'] 
    allBackgroundsGreen = arguments['allBackgroundsGreen']
    
    # Initalize empty background names list.
    backgroundNames = []
    for x in xrange(len(coordFileNames)):
        backgroundNames.append(None)
           
    #allDistances = []
    #allSlopes = []
    #allNormedIntensities = [] 
    unProcessedFiles = []
    shapeDiscrepencenyFiles = []
    #normedIntensities = []
    #normSlopes = []
    #allDistances = []
    
    # Names of arguments to pass to analyzeRedAndGreenSet
    keys = ['cordName', 'redName', 'greenName', 'resultsName',\
            'redBackgroundName', 'greenBackgroundName']
    # group by by data set and analyze.
    for  setInputs in zip(coordFileNames, redImageNames, greenImageNames, \
                      ResultsNames, redBackgroundsRed, allBackgroundsGreen): 
        try:
            argIn = dict(zip(keys, setInputs))
            argIn['covertedCenter'] = covertedCenter
            data = analyzeRedAndGreenSet(argIn)
        except _missingFileError as fName:
            unProcessedFiles.append(fName)   
        except analyzeRedAndGreenSet as imageNamges: 
            shapeDiscrepencenyFiles.append(imageNamges)       
        else:
            pass
            #distances = data['distances']
            #normedIntensities = data['normedIntensities']
            #normSlopes = data['normSlopes']
            #allDistances = allDistances + distances
            #allSlopes = allSlopes + normSlopes
            #allNormedIntensities = allNormedIntensities + normedIntensities
    
    #allResultsName = dataDirectory + '/allSlopes.xlsx'
    
    # Print warning if there were missing files.
    if len(unProcessedFiles) > 0:  
        warnStr = 'The following files do not exists. These data sets were \
        ignored.'
        for f in unProcessedFiles:
            warnStr = warnStr + '\n' + f.fileName
        print(warnStr)
    
    # Print warning if there were data sets with images of different sizes
    if len(shapeDiscrepencenyFiles) > 0:
        warnStrShape = 'The following files contained images with different \
        shapes. These data sets were ignored.'
        for f in shapeDiscrepencenyFiles:
            warnStrShape = warnStrShape + '\n' + f.files
        print(warnStrShape)
    """    
    if len(allNormedIntensities) > 0: # If at least some of the files where 
                                      # succesfully processed. 
        workbook = xlsxwriter.Workbook(allResultsName)
        writeAllDataToFile(workbook, allDistances, allSlopes)
        plotByDistances(allNormedIntensities, allDistances, dataDirectory, 
        fileNumber = None)
        workbook.close()
    """
        
def runBatch(fileNumbers, directory, rootName, covertedCenter,\
            convertedRadius = 50, background = True,\
            bleachingFile = None):        
    """ **This is the primary function to use for analyzing an experiment.*** 
        runBatch performs  a complete analysis on multiple sets of data. Each 
        data file needs to be named with the following convention 
        for run batch to work: 
            
        1. rootname + red + number + .tif
        2. rootname + green + number + .tif
        3. rootname + before + number + .tif
        4. rootname + .txt 
        
        
        Where rootname is a name that is the same for all files in the batch and 
        number is unique. See the explanation at the top of ulfexchange.py and 
        Robert *et al.* for a complete explanation of ULF intensity analysis.
         
    """ 
    
    allCoords = []
    allRed = []
    allGreens = []
    allResults = []
    allBackgrounds = []
    
    for fileNumber in fileNumbers:
        CoordFileName = join(directory, (rootName + fileNumber + '.txt'))
        redImageName = join(directory, (rootName + fileNumber + 'red.tif'))
        greenImageName =join(directory, (rootName + fileNumber + 'green.tif'))
        resultsName = join(directory, (rootName + '_results' + fileNumber + \
        '.xlsx'))
        backgroundName = join(directory, (rootName + fileNumber + 'before.tif'))
        
        allCoords.append(CoordFileName)
        allRed.append(redImageName)
        allGreens.append(greenImageName)
        allResults.append(resultsName)
        allBackgrounds.append(backgroundName)
       
    arguments = {}     
    arguments['coordFileNames'] = allCoords
    arguments['redImageNames'] = allRed
    arguments['greenImageNames'] = allGreens
    arguments['covertedCenter'] = covertedCenter
    arguments['convertedRadius'] = convertedRadius
    arguments['dataDirectory'] = directory
    arguments['FileNumbers'] = fileNumbers
    arguments['ResultsNames'] = allResults
    arguments['bleachingFileName'] = bleachingFile
    
    if background:
        backgroundNames = allBackgrounds
    else:
        backgroundNames = None
        
    analyzeFilesTogether(arguments, backgroundNames)
    
def runBatchRedAndGreen(fileNumbers, directory, rootName, covertedCenter,\
            convertedRadius):        
    """ Red and Green channel comparison. 

        Creates a plot and Excel file to show bleaching in red and green channels.
            
        1. rootname + red + number + .tif
        2. rootname + green + number + .tif
        3. rootname + before + number + .tif
        4. rootname + .txt    
    """ 
    
    allCoords = []
    allRed = []
    allGreens = []
    allResults = []
    allBackgroundsRed = []
    allBackgroundsGreen = []
    
    for fileNumber in fileNumbers:
        CoordFileName = join(directory, (rootName + fileNumber + '.txt'))
        redImageName = join(directory, (rootName + fileNumber + 'red.tif'))
        greenImageName =join(directory, (rootName + fileNumber + 'green.tif'))
        resultsName = join(directory, (rootName + '_redgreenresults' + fileNumber + \
        '.xlsx'))
        redBackgroundName = join(directory, (rootName + fileNumber + 'before.tif'))
        greenBackgroundName = join(directory, 'background green.tif') # Always the same name
        
        allCoords.append(CoordFileName)
        allRed.append(redImageName)
        allGreens.append(greenImageName)
        allResults.append(resultsName)
        allBackgroundsRed.append(redBackgroundName)
        allBackgroundsGreen.append(greenBackgroundName)
       
    arguments = {}     
    arguments['coordFileNames'] = allCoords
    arguments['redImageNames'] = allRed
    arguments['greenImageNames'] = allGreens
    arguments['covertedCenter'] = covertedCenter
    arguments['convertedRadius'] = convertedRadius
    arguments['dataDirectory'] = directory
    arguments['FileNumbers'] = fileNumbers
    arguments['ResultsNames'] = allResults
    arguments['allBackgroundsRed'] = allBackgroundsRed
    arguments['allBackgroundsGreen'] = allBackgroundsGreen
        
    analyzeRedAndGreen(arguments)
    
