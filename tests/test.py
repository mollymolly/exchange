from numpy import ones, mean
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from ulfexchange import  getAverageIntensityInCircle,  \
                        getDistancesToFirstPointInTrack,\
                        removeIntensitiesInConvertedRegion,\
                        normalizeIntensities, \
                        findSlopes, readInCoordinatesFile, getIntensities, \
                        distance, getAverageIntensityInCircleOverTime, \
                        measureBleachingOverEntireMovie, \
                        getBackgroundULFInensities, subtractULFBackground, \
                        findBleachFraction, bleachCorrectInts, ratioBleachCorrect
                        
homeDir = os.path.dirname(os.path.realpath(__file__))

# Test ReadInTi
# TODO

# Test readInCoordinatesFile
tracks = readInCoordinatesFile(homeDir + '/testCoords.txt', 6)
if tracks[0][0][0] == 289.003846:
    print 'pass readInCoordinatesFile'
else:
    print 'fail readInCoordinatesFile'
   
# Test distance
d = distance((0,0),(3,4))
if d == 5:
    print 'pass test distance'
else:
    print 'fail test distance'
d = distance((5,5),(8,9))
if d == 5:
    print 'pass test distance'
else:
    print 'fail test distance'

# Test getAverageIntensityInCircle
frame = ones((10,10))
frame[8,5]= 100
intensity1 = getAverageIntensityInCircle(frame, (5,8), 1)
intensity2 = getAverageIntensityInCircle(frame, (7,7), 2) 
if intensity1 == mean([1.0, 1.0, 100.0, 1.0, 1.0]):
    print 'pass getAverageIntensityInCircle'
else:
    print 'fail getAverageIntensityInCircle'
if intensity2 == mean([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]):
    print 'pass getAverageIntensityInCircle'
else:
    print 'fail getAverageIntensityInCircle'

# Test getAverageIntensityInCircleOverTime
imMat = ones((10,10,3))
imMat[5,5,:] = 100
imMat[7,6,1] = 2
imMat[7,6,2] = 3 
intensities1 = getAverageIntensityInCircleOverTime(imMat, (5,5), 1)
intensities2 = getAverageIntensityInCircleOverTime(imMat, (6,7), 1)
intensities3 = getAverageIntensityInCircleOverTime(imMat[:,:,0], (5,5), 1)
i1 = mean([1.0, 1.0, 100.0, 1.0, 1.0])
i2a = 1.0
i2b = mean([1.0, 1.0, 2, 1.0, 1.0])
i2c = mean([1.0, 1.0, 3, 1.0 ,1.0])

if intensities1 == [i1, i1, i1]:
    print 'pass getAverageIntensityInCircleOverTime'
else: 
    print 'fail getAverageIntensityInCircleOverTime'
    
if intensities2 == [i2a, i2b, i2c]:
    print 'pass getAverageIntensityInCircleOverTime'
else: 
    print 'fail getAverageIntensityInCircleOverTime'
    
if intensities3 == [i1]:
    print 'pass getAverageIntensityInCircleOverTime'
else: 
    print 'fail getAverageIntensityInCircleOverTime'
    
    
# test processData

data = ones((100,100,6))
data[:,:,1] = 2
data[:,:,2] = 3
data[:,:,3] = 4
data[:,:,4] = 5
data[0:50,:,5] = 6
data[51:-1,:,5] = 100

# set up some tracks for testing
x1 = [10, 10, 10, 10, 10, 10]
y1 = [20, 20, 20, 20, 20, 20]
x2 = [70, 70, 70, 70, 70, 70]
y2 = [10, 10, 10, 10, 10, 10]
x3 = [50, 100, 100, 100, 100, 100]
y3 = [50, 100, 100, 100, 100, 100]
tracks = [(x1,y1),(x2,y2),(x3,y3)]

# Test getDistancesToFirstPointInTrack
dists = getDistancesToFirstPointInTrack(tracks, (50,50))

if dists == [50.0, 44.721359549995796, 0.0]:
    print 'pass getDistancesToFirstPointInTrack'
else: 
    print 'fail getDistancesToFirstPointInTrack'
    
# test normalizeIntensities
intensities = [[1,2,3,4,5],[100,200,300,400,500]]

normedIntensities = normalizeIntensities(intensities, 10)
if [[0.1, 0.2, 0.3, 0.4, 0.5], [10.0, 20.0, 30.0, 40.0, 50.0]] == normedIntensities:
    print 'pass normalizeIntensities'
else:
    print 'fail normalizeIntensities'


# test removeIntensitiesInConvertedRegion(
intensities = [[1,2,3], [4,5,6], [7,8,9], [1,4,5]]
distances = [1, 6, 7, 2]
(newIntensities, newDistances) = removeIntensitiesInConvertedRegion(5, intensities, distances)
if (newIntensities == [[4, 5, 6], [7, 8, 9]]) and (newDistances == [6, 7]):
    print 'pass removeIntensitiesInConvertedRegion'
else:
    print 'fail removeIntensitiesInConvertedRegion'


# test find slopes Only tests slopes
intensities = [[1,2,3,4,5],[100,200,300,400,500]]
(dataSlopes, dataIntercepts, r_values) = findSlopes(intensities)
(dataSlopes2, dataIntercepts2, r_value2s) = findSlopes(intensities, times = [0,2,4,6,8])
if (dataSlopes == [1.0, 100.0]):
    print 'pass findSlopes'
else:
    print 'fail findSlopes'
    
if (dataSlopes2 == [0.5, 50.0]):
        print 'pass findSlopes'
else:
    print 'fail findSlopes'
    
# test getIntensities(mat, tracks, npoints, radius = 6)
data2 = ones((100,100,3))
data2[10:20, 50:60, 0] = 3
data2[20:30, 60:70, 1] = 4
data2[30:40, 70:80, 2] = 5
x1 = [55, 65, 75]
y1 = [15, 25, 35]
tracks1 = [(x1,y1)]
x2 = [50, 51, 52]
y2 = [30, 31,32]
tracks2 = [(x2,y2)]

#ax = subplot(111)
#imshow(data2[:,:,1])
#ax.set_ylim([0,100])
#ax.set_xlim([0,100])

#hold(True)
#plot(x1,y1,'go')
#plot(x2,y2,'y^')

intensities1 = getIntensities(data2, tracks1, radius = 3)
if intensities1 == [[3.0, 4.0, 5.0]]:
    print 'pass getIntensities'
else:
    print 'fail getIntensities'
intensities2 = getIntensities(data2, tracks2, radius = 3)

if intensities2 == [[1.0, 1.0, 1.0]]:
    print 'pass getIntensities'
else:
    print 'fail getIntensities'
    
data3 = ones((100,100,3))
data3[:,:,1] = 2
data3[:,:,2] = 3
frameInts = measureBleachingOverEntireMovie(data3)
if frameInts == [1.0, 2.0, 3.0]:
    print 'pass measureBleachingOverEntireMovie'
else:
    print 'fail measureBleachingOverEntireMovie'

# set up some tracks for testing
x1 = [33, 10, 10, 10, 10, 10]
y1 = [23, 20, 20, 20, 20, 20]
x2 = [75, 70, 70, 70, 70, 70]
y2 = [15, 10, 10, 10, 10, 10]
x3 = [50, 100, 100, 100, 100, 100]
y3 = [50, 100, 100, 100, 100, 100]
tracks = [(x1,y1),(x2,y2),(x3,y3)]

bk = ones((100,90))
bk[23,33] = 100
bk[15,75] =22
bkInts = getBackgroundULFInensities(bk, tracks)
if bkInts == [1.8761061946902655, 1.1858407079646018, 1.0]:
    print 'pass getBackgroundULFIntensities'
else:
    print 'fail getBackgroundULFIntensities'
    
intensities = [[1,2,3], [10,10, 10], [4,4,4], [9, 8, 7]]
ULFBackgrounds = [1,3,2,10]    
correctedInts = subtractULFBackground(ULFBackgrounds, intensities)
if correctedInts == [[0, 1, 2], [7, 7, 7], [2, 2, 2], [0, 0, 0]]:
    print 'pass subtractULFBackground'
else:
    print 'fail subtractULFBackground'
    
movie = ones((10,10,5))
movie[:,:,0] = 10
movie[:,:,1] = 8
movie[:,:,2] = 7
movie[:,:,3] = 5
movie[:,:,4] = 4

i = findBleachFraction(movie)
if i ==  [1.0, 0.80000000000000004, 0.69999999999999996, 0.5, 0.40000000000000002]:
    print 'pass findPercentBleached'
else:
    print 'fail findPercentBleached'
    
bints = bleachCorrectInts([[100, 90, 80, 70, 50]], [1, 0.8, 0.8, 0.5, 0.2])

redMat = ones((10,10,5))
redMat[:,:,0] = 5
redMat[:,:,1] = 4
redMat[:,:,2] = 3
ratioBleachCorrect(redMat)

solution = 5*ones((10,10,5))
'''
if all(redMat) == all(solution):
    print 'pass ratioBleachCorrect'
else:
    print 'fail ratioBleachCorrect' '''  
redMat = ones((10,10,5))
redMat[:,:,0] = 5
redMat[5,5,0] = 6
redMat[:,:,1] = 4
redMat[5,5,1] = 10
redMat[:,:,2] = 3
ratioBleachCorrect(redMat)

