import unittest
import os
import sys
# Import ulfexchange from one directory up
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import ulfexchange as ulfx
import numpy as np

class TestSequenceFunctions(unittest.TestCase):
    
    def setUp(self):
        self.homeDir = os.path.dirname(os.path.realpath(__file__))
    
    ############
    # File I/O #
    ############  
    def test_readInCoordinatesFiles(self):
        # Test reading in coordinates files
        nFrames = 6; # only read in the first 6 lines
        tracks = ulfx.readInCoordinatesFile(self.homeDir + '/testCoords.txt', \
                nFrames)
        self.assertEqual(tracks[0][0][0], 289.003846)
        
    def test_readInTif(self):
        fname = 'test.tif'
        data = ulfx.readInTif(fname)
        
    ######################
    # Basic Calculations #
    ######################
    def test_distance(self):
        # Test calculating distance between two points
        d1 = ulfx.distance((0,0),(3,4))
        self.assertEqual(d1,5)
        d2 = ulfx.distance((5,5),(8,9))
        self.assertEqual(d2,5)
        d3 = ulfx.distance((0,0),(0,0))
        self.assertEqual(d3,0)
        
    def test_getAverageIntensityInCircle(self):
        
        frame = np.ones((10,10)) # frame of an image
        frame[8,5]= 100
        intensity1 = ulfx.getAverageIntensityInCircle(frame, (5,8), 1)
        intensity2 = ulfx.getAverageIntensityInCircle(frame, (7,7), 2) 
        self.assertEqual(intensity1, np.mean([1.0, 1.0, 100.0, 1.0, 1.0]))
        self.assertEqual(intensity2, \
        np.mean([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, \
        1.0, 1.0, 1.0, 1.0, 1.0]))
        
    def test_getAverageIntensityInCircleOverTime(self):
    
        imMat = np.ones((10,10,3))
        imMat[5,5,:] = 100
        imMat[7,6,1] = 2
        imMat[7,6,2] = 3 
        intensities1 = ulfx.getAverageIntensityInCircleOverTime(imMat, (5,5), 1)
        intensities2 = ulfx.getAverageIntensityInCircleOverTime(imMat, (6,7), 1)
        intensities3 = ulfx.getAverageIntensityInCircleOverTime(imMat[:,:,0], (5,5), 1)
        i1 = np.mean([1.0, 1.0, 100.0, 1.0, 1.0])
        i2a = 1.0
        i2b =np. mean([1.0, 1.0, 2, 1.0, 1.0])
        i2c = np.mean([1.0, 1.0, 3, 1.0 ,1.0])

        self.assertEqual(intensities1, [i1, i1, i1])
        self.assertEqual(intensities2, [i2a, i2b, i2c])
        self.assertEqual(intensities3, [i1])
        
    def test_getDistancesToFirstPointInTrack(self):
        # set up some tracks for testing
        x1 = [10, 10, 10, 10, 10, 10]
        y1 = [20, 20, 20, 20, 20, 20]
        x2 = [70, 70, 70, 70, 70, 70]
        y2 = [10, 10, 10, 10, 10, 10]
        x3 = [50, 100, 100, 100, 100, 100]
        y3 = [50, 100, 100, 100, 100, 100]
        tracks = [(x1,y1),(x2,y2),(x3,y3)]

        dists = ulfx.getDistancesToFirstPointInTrack(tracks, (50,50))
        
        self.assertEqual(dists,[50.0, 44.721359549995796, 0.0])
        
    def test_normalizeIntensities(self):
        # test normalizeIntensities
        intensities = [[1,2,3,4,5],[100,200,300,400,500]]

        normedIntensities = ulfx.normalizeIntensities(intensities, 10)
        
        self.assertEqual(normedIntensities,[[0.1, 0.2, 0.3, 0.4, 0.5], \
        [10.0, 20.0, 30.0, 40.0, 50.0]] )
        
    def test_removeIntensitiesInConvertedRegion(self):
        intensities = [[1,2,3], [4,5,6], [7,8,9], [1,4,5]]
        distances = [1, 6, 7, 2]
        (newIntensities, newDistances) = \
        ulfx.removeIntensitiesInConvertedRegion(5, intensities, distances)
        self.assertEqual(newIntensities,[[4, 5, 6], [7, 8, 9]])
        self.assertEqual(newDistances, [6, 7])
        
    def test_findSlopes(self):       
        intensities = [[1,2,3,4,5],[100,200,300,400,500]]
        (dataSlopes, dataIntercepts, r_values) = ulfx.findSlopes(intensities)
        (dataSlopes2, dataIntercepts2, r_value2s) = ulfx.findSlopes(intensities, times = [0,2,4,6,8])
        self.assertEqual(dataSlopes, [1.0, 100.0])
        self.assertEqual(dataSlopes2, [0.5, 50.0])
        
    def test_getIntensities(self):
        data2 = np.ones((100,100,3))
        data2[10:20, 50:60, 0] = 3
        data2[20:30, 60:70, 1] = 4
        data2[30:40, 70:80, 2] = 5
        x1 = [55, 65, 75]
        y1 = [15, 25, 35]
        tracks1 = [(x1,y1)]
        x2 = [50, 51, 52]
        y2 = [30, 31,32]
        tracks2 = [(x2,y2)]
        
        intensities1 = ulfx.getIntensities(data2, tracks1, radius = 3)
        intensities2 = ulfx.getIntensities(data2, tracks2, radius = 3)
        
        self.assertEqual(intensities1, [[3.0, 4.0, 5.0]])
        self.assertEqual(intensities2, [[1.0, 1.0, 1.0]])
        
    def test_measureBleachingOverEntireMovie(self):
        data3 = np.ones((100,100,3))
        data3[:,:,1] = 2
        data3[:,:,2] = 3
        frameInts = ulfx.measureBleachingOverEntireMovie(data3)
        self.assertEqual(frameInts,[1.0, 2.0, 3.0])
        
    def test_getBackgroundULFInensities(self):
        # set up some tracks for testing
        x1 = [33, 10, 10, 10, 10, 10]
        y1 = [23, 20, 20, 20, 20, 20]
        x2 = [75, 70, 70, 70, 70, 70]
        y2 = [15, 10, 10, 10, 10, 10]
        x3 = [50, 100, 100, 100, 100, 100]
        y3 = [50, 100, 100, 100, 100, 100]
        tracks = [(x1,y1),(x2,y2),(x3,y3)]

        bk = np.ones((100,90))
        bk[23,33] = 100
        bk[15,75] =22
        bkInts = ulfx.getBackgroundULFInensities(bk, tracks)
        
        self.assertEqual(bkInts, [1.8761061946902655, 1.1858407079646018, 1.0])
        
    def test_subtractULFBackground(self):
        
        intensities = [[1,2,3], [10,10, 10], [4,4,4], [9, 8, 7]]
        ULFBackgrounds = [1,3,2,10]    
        correctedInts = ulfx.subtractULFBackground(ULFBackgrounds, \
                        intensities)
                        
        self.assertEqual(correctedInts, \
                         [[0, 1, 2], [7, 7, 7], [2, 2, 2], [0, 0, 0]])
            

suite = unittest.TestLoader().loadTestsFromTestCase(TestSequenceFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)