import os
import math
import time
import unittest
import vtk, qt, ctk, slicer
import numpy
from slicer.ScriptedLoadableModule import *
import logging # NOTE: logging is slow!

#
# CM_Tracking_Analysis
#

class CM_Tracking_Analysis(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "CM_Tracking_Analysis"
    self.parent.categories = ["Continuous Monitoring"]
    self.parent.dependencies = []
    self.parent.contributors = ["Sarah Frisken (Radiology, Brighham and Women's Hospital)"]
    self.parent.helpText = """
    This is a scripted loadable module bundled into the Continuous Monitoring Extension.
    It supports analyzing and visualizing data from tracked of surgical instruments.
    """
    self.parent.helpText += self.getDefaultModuleDocumentationLink()
    self.parent.acknowledgementText = """
    This file was originally developed by Sarah Frisken, Radiology, BWH and was partially funded
    by NIH grant R01EB027134-01.
    """

#
# CM_TrackingAnalysisWidget
#
class CM_Tracking_AnalysisWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)
    # Instantiate and connect widgets ...

    #
    # Tracking Data
    #
    trackingCollapsibleButton = ctk.ctkCollapsibleButton()
    trackingCollapsibleButton.text = "Tracking Data"
    self.layout.addWidget(trackingCollapsibleButton)
    trackingFormLayout = qt.QFormLayout(trackingCollapsibleButton)

    # Loading point lists
    loadPointsButton = qt.QPushButton("Load Points")
    loadPointsButton.setToolTip("Load a new point list")
    loadPointsButton.connect('clicked(bool)', self.onLoadPointsButton)
    trackingFormLayout.addWidget(loadPointsButton)

    #
    # Resection Volume
    #
    resectionCollapsibleButton = ctk.ctkCollapsibleButton()
    resectionCollapsibleButton.text = "Resection Mapping"
    self.layout.addWidget(resectionCollapsibleButton)
    resectionLayout = qt.QVBoxLayout(resectionCollapsibleButton)

    # Resection type
    resTypeGroupBox = qt.QGroupBox("Type")
    resTypeLayout = qt.QVBoxLayout()
    self.resTypePointButton = qt.QRadioButton("Points")
    self.resTypePointButton.setChecked(True)
    self.resTypeLineButton = qt.QRadioButton("Lines")
    self.resTypeCurveButton = qt.QRadioButton("Curves")
    resTypeLayout.addWidget(self.resTypePointButton)
    resTypeLayout.addWidget(self.resTypeLineButton)
    resTypeLayout.addWidget(self.resTypeCurveButton)
    resTypeGroupBox.setLayout(resTypeLayout)
    resectionLayout.addWidget(resTypeGroupBox)

    # Time range for the resection
    resTimeRangeGroupBox = qt.QGroupBox("Time Range")
    resTimeRangeLayout = qt.QHBoxLayout()
    self.resStartTimeSpinBox = qt.QDoubleSpinBox()
    self.resStartTimeSpinBox.setRange(0.0, 1.0)
    self.resStartTimeSpinBox.setSingleStep(0.1)
    self.resEndTimeSpinBox = qt.QDoubleSpinBox()
    self.resEndTimeSpinBox.setRange(0.0, 1.0)
    self.resEndTimeSpinBox.setSingleStep(0.1)
    self.resEndTimeSpinBox.setValue(1.0)
    resTimeRangeLayout.addWidget(qt.QLabel("Start Time:"))
    resTimeRangeLayout.addWidget(self.resStartTimeSpinBox)
    resTimeRangeLayout.addStretch()
    resTimeRangeLayout.addWidget(qt.QLabel("End Time:"))
    resTimeRangeLayout.addWidget(self.resEndTimeSpinBox)
    resTimeRangeLayout.addStretch()
    resTimeRangeGroupBox.setLayout(resTimeRangeLayout)
    resectionLayout.addWidget(resTimeRangeGroupBox)

    # Resection volume name
    resFilesGroupBox = qt.QGroupBox("Files")
    resFilesLayout = qt.QFormLayout()
    self.resNameLineEdit = qt.QLineEdit("Resection")
    self.resNameLineEdit.setToolTip("Enter a name for the resection volume")
    resFilesLayout.addRow("Name Output Volume:", self.resNameLineEdit)

    # Mask volume for points contributing to resection
    self.resMaskSelector = slicer.qMRMLNodeComboBox()
    self.resMaskSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.resMaskSelector.selectNodeUponCreation = True
    self.resMaskSelector.addEnabled = False
    self.resMaskSelector.removeEnabled = False
    self.resMaskSelector.noneEnabled = True
    self.resMaskSelector.showHidden = False
    self.resMaskSelector.showChildNodeTypes = False
    self.resMaskSelector.setMRMLScene( slicer.mrmlScene )
    self.resMaskSelector.setToolTip( "Select a volume to mask resection data. Points where the value in the mask volume is zero will not be used." )
    resFilesLayout.addRow("Mask Volume:", self.resMaskSelector)
    resFilesGroupBox.setLayout(resFilesLayout)
    resectionLayout.addWidget(resFilesGroupBox)

    # Other resection volume parameters
    resParamGroupBox = qt.QGroupBox("Resection Parameters")
    resParamLayout = qt.QHBoxLayout()
    self.resVoxelSizeSpinBox = qt.QDoubleSpinBox()
    self.resVoxelSizeSpinBox.setValue(1.0)
    self.resToolRadiusSpinBox = qt.QDoubleSpinBox()
    self.resToolRadiusSpinBox.setValue(1.0)
    self.resPaddingSpinBox = qt.QSpinBox()
    resParamLayout.addWidget(qt.QLabel("Voxel Size:"))
    resParamLayout.addWidget(self.resVoxelSizeSpinBox)
    resParamLayout.addStretch()
    resParamLayout.addWidget(qt.QLabel("Tool Radius:"))
    resParamLayout.addWidget(self.resToolRadiusSpinBox)
    resParamLayout.addStretch()
    resParamLayout.addWidget(qt.QLabel("Volume Padding:"))
    resParamLayout.addWidget(self.resPaddingSpinBox)
    resParamLayout.addStretch()
    resParamGroupBox.setLayout(resParamLayout)
    resectionLayout.addWidget(resParamGroupBox)

    # Create resection volume
    createResVolButton = qt.QPushButton("Create Resection Volume")
    createResVolButton.setToolTip("Create a resection volume of the specified name and type from the current point list")
    createResVolButton.connect('clicked(bool)', self.onCreateResVolButton)
    resectionLayout.addWidget(createResVolButton)

    #
    # Analysis
    #
    self.analysisCollapsibleButton = ctk.ctkCollapsibleButton()
    self.analysisCollapsibleButton.text = "Analysis"
    self.layout.addWidget(self.analysisCollapsibleButton)
    analysisLayout = qt.QVBoxLayout(self.analysisCollapsibleButton)

    # Segmented resection volume 
    segMaskGroupBox = qt.QGroupBox("Input")
    segMaskLayout = qt.QFormLayout(segMaskGroupBox)
    self.segMaskSelector = slicer.qMRMLNodeComboBox()
    self.segMaskSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.segMaskSelector.selectNodeUponCreation = True
    self.segMaskSelector.addEnabled = False
    self.segMaskSelector.removeEnabled = False
    self.segMaskSelector.noneEnabled = True
    self.segMaskSelector.showHidden = False
    self.segMaskSelector.showChildNodeTypes = False
    self.segMaskSelector.setMRMLScene( slicer.mrmlScene )
    self.segMaskSelector.setToolTip( "Select the segmented resection volume mask. Analyze tracked points that intersect this mask vs. tracked points outside the mask." )
    segMaskLayout.addRow("Segmented Resection Volume:", self.segMaskSelector)
    analysisLayout.addWidget(segMaskGroupBox)

    # Perform the analysis
    analyzeButton = qt.QPushButton("Analyze")
    analyzeButton.setToolTip("Compare the tracked resection volume to the resection mask.")
    analysisLayout.addWidget(analyzeButton)
    analyzeButton.connect('clicked(bool)', self.onAnalyzeButton)

    # Results
    self.analysisResultsGroupBox = qt.QGroupBox("Results")
    self.analysisResultsGroupBox.setHidden(True)
    analysisResultsMaskLayout = qt.QFormLayout(self.analysisResultsGroupBox)
    self.analysisCoverage = qt.QLabel("")
    self.analysisOutside = qt.QLabel("")
    analysisResultsMaskLayout.addRow("Percent coverage:", self.analysisCoverage)
    analysisResultsMaskLayout.addRow("Percent points outside:", self.analysisOutside)
    analysisLayout.addWidget(self.analysisResultsGroupBox)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Set the initial state
    self.initialize()
    self.disableAnalysis()

  def cleanup(self):
    self.points = []
    self.resPoints = []
    self.pointsBounds = []
    if self.pointsHistogram:
      tableNodeID = self.pointsHistogram.GetNthPlotSeriesNode(0).GetTableNodeID()
      slicer.mrmlScene.RemoveNode(slicer.mrmlScene.GetNodeByID(tableNodeID))
      slicer.mrmlScene.RemoveNode(self.pointsHistogram.GetNthPlotSeriesNode(0))
      slicer.mrmlScene.RemoveNode(self.pointsHistogram)
    if self.trackedResDistVol:
      slicer.mrmlScene.RemoveNode(self.trackedResDistVol)

  def initialize(self):
    self.logic = CM_Tracking_AnalysisLogic()

    # Points
    self.points = []
    self.resPoints = []
    self.pointsBounds = []
    self.pointsHistogram = None
    self.pointsFile = None
    self.pointsDirectory = os.path.join(os.path.dirname(__file__))

    #  Resection volume
    self.trackedResDistVol = None

  def onLoadPointsButton(self):
    self.pointsFile = qt.QFileDialog.getOpenFileName(None, "Point Data File", self.pointsDirectory, '*.txt')
    if not self.pointsFile: return
    self.pointsDirectory = os.path.dirname(self.pointsFile)

    # read points from the points file
    self.points = self.readPoints(self.pointsFile)

    # process points and display point data
    self.pointsBounds = self.getPointsBounds(self.points)
    if self.pointsHistogram:
      tableNodeID = self.pointsHistogram.GetNthPlotSeriesNode(0).GetTableNodeID()
      slicer.mrmlScene.RemoveNode(slicer.mrmlScene.GetNodeByID(tableNodeID))
      slicer.mrmlScene.RemoveNode(self.pointsHistogram.GetNthPlotSeriesNode(0))
      slicer.mrmlScene.RemoveNode(self.pointsHistogram)
    self.pointsHistogram = self.getHistogramFromPoints(self.points)

  def readPoints(self, filename):
    # Read in a list of points. Each point p is an array of (time, x, y, z).
    points = []
    file = open(filename, 'r')
    for line in file:
      columns = line.split()
      p = [float(columns[1]), float(columns[3]), float(columns[4]), float(columns[5])]
      points.append(p)

    # scale times from 0 to 1
    pointsArray = numpy.array(points)
    times = pointsArray[:,0]
    t0 = numpy.amin(times)
    t1 = numpy.amax(times)
    offset = t0
    scale = 1.0 / (t1 - t0)
    for p in points:
      p[0] = (p[0] - offset) * scale
    return points

  def getPointsBounds(self, pointList):
    if len(pointList) == 0: return None
    xmin = pointList[0][1]
    ymin = pointList[0][2]
    zmin = pointList[0][3]
    xmax = xmin
    ymax = ymin
    zmax = zmin
    for point in pointList:
      if point[1] < xmin: xmin = point[1]
      elif point[1] > xmax: xmax = point[1]
      if point[2] < ymin: ymin = point[2]
      elif point[2] > ymax: ymax = point[2]
      if point[3] < zmin: zmin = point[3]
      elif point[3] > zmax: zmax = point[3]
    bounds = [xmin, ymin, zmin, xmax, ymax, zmax]
    return bounds
      
  def getHistogramFromPoints(self, pointList):
    # get numpy array of times
    times = [point[0] for point in pointList]
    timeArray = numpy.array(times)

    # crete numpy histogram from timeArray with N bins
    bins = 100
    histogram = numpy.histogram(timeArray, bins)

    # create a chart from the histogram that can be plotted in QT
    histogramChart = slicer.util.plot(histogram, xColumnIndex = 1)
    histogramChart.SetTitle("Points per time interval")
    histogramChart.GetNthPlotSeriesNode(0).SetPlotType(slicer.vtkMRMLPlotSeriesNode.PlotTypeScatterBar)
    return histogramChart
    
  def onCreateResVolButton(self):
    if len(self.points) == 0: 
      logging.info("Warning: Load tracked points before computing the resection volume.")
      return

    # Restrict points by time range and resection mask
    t0 = self.resStartTimeSpinBox.value
    t1 = self.resEndTimeSpinBox.value
    self.resPoints = self.getPointsInTimeRange(self.points, t0, t1)
    maskVol = self.resMaskSelector.currentNode()
    if maskVol:
      self.resPoints = self.getPointsInMask(self.resPoints, maskVol)
    if len(self.resPoints) == 0: 
      logging.info("Warning: No points in selected time range and resection mask.")
      return

    if self.trackedResDistVol:
      slicer.mrmlScene.RemoveNode(self.trackedResDistVol)
      self.trackedResDistVol = None
      
    voxelSize = self.resVoxelSizeSpinBox.value
    padding = self.resPaddingSpinBox.value
    name = self.resNameLineEdit.text
    self.trackedResDistVol = self.createVolume(self.pointsBounds, padding, voxelSize, name)

    # Add the resection volume to the scene
    slicer.mrmlScene.AddNode(self.trackedResDistVol)
    displayNode = slicer.vtkMRMLScalarVolumeDisplayNode()
    slicer.mrmlScene.AddNode(displayNode)
    colorNode = slicer.util.getNode('Grey')
    displayNode.SetAndObserveColorNodeID(colorNode.GetID())
    self.trackedResDistVol.SetAndObserveDisplayNodeID(displayNode.GetID())

    # Compute the binary resection volume by adding the resection points 
    # to the empty volume. (Note this only works after the volume has been 
    # added to the scene) Then compute the distance field of the resection volume.
    self.addPointsToVolume(self.trackedResDistVol, self.resPoints)
    self.logic.convertToDistanceField(self.trackedResDistVol)

    # Enable analysis of the resection volume 
    self.enableAnalysis()
    
  def createVolume(self, bounds, pad, voxelSize, name = 'volume'):
    origin = [bounds[0] - pad, bounds[1] - pad, bounds[2] - pad]
    spacing = [voxelSize, voxelSize, voxelSize]
    size = [
      (bounds[3] - bounds[0] + 2 * pad) / voxelSize + 1, 
      (bounds[4] - bounds[1] + 2 * pad) / voxelSize + 1, 
      (bounds[5] - bounds[2] + 2 * pad) / voxelSize + 1 ]
    if (size[0] is 0) or (size[1] is 0) or (size[2] is 0) :
      return

    # Create an empty volume that contains all the points
    imageData = vtk.vtkImageData()
    imageData.SetDimensions(int(math.ceil(size[0])), int(math.ceil(size[1])), int(math.ceil(size[2])))
    imageData.AllocateScalars(vtk.VTK_SHORT, 1)
    thresholder=vtk.vtkImageThreshold()
    thresholder.SetInputData(imageData)
    thresholder.SetInValue(0)
    thresholder.SetOutValue(0)
    volume = slicer.vtkMRMLScalarVolumeNode()
    volume.SetOrigin(origin)
    volume.SetSpacing(spacing)
    volume.SetName(name)
    volume.CreateDefaultStorageNode()
    volume.SetImageDataConnection(thresholder.GetOutputPort())
    return volume

  def addPointsToVolume(self, volume, pointList):
    # Get transforms between world and volume coordinates
    matrixRASToIJK = vtk.vtkMatrix4x4()
    volume.GetRASToIJKMatrix(matrixRASToIJK)
    spacing = volume.GetSpacing()
    imageData=volume.GetImageData()
    volDims = imageData.GetDimensions()
    volArray = slicer.util.arrayFromVolume(volume)

    # Add points and point connections depending on fill type
    fillType = "Points"
    if self.resTypeLineButton.isChecked() == True: 
      fillType = "Lines"
    elif self.resTypeCurveButton.isChecked() == True: 
      fillType = "Curves"

    prevPrevIJK = None
    prevIJK = None
    for point in pointList:
      ras = [point[1], point[2], point[3], 1]
      pointIJK = matrixRASToIJK.MultiplyPoint(ras)
      if fillType == "Points":
        self.addPointToDistArray(pointIJK, volArray)
      elif fillType == "Lines":
        self.addLineToDistArray(prevIJK, pointIJK, volArray)
        prevIJK = pointIJK

  def getPointsInTimeRange(self, pointList, t0, t1):
    # Returns a copy of pointList without points outside [t0,t1]
    points = []
    for point in pointList:
      time = point[0]
      if (time >= t0) and (time <= t1):
        points.append(point)
    return points

  def getPointsInMask(self, pointList, maskVol):
    # Returns a copy of pointList without points that do not lie inside non-zero voxels of maskVol
    maskRASToIJK = vtk.vtkMatrix4x4()
    maskVol.GetRASToIJKMatrix(maskRASToIJK)
    maskDims = maskVol.GetImageData().GetDimensions()
    maskArray = slicer.util.arrayFromVolume(maskVol)
    points = []
    for point in pointList:
      ras = [point[1], point[2], point[3], 1]
      pMaskIJK = maskRASToIJK.MultiplyPoint(ras)
      if self.isPointInMask(pMaskIJK, maskArray, maskDims):
        points.append(point)
    return points

  def isPointInMask(self, p, mask, maskDims):
    if (p[0] < 0) or (p[0] > maskDims[0]-1) or \
       (p[1] < 0) or (p[1] > maskDims[1]-1) or \
       (p[2] < 0) or (p[2] > maskDims[2]-1): 
       return False
    maskValue = mask[int(p[2])][int(p[1])][int(p[0])]
    if maskValue == 0: return False
    else: return True

  def addPointToDistArray(self, p, distArray):
    # Assumes that the point is inside the volume (no bounds checking -- for speed)
    distArray[int(p[2])][int(p[1])][int(p[0])] = 1  

  def addLineToDistArray(self, prev, p, distArray):
    # Assumes that both points are inside the volume (no bounds checking -- for speed)
    if prev == None:
      return

    # Use the DDA (digital differential analyzer) algorithm to draw the line in the volume
    dx = p[0] - prev[0]
    dy = p[1] - prev[1]
    dz = p[2] - prev[2]
    xLen = abs(dx)
    yLen = abs(dy)
    zLen = abs(dz)
    if xLen == 0 and yLen == 0 and zLen == 0:
      return
    if xLen >= yLen and xLen >= zLen:
      len = xLen
    elif yLen >= xLen and yLen >= zLen:
      len = yLen
    elif zLen >= xLen and zLen >= yLen:
      len = zLen
    dx = dx / len
    dy = dy / len
    dz = dz / len

    x = prev[0]
    y = prev[1]
    z = prev[2]
    step = 0
    while step <= len:
      distArray[int(z)][int(y)][int(x)] = 1
      x = x + dx
      y = y + dy
      z = z + dz
      step = step + 1

  def onAnalyzeButton(self):
    # Get the segmentation mask and other analysis parameters
    segMaskVol = self.segMaskSelector.currentNode()
    numBins = 10 # More bins are not meaningful for the approximate distance function used

    # Compute coverage of the interior of the segmentation mask by the tracked resection
    coverage = CM_Coverage(self.trackedResDistVol, segMaskVol, numBins)
    coverage.commaDelimitedOutput()

    # Compute overlap of the tracked resection with the exterior of the segmentation mask
    overlap = CM_Overlap(self.resPoints, segMaskVol, numBins)
    overlap.commaDelimitedOutput()

    # Display the output of the analysis 
    toolRadius = self.resToolRadiusSpinBox.value
    pctCoverage = coverage.getPercentCoverage(toolRadius)
    self.analysisCoverage.setText(str(pctCoverage))
    pctOverlap = overlap.getPercentOverlap(toolRadius)
    self.analysisOutside.setText(str(pctOverlap))
    self.analysisResultsGroupBox.setHidden(False)

  def onChangeParameter(self):
    self.disableAnalysis()

  def enableAnalysis(self):
    self.analysisCoverage.setText("")
    self.analysisOutside.setText("")
    self.analysisResultsGroupBox.setHidden(True)
    self.analysisCollapsibleButton.collapsed = False
    self.analysisCollapsibleButton.setEnabled(True)
  
  def disableAnalysis(self):
    self.analysisResultsGroupBox.setHidden(True)
    self.analysisCollapsibleButton.collapsed = True
    self.analysisCollapsibleButton.setEnabled(False)

#
# CM_Coverage
#
class CM_Coverage():
  # Computes coverage of the specified labelmap by a shape defined by the given distance field,
  # where distances are negative inside and positive outside the shape. Records coverage by 
  # counting the number of non-zero voxels in the labelmap that lie *outside* the shape. Points 
  # closer than binMinDist to shape are put in bin[0].
  def __init__(self, distVolume, labelmap, numBins = 10, binMinDist = 0, binMaxDist = 10.0):
    self.numNonzeroVoxels = 0
    self.numVoxelsInBin = []

    if labelmap is None :
      logging.info("Cannot compute coverage. Missing labelmap volume")
      return
    else:
      labelmapIJKToRAS = vtk.vtkMatrix4x4()
      labelmap.GetIJKToRASMatrix(labelmapIJKToRAS)
      labelmapDims = labelmap.GetImageData().GetDimensions()
      labelmapArray = slicer.util.arrayFromVolume(labelmap)

    if distVolume is None :
      logging.info("Cannot compute coverage. Missing distance volume")
      return
    else :
      distRASToIJK = vtk.vtkMatrix4x4()
      distVolume.GetRASToIJKMatrix(distRASToIJK)
      distDims = distVolume.GetImageData().GetDimensions()
      distArray = slicer.util.arrayFromVolume(distVolume)

    # Initialize the binned points
    self.binMinDist = binMinDist
    self.binMaxDist = binMaxDist
    numBins = 1 + min(10000, max(1, numBins))
    self.binSize = (binMaxDist - binMinDist) / (numBins - 1)
    self.numVoxelsInBin = numpy.zeros(numBins)

    ras = [0,0,0,1]
    distIJK = [0,0,0,1]
    self.numSegVoxels = 0
    for k in range (0, labelmapDims[2]-1):
      for j in range(0, labelmapDims[1]-1):
        for i in range(0, labelmapDims[0]-1):
          if labelmapArray[k][j][i] > 0 :
            self.numNonzeroVoxels = self.numNonzeroVoxels + 1

            # Check the distance to the shape at this voxel
            labelmapIJKToRAS.MultiplyPoint([float(i), float(j), float(k), 1.0], ras)
            distRASToIJK.MultiplyPoint(ras, distIJK)
            di = int(distIJK[0])
            dj = int(distIJK[1])
            dk = int(distIJK[2])
            if (di>0) and (di<distDims[0]) and (dj>0) and (dj<distDims[1]) and (dk>0) and (dk<distDims[2]): 
              dist = distArray[dk][dj][di]
              if dist < self.binMinDist:
                self.numVoxelsInBin[0] = self.numVoxelsInBin[0] + 1
              else:
                bin = min(1 + int((dist - self.binMinDist) / self.binSize), numBins - 1)
                self.numVoxelsInBin[bin] = self.numVoxelsInBin[bin] + 1

  def __str__(self):
    # For printing this class
    return  str(self.__class__) + '\n'+ '\n'.join(('  {} = {}'.format(item, self.__dict__[item]) for item in self.__dict__))
 
  def getPercentCoverage(self, maxDist):
    if self.numNonzeroVoxels is 0:
      logging.info("Warning: no non-zero voxels in labelmap")
      return 0

    maxBin = 1 + int((maxDist - self.binMinDist) / self.binSize)
    numVoxelsCovered = 0
    for i in range (0, maxBin):
      numVoxelsCovered = numVoxelsCovered + self.numVoxelsInBin[i]
    return 100 * numVoxelsCovered / self.numNonzeroVoxels

  def commaDelimitedOutput(self):
    labelStr = "binMin, binMax, binSize, total #voxels, "
    valueStr = str(self.binMinDist) + ", " + str(self.binMaxDist) + ", " + str(self.binSize) + ", " + str(self.numNonzeroVoxels) + ", "
    for i in range(0, len(self.numVoxelsInBin)) :
      labelStr = labelStr + "#voxels[" + str(i) + "], "
      valueStr = valueStr + str(self.numVoxelsInBin[i]) + ", "
    logging.info(labelStr)
    logging.info(valueStr)


#
# CM_Overlap
#
class CM_Overlap():
  # Computes the overlap between points in a point list and the *outside* of the given labelmap. 
  # Records the overlap by counting the number of points that lie within binned distance ranges. 
  # Points closer than binMinDist to labeled regions are put in bin[0].
  def __init__(self, pointList, labelmap, numBins = 10, binMinDist = 0, binMaxDist = 10.0):
    self.numPoints = 0
    self.numPointsInBin = []

    if labelmap is None :
      logging.info("Cannot compute overlap. Missing labelmap volume")
      return

    # Create a distance map from the labelmap volume. Distance values will be positive outside 
    # non-zero pixels in the labelmap and negative inside.
    logic = CM_Tracking_AnalysisLogic()
    distVol = slicer.modules.volumes.logic().CloneVolume(slicer.mrmlScene, labelmap, "Distance Volume")
    logic.convertToDistanceField(distVol)

    # Initialize the binned points
    self.binMinDist = binMinDist
    self.binMaxDist = binMaxDist
    numBins = 1 + min(10000, max(1, numBins))
    self.binSize = (binMaxDist - binMinDist) / (numBins - 1)
    self.numPointsInBin = numpy.zeros(numBins)

    # Count points outside the labelmap and bin the results by distance to the closest non-zero label
    distRASToIJK = vtk.vtkMatrix4x4()
    distVol.GetRASToIJKMatrix(distRASToIJK)
    distDims = distVol.GetImageData().GetDimensions()
    distArray = slicer.util.arrayFromVolume(distVol)
    self.numPoints = len(pointList)
    for point in pointList:
      ras = [point[1], point[2], point[3], 1]
      distIJK = distRASToIJK.MultiplyPoint(ras)
      di = int(distIJK[0])
      dj = int(distIJK[1])
      dk = int(distIJK[2])
      if (di>0) and (di<distDims[0]) and (dj>0) and (dj<distDims[1]) and (dk>0) and (dk<distDims[2]): 
        dist = distArray[dk][dj][di]
        if dist < self.binMinDist:
          self.numPointsInBin[0] = self.numPointsInBin[0] + 1
        else:
          bin = min(1 + int((dist - self.binMinDist) / self.binSize), numBins - 1)
          self.numPointsInBin[bin] = self.numPointsInBin[bin] + 1
   
    # Clean up
    slicer.mrmlScene.RemoveNode(distVol)

  def __str__(self):
    # For printing this class
    return  str(self.__class__) + '\n'+ '\n'.join(('  {} = {}'.format(item, self.__dict__[item]) for item in self.__dict__))
 
  def getPercentOverlap(self, maxDist):
    if self.numPoints is 0:
      logging.info("Warning: no points in overlap")
      return 0

    maxBin = 1 + int((maxDist - self.binMinDist) / self.binSize)
    numPointsInside = 0
    for i in range (0, maxBin):
      numPointsInside = numPointsInside + self.numPointsInBin[i]
    return 100 * (self.numPoints - numPointsInside) / self.numPoints

  def commaDelimitedOutput(self):
    labelStr = "binMin, binMax, binSize, total #points, "
    valueStr = str(self.binMinDist) + ", " + str(self.binMaxDist) + ", " + str(self.binSize) + ", " + str(self.numPoints) + ", "
    for i in range(0, len(self.numPointsInBin)) :
      labelStr = labelStr + "#points[" + str(i) + "], "
      valueStr = valueStr + str(self.numPointsInBin[i]) + ", "
    logging.info(labelStr)
    logging.info(valueStr)
  
#
# CM_TrackingAnalysisLogic
#
class CM_Tracking_AnalysisLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual computation done by your module.  
  The interface should be such that other python code can import this class and
  make use of the functionality without requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  def convertToDistanceField(self, volume):
    # The SignedMaurerDistanceMapImageFilter returns an image whose voxels contain the
    # distance to the closest non-background voxel. By default, background voxels are zero.
    import SimpleITK as sitk
    import sitkUtils
    image = sitk.ReadImage(sitkUtils.GetSlicerITKReadWriteAddress(volume))
    f = sitk.SignedMaurerDistanceMapImageFilter()
    f.SetSquaredDistance(False)
    f.SetInsideIsPositive(False)
    f.SetUseImageSpacing(True)
    distImage = f.Execute(image)
    sitk.WriteImage(distImage, sitkUtils.GetSlicerITKReadWriteAddress(volume))
    volume.GetImageData().Modified()



class CM_Tracking_AnalysisTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test1()

  def test1(self):
    self.delayDisplay("Starting the test")
    self.delayDisplay('Test passed!')
