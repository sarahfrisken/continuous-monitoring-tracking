import os
import math
import time
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging # NOTE: logging is slow!

#
# CM_Tracking_Simple
#

class CM_Tracking_Simple(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "CM_Tracking_Simple"
    self.parent.categories = ["Continuous Monitoring"]
    self.parent.dependencies = []
    self.parent.contributors = ["Sarah Frisken (Radiology, Brighham and Women's Hospital)"]
    self.parent.helpText = """
This is a scripted loadable module bundled into the Continuous Monitoring Extension.
It supports the tracking of surgical instruments and writing data to a file as it is collected.
"""
    self.parent.helpText += self.getDefaultModuleDocumentationLink()
    self.parent.acknowledgementText = """
This file was originally developed by Sarah Frisken, Radiology, BWH and was partially funded
by NIH grant R01EB027134-01.
"""

#
# CM_Tracking_SimpleWidget
#
class CM_Tracking_SimpleWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...

    #
    # Tracked Instruments
    #
    trackedCollapsibleButton = ctk.ctkCollapsibleButton()
    trackedCollapsibleButton.text = "Tracked Instruments"
    self.layout.addWidget(trackedCollapsibleButton)
    trackedFormLayout = qt.QFormLayout(trackedCollapsibleButton)

    #
    # Instrument list
    #
    self.instrumentsLayout = qt.QVBoxLayout()
    trackedFormLayout.addRow(self.instrumentsLayout)

    # New instrument button
    self.newButton = qt.QPushButton("New Instrument")
    self.newButton.setToolTip("Set up data collection for a new instrument")
    self.newButton.connect('clicked(bool)', self.onNewInstrumentButton)
    trackedFormLayout.addWidget(self.newButton)
    
    #
    # Saving and Loading
    #
    # - set save directory
    #   - consider auto saving -- what happens when save directory is changed?
    # - save/load volume, ROI, instruments, tracking data
    # - when load, consider automatically generating resection volumes
    saveLoadCollapsibleButton = ctk.ctkCollapsibleButton()
    saveLoadCollapsibleButton.text = "Saving and Loading"
    self.layout.addWidget(saveLoadCollapsibleButton)
    saveLoadLayout = qt.QVBoxLayout(saveLoadCollapsibleButton)

    # Save button
    self.saveButton = qt.QPushButton('Save')
    self.saveButton.setToolTip("Save the instrument data to respective save files")
    self.saveButton.connect('clicked(bool)', self.onSaveInstrumentData)
    saveLoadLayout.addWidget(self.saveButton)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Set the initial state
    self.initialize()

  def cleanup(self):
    pass

  def initialize(self):
    self.saveDirectory = os.path.join(os.path.dirname(__file__))
    self.instruments = []

  def onNewInstrumentButton(self):
    d = CM_NewInstrumentDialog()
    d.setModal(True)
    d.saveDirectory = self.saveDirectory
    if d.exec_():
      if d.instrument:
        self.instruments.append(d.instrument)
        self.addInstrumentPanel(d.instrument) 
      if d.saveFile:
        self.saveFile = d.saveFile
        self.saveDirectory = os.path.dirname(d.saveFile)

  def onSaveInstrumentData(self):
    for instrument in self.instruments:
      if instrument.dataFile:
        instrument.dataFile.close()
        instrument.dataFile = open(instrument.dataFilePath, 'a')

  def addInstrumentPanel(self, instrument): 
    panelLayout = qt.QHBoxLayout()
    setActiveButton = instrument.setActiveButton
    setActiveButton.setChecked(instrument.isActive)
    setActiveButton.toggled.connect(lambda:self.activateInstrument(instrument))
    panelLayout.addWidget(setActiveButton)
    panelLayout.addWidget(instrument.locationLabel)
    panelLayout.addWidget(instrument.saveFileLabel)
    panelLayout.addStretch(1)
    self.instrumentsLayout.addLayout(panelLayout)
    
  def activateInstrument(self, instrument):
    instrument.setActive(instrument.setActiveButton.isChecked())
    logging.info("instrument = {}, buttonState = {}, isActive = {}".format(
      instrument.transform.GetName(), instrument.setActiveButton.isChecked(), instrument.isActive))

#
# CM_NewInstrumentDialog
#
class CM_NewInstrumentDialog(qt.QDialog):
  def __init__(self, parent = None):
    super(CM_NewInstrumentDialog, self).__init__(parent)

    # Initialize variables
    self.instrument = None
    self.saveFile = None
    self.saveDirectory = None

    layout = qt.QVBoxLayout(self)

    #
    form = qt.QGroupBox()
    formLayout = qt.QFormLayout()
    form.setLayout(formLayout)
    layout.addWidget(form)

    # Add tracked instruments
    self.instrumentSelector = slicer.qMRMLNodeComboBox()
    self.instrumentSelector.nodeTypes = ["vtkMRMLLinearTransformNode"]
    self.instrumentSelector.selectNodeUponCreation = True
    self.instrumentSelector.addEnabled = False
    self.instrumentSelector.removeEnabled = False
    self.instrumentSelector.noneEnabled = True
    self.instrumentSelector.showHidden = False
    self.instrumentSelector.showChildNodeTypes = False
    self.instrumentSelector.setMRMLScene( slicer.mrmlScene )
    self.instrumentSelector.setToolTip( "Add an instrument to track" )
    formLayout.addRow("Instrument Transform: ", self.instrumentSelector)

    # Define volume of interest for tracking
    self.roiSelector = slicer.qMRMLNodeComboBox()
    self.roiSelector.nodeTypes = ["vtkMRMLAnnotationROINode"]
    self.roiSelector.selectNodeUponCreation = True
    self.roiSelector.addEnabled = True
    self.roiSelector.baseName = "Tracking ROI"
    self.roiSelector.removeEnabled = True
    self.roiSelector.noneEnabled = False
    self.roiSelector.showHidden = False
    self.roiSelector.showChildNodeTypes = False
    self.roiSelector.setMRMLScene( slicer.mrmlScene )
    self.roiSelector.setToolTip( 
      "Pick a region of interest to bound tracked data. "
      "The current ROI defines the size of the tracking volume "
      "for a particular tracked instrument when it is created." )
    formLayout.addRow("Tracking Volume: ", self.roiSelector)

    #
    buttons = qt.QGroupBox()
    buttonsLayout = qt.QVBoxLayout()
    buttons.setLayout(buttonsLayout)
    layout.addWidget(buttons)

    # Define the data file for this instrument
    setSaveFileButton = qt.QPushButton("Set Save Filename")
    setSaveFileButton.setToolTip("Set the name of the file for saving instrument data")
    setSaveFileButton.connect('clicked(bool)', self.onSetSaveFileButton)
    self.saveFileTextLabel = qt.QLabel("Save File: ")
    buttonsLayout.addWidget(setSaveFileButton)
    buttonsLayout.addWidget(self.saveFileTextLabel)

    # Create the tracked instrument object
    createInstrumentButton = qt.QPushButton("Create Instrument")
    createInstrumentButton.setToolTip("Create the instrument")
    createInstrumentButton.connect('clicked(bool)', self.onCreateInstrumentButton)
    buttonsLayout.addWidget(createInstrumentButton)

  def onSetSaveFileButton(self):
    self.saveFile = qt.QFileDialog.getSaveFileName(None, "Instrument Data File", 
      self.saveDirectory, '*.txt')
    self.saveFileTextLabel.setText("Save File: " + self.saveFile)

  def onCreateInstrumentButton(self):
    transform = self.instrumentSelector.currentNode()
    roi = self.roiSelector.currentNode()
    if self.saveFile and transform and roi :
      self.instrument = CM_TrackingDevice(transform, roi, self.saveFile)
      if self.instrument:
        self.accept()
    else :
      self.showWarningMessage()

  def showWarningMessage(self):
    msg = qt.QMessageBox()
    msg.setIcon(qt.QMessageBox.Warning)
    msg.setText("Please set the instrument, tracking volume and save file")
    msg.setWindowTitle("Warning")
    msg.setStandardButtons(qt.QMessageBox.Ok)
    msg.exec_()

#
# CM_TrackingDevice
#
class CM_TrackingDevice():
  def __init__(self, trackingTransform, roi, savefile):
    # initialize the device
    self.doEditVolume = False #True #False
    self.transform = trackingTransform
    self.isActive = False
    self.setActiveButton = qt.QCheckBox(self.transform.GetName())
    self.locationLabel = qt.QLabel()
    self.saveFileLabel = qt.QLabel(savefile)
    self.observerTag = None
    self.dataFilePath = savefile
    self.dataFile = open(savefile, 'w')
    self.dataFile.close()
    self.dataFile = None
    voxelSize = [0.5, 0.5, 0.5]
    self.createDataVolume(roi, voxelSize, trackingTransform.GetName())

  def __del__(self):
    self.setActive(False)

  def setActive(self, isActive): 
    if isActive is True:
      self.dataFile = open(self.dataFilePath, 'a')
      self.observerTag = self.transform.AddObserver(
      slicer.vtkMRMLTransformNode.TransformModifiedEvent, self.onDeviceMovedEvent)
    else:
      self.transform.RemoveObserver(self.observerTag)
      self.observerTag = None
      self.dataFile.close()
      self.dataFile = None
    self.isActive = isActive

  def createDataVolume(self, roi, voxelSize, name):
    self.dataVolume = None
    if roi is None:
      logging.info('No data volume created: ROI not specified')
      return

    # Determine volume parameters
    dataExtent = [0, -1, 0, -1, 0, -1]
    roi.GetBounds(dataExtent)
    logging.info('data extent = {}'.format(dataExtent))
    origin = [dataExtent[0], dataExtent[2], dataExtent[4]]
    spacing = voxelSize
    size = [
      (dataExtent[1] - dataExtent[0]) / voxelSize[0] + 1, 
      (dataExtent[3] - dataExtent[2]) / voxelSize[1] + 1, 
      (dataExtent[5] - dataExtent[4]) / voxelSize[2] + 1 ]
    if (size[0] is 0) or (size[1] is 0) or (size[2] is 0) :
      logging.info('No data volume created: ROI size is zero')
      return
    logging.info('Create volume with origin {}, size {} and voxel size {}'.format(origin, size, spacing))

     # Create an empty volume that contains the ROI
    imageData = vtk.vtkImageData()
    imageData.SetDimensions(int(math.ceil(size[0])), int(math.ceil(size[1])), int(math.ceil(size[2])))
    imageData.AllocateScalars(vtk.VTK_UNSIGNED_CHAR, 1)
    thresholder=vtk.vtkImageThreshold()
    thresholder.SetInputData(imageData)
    thresholder.SetInValue(0)
    thresholder.SetOutValue(0)
    self.dataVolume = slicer.vtkMRMLScalarVolumeNode()
    self.dataVolume.SetOrigin(origin)
    self.dataVolume.SetSpacing(spacing)
    self.dataVolume.SetName(name + 'Volume')
    self.dataVolume.CreateDefaultStorageNode()
    self.dataVolume.SetImageDataConnection(thresholder.GetOutputPort())

    # Add it to the scene
    slicer.mrmlScene.AddNode(self.dataVolume)
    displayNode = slicer.vtkMRMLScalarVolumeDisplayNode()
    slicer.mrmlScene.AddNode(displayNode)
    colorNode = slicer.util.getNode('Grey')
    displayNode.SetAndObserveColorNodeID(colorNode.GetID())
    self.dataVolume.SetAndObserveDisplayNodeID(displayNode.GetID())

    # eventually this is a UI checkbox
    self.setActive(True)

  def onDeviceMovedEvent(self, caller, event):
    if self.isActive is False:
      return
      
    p = [0.0,0.0,0.0]
    isInROI = False
    if self.transform:
      matrix = self.transform.GetMatrixTransformToParent()
      p[0] = matrix.GetElement(0, 3)
      p[1] = matrix.GetElement(1, 3)
      p[2] = matrix.GetElement(2, 3)
      time = matrix.GetMTime()
      if self.dataVolume:
        if self.addPointToVolume(p[0], p[1], p[2], self.doEditVolume):
          isInROI = True
          self.dataFile.write("time: {} location: {} {} {}\n".format(time, p[0], p[1], p[2]))

    # Update the GUI of the instrument location. Text is blue if the instrument
    # is inside the ROI
    positionText = "({p_x:3.1f}, {p_y:3.1f}, {p_z:3.1f}) ".format(p_x= p[0], p_y= p[1], p_z= p[2])
    self.locationLabel.text = positionText
    if isInROI is True:
      self.locationLabel.setStyleSheet('color: blue')
    else:
      self.locationLabel.setStyleSheet('color: black')
    
  def addPointToVolume(self, x, y, z, doEditVolume):
    # keep track of intersection with volue
    isPointInVolume = False

    # Filter constants
    minValue = 0
    maxValue = 255
    toolRadius = 1.0
    filterRadius = 2.0
    filterScale = maxValue / (2 * filterRadius)

    # Get transforms between world and volume coordinates
    matrixRASToIJK = vtk.vtkMatrix4x4()
    self.dataVolume.GetRASToIJKMatrix(matrixRASToIJK)
    
    # Transform (x,y,z) to image coordinates and get closest raster position
    ras = [x, y, z, 1]
    ijk = matrixRASToIJK.MultiplyPoint(ras)
    i = int(ijk[0])
    j = int(ijk[1])
    k = int(ijk[2])
    
    # Get bounds around (i, j, k)
    # NEED TO DO BETTER -- take spacing into account
    imageData=self.dataVolume.GetImageData()
    dims = imageData.GetDimensions()
    filterExtent = toolRadius + filterRadius
    iMin = int(max(0, i-filterExtent))
    iMax = int(min(dims[0]-1, i+filterExtent))
    jMin = int(max(0, j-filterExtent))
    jMax = int(min(dims[1]-1, j+filterExtent))
    kMin = int(max(0, k-filterExtent))
    kMax = int(min(dims[2]-1, k+filterExtent))

    # Update the distance field near (i, j, k)
    for k in range(kMin, kMax):
      for j in range(jMin, jMax):
        for i in range(iMin, iMax):
          isPointInVolume = True
          if doEditVolume:
            distToPt = toolRadius - math.sqrt((ijk[0]-i)**2 + (ijk[1]-j)**2 + (ijk[2]-k)**2)
            distToPt = (distToPt + filterRadius) * filterScale
            distToPt = max(min(distToPt, maxValue), minValue)
            dist = imageData.GetScalarComponentAsFloat(i, j, k, 0)            # SLOW, pointer access to data would be better, but can't be done in Python
            if dist < distToPt:
              imageData.SetScalarComponentFromFloat(i, j, k, 0, distToPt)  # SLOW, pointer access to data would be better, but can't be done in Python
    if isPointInVolume and doEditVolume:
      imageData.Modified()
    return isPointInVolume

#
# CM_Tracking_SimpleLogic
#
class CM_Tracking_SimpleLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual computation done by your module.  
  The interface should be such that other python code can import this class and
  make use of the functionality without requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def hasImageData(self,volumeNode):
    """This is an example logic method that returns true if the passed in volume
    node has valid image data
    """
    if not volumeNode:
      logging.debug('hasImageData failed: no volume node')
      return False
    if volumeNode.GetImageData() is None:
      logging.debug('hasImageData failed: no image data in volume node')
      return False
    return True


class CM_Tracking_SimpleTest(ScriptedLoadableModuleTest):
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
    self.test_CM_Tracking1()

  def test_CM_Tracking1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import urllib
    downloads = (
        ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
        )

    for url,name,loader in downloads:
      filePath = slicer.app.temporaryPath + '/' + name
      if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
        logging.info('Requesting download %s from %s...\n' % (name, url))
        urllib.urlretrieve(url, filePath)
      if loader:
        logging.info('Loading %s...' % (name,))
        loader(filePath)
    self.delayDisplay('Finished with download and loading')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = CM_Tracking_SimpleLogic()
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
