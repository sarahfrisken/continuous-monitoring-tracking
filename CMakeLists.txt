cmake_minimum_required(VERSION 2.8.9)

project(ContinuousMonitoring)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/ContinuousMonitoring")
set(EXTENSION_CATEGORY "Examples")
set(EXTENSION_CONTRIBUTORS "Sarah Frisken (Radiology, Brigham and Women's Hospital)")
set(EXTENSION_DESCRIPTION "A module for experimental use in AMIGO at Brigham and Women's Hospital")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/ContinuousMonitoring.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/ContinuousMonitoring/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a space separated string, a list or 'NA' if any

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(CM_Ultrasound)
add_subdirectory(CM_Tracking)
add_subdirectory(CM_Tracking_Analysis)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
