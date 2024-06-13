##
#  CMake script 
##

# Cmake version
CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0)

# Find dealii
FIND_PACKAGE(deal.II 9.3.0 HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../../ $ENV{HOME}/Desktop/d2)

IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

# Include madeal libraries
INCLUDE_DIRECTORIES(${MADEAL_DIRS} ../../include/)

# Set the target folder
SET(TARGET "fibre_composite_2d")

# Declare the main file
SET(TARGET_SRC main.cc)

# Define project
DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()
