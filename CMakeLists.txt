# $Id$
# 
# GMAT: General Mission Analysis Tool.
# 
# CMAKE script file for GMAT Project
# This file must be installed in the main GMAT directory.
# That is, we should have the directory structure:
#   ./src
#   ./plugins
#   ./depends
#   ./application
#   etc...
#
# Original Author: Joris T. Olympio
# Modified For Production By: Josh Fisher
# Modified for CMake standardization by: Ravi Mathur
#
# DO NOT MODIFY THIS FILE UNLESS YOU KNOW WHAT YOU ARE DOING!
#

# Require CMake minimum version
cmake_minimum_required(VERSION 3.14.0)
CMAKE_POLICY(VERSION 3.14)

# Require OSX minimum version
IF(APPLE)
  SET(CMAKE_OSX_DEPLOYMENT_TARGET 10.15 CACHE STRING
    "Minimum OSX version to target for deployment (at runtime). Set to empty string to target current MacOS version.")
ENDIF()

# Set Cmake project
PROJECT(GMAT C CXX)
SET(GMAT_RELEASE_NAME "R2022a" CACHE STRING "GMAT version name")
SET(GMAT_VERSION ${GMAT_RELEASE_NAME})
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)

# Enforce C++ Standard. C++17 is preferred.
SET(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ Standard used for the build. C++17 is preferred")
SET_PROPERTY(CACHE CMAKE_CXX_STANDARD PROPERTY STRINGS 17 14 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_EXTENSIONS OFF)

# Path to custom Find****.cmake files
SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/build/cmake_modules")

# Detect 32/64 bit compiler
if(CMAKE_CL_64 OR (CMAKE_SIZEOF_VOID_P EQUAL 8))
  SET(GMAT_64_BIT ON)
  MESSAGE(STATUS "*** Creating 64-bit GMAT Build System. Please ensure that 64-bit versions of all dependencies are available. ***")
else()
  SET(GMAT_64_BIT OFF)
  MESSAGE(STATUS "*** Creating 32-bit GMAT Build System. Please ensure that 32-bit versions of all dependencies are available. ***")
endif()

# Default to Release configuration on single-config generators
# e.g. make, but not VisualStudio or XCode
# Note that CMAKE_BUILD_TYPE is initialized to "" on the initial cmake-gui
# run, so we check for that and set the default appropriately
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "Build Configuration (Release, Debug, RelWithDebInfo, or MinSizeRel" FORCE)
else()
  SET(CMAKE_BUILD_TYPE Debug CACHE STRING "Build Configuration (Release, Debug, RelWithDebInfo, or MinSizeRel" FORCE)
endif()

# Set default GMAT installation location
# Can be overridden by using "-DCMAKE_INSTALL_PREFIX=/foo/bar" at command line
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if(GMAT_64_BIT)
    SET(GMAT_RELEASE_TYPE "x64")
  else()
    SET(GMAT_RELEASE_TYPE "x86")
  endif()

  SET(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/GMAT-${GMAT_RELEASE_NAME}-${CMAKE_SYSTEM_NAME}-${GMAT_RELEASE_TYPE}" CACHE PATH "GMAT top-directory installation location" FORCE )
endif()

# Set path to app bundle on Mac, relative to CMAKE_INSTALL_PREFIX
# This may need to be a user-selectable cache variable in the future
if(APPLE)
  if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    SET(GMAT_MAC_APPBUNDLE_PATH "bin/GMATd-${GMAT_RELEASE_NAME}_Beta.app/Contents")
  else()
    SET(GMAT_MAC_APPBUNDLE_PATH "bin/GMAT-${GMAT_RELEASE_NAME}_Beta.app/Contents")
  endif()
endif()

# Additional Linux compile/link flags
IF(UNIX AND NOT APPLE)
  SET(GMAT_RUN_ADDRESS_SANITIZER OFF CACHE BOOL "Check Memory using Address Sanitizer")

  # Tell CMake Find* commands to search for lib64 paths on Linux
  SET_PROPERTY(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS ON)

  # Generate link errors for undefined symbols 
  # This is enabled by default on OSX but not Linux
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")
  IF (GMAT_RUN_ADDRESS_SANITIZER)
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
    SET(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -lasan -fno-omit-frame-pointer -fsanitize=address")
    SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -lasan -fno-omit-frame-pointer -fsanitize=address")
  ENDIF()
ENDIF()

# ====================================================================
# Allow user to enable GMAT components

INCLUDE(CMakeDependentOption)
OPTION(GMAT_INCLUDE_CSALT "Build CSALT with GMAT" ON)
CMAKE_DEPENDENT_OPTION(GMAT_INCLUDE_CSALT_TESTPROGRAM "Build CSALT test program" ON "GMAT_INCLUDE_CSALT" ON)
CMAKE_DEPENDENT_OPTION(GMAT_INCLUDE_CSALT_INSTALL_SNOPT "Install SNOPT libraries with CSALT" ON "GMAT_INCLUDE_CSALT" ON)

# ====================================================================
# Enable boost::variant as needed

# Use if on OSX before 10.15
if(APPLE AND (CMAKE_OSX_DEPLOYMENT_TARGET VERSION_LESS 10.15))
   SET(GMAT_USE_BOOST_VARIANT ON)
else()
   SET(GMAT_USE_BOOST_VARIANT OFF)
endif()

# Use if C++17 not available
if(NOT ${CMAKE_CXX_STANDARD} MATCHES "17")
   SET(GMAT_USE_BOOST_VARIANT ON)
endif()

# ====================================================================
# Mark Boost as a dependency if needed

if (GMAT_USE_BOOST_VARIANT OR GMAT_INCLUDE_CSALT)
   SET(GMAT_USE_BOOST ON)
   
   # On macOS, don't show warnings related to Boost deprecated declarations
   if(APPLE)
      SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
   endif()
else()
   SET(GMAT_USE_BOOST OFF)
endif()

# Additional libraries needed to use on Mac/Linux
if(UNIX)
  if(APPLE)
    FIND_LIBRARY(COREFOUNDATION_LIBRARY CoreFoundation)
    FIND_LIBRARY(CORESERVICES_LIBRARY CoreServices)
  else()
    FIND_PACKAGE(Threads)
  endif()
endif()

# ====================================================================
# Setup CMake-managed dependencies
ADD_SUBDIRECTORY(depends)

# ====================================================================
# ------------------------------------------------------------------ #
# CSALT Dependencies: Only required if building the CSALT subsystem  #
# ------------------------------------------------------------------ #
if (GMAT_INCLUDE_CSALT)
   SET(CSALT_DEPENDS_FAIL FALSE)

   # Find SNOPT in depends folder
   SET(SNOPT_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/depends/snopt7" CACHE PATH "SNOPT Root Directory")
   SET(SNOPT_LIB_DIR "${SNOPT_ROOT_DIR}/lib" CACHE PATH "SNOPT Library Directory")
   SET(SNOPT_INCLUDE_DIR ${SNOPT_ROOT_DIR}/include)
   FIND_LIBRARY(SNOPT_CPPLIBRARY snopt7_cpp PATHS ${SNOPT_LIB_DIR} NO_DEFAULT_PATH)
   FIND_LIBRARY(SNOPT_LIBRARY snopt7 PATHS ${SNOPT_LIB_DIR} NO_DEFAULT_PATH)
   if(NOT SNOPT_CPPLIBRARY)
      MESSAGE(WARNING "SNOPT CPP library not found. Check that SNOPT_ROOT_DIR contains a valid SNOPT installation.")
      SET(CSALT_DEPENDS_FAIL TRUE)
   endif()
   if(NOT SNOPT_LIBRARY)
      MESSAGE(WARNING "SNOPT library not found. Check that SNOPT_ROOT_DIR contains a valid SNOPT installation.")
      SET(CSALT_DEPENDS_FAIL TRUE)
   endif()

   # Fail Configuration if any dependencies are missing
   if(CSALT_DEPENDS_FAIL)
      MESSAGE(FATAL_ERROR "Missing required CSALT dependencies: aborting build.")
   else()
      MESSAGE(STATUS "Found SNOPT: " ${SNOPT_ROOT_DIR})
   endif()
endif()

# ====================================================================
# Common build options

if(WIN32)
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W2")
endif()

# Enable parallel builds on VisualStudio
if(MSVC)
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
endif()

if(UNIX)
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-strict-aliasing")
endif()

# check fPIC flag on GCC compilers
include(CheckCXXCompilerFlag)
if(CMAKE_COMPILER_IS_GNUCC)
  CHECK_CXX_COMPILER_FLAG("-fPIC" GCC_PIC)
  if(GCC_PIC)
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  endif()
endif()

# flags to append to CMAKE_CXX_FLAGS according to CMAKE_BUILD_TYPE 
if(UNIX)
  SET(CMAKE_CXX_FLAGS_RELEASE "-O3")
  SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
  SET(CMAKE_CXX_FLAGS_MINSIZEREL "-O3")
endif()

# Common definitions
ADD_DEFINITIONS(-DNO_GCC_PRAGMA)
ADD_DEFINITIONS(-DUNICODE -D_UNICODE)

if(GMAT_USE_BOOST_VARIANT)
   ADD_DEFINITIONS(-DGMAT_USE_BOOST_VARIANT)
endif()

# Mark debug libraries with a "d" suffix
SET(CMAKE_DEBUG_POSTFIX "d")

# Warning suppression definitions
# Note that warning C4430 presents as an error in the MatlabInterface plugin
if(WIN32)
  ADD_DEFINITIONS("/wd4430 /wd4251 /wd4231 /wd4996 /nologo")
endif()

# Set 32/64-bit specific definitions and directories
if(GMAT_64_BIT)
  ADD_DEFINITIONS("-DUSE_64_BIT_LONGS")
else()
endif()

# Additional *NIX-specific definitions
if(UNIX)
  ADD_DEFINITIONS("-DLINUX_MAC")

  if(APPLE)
    # GMAT GUI has many instances of text like "XX"GUI_ACCEL_KEY"YY"
    # This is incompatible with the C++11 standard, and should be fixed
    # For now, disable this compile error (-Wreserved-user-defined-literal)
    # [UPDATE] Removed since this no longer seems to create build errors.
    #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-reserved-user-defined-literal")
    
    # Silence OpenGL deprecation warnings
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGL_SILENCE_DEPRECATION")
  endif()
endif()

# Enable RPATH on Mac/Linux so that GMAT shared libraries can be found
# relative to each other. This is necessary to distribute GMAT and use it
# from 3rd-party apps such as Matlab.
if(UNIX)
  if(APPLE)
    SET(CMAKE_MACOSX_RPATH TRUE) # Enable Mac RPATH
  else()
    # Linux searches the RPATH before the system LD_LIBRARY_PATH, so enable
    # the use of RUNPATH which is searched after LD_LIBRARY_PATH. This makes
    # the Linux search order similar to Mac/Windows.
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--enable-new-dtags")
    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--enable-new-dtags")
  endif()
endif()

# ====================================================================
# Define macro that sets the output directory for build products,
# e.g. executables or shared libraries. By default these are placed under
# the <GMAT>/application/ directory so that devs can easily test without
# having to install the build (via make install or VS INSTALL).
# This default location can be changed by the user.
SET(GMAT_BUILDOUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/application" CACHE PATH "Base directory for GMAT binaries. Should contain bin, plugins, etc.")
SET(GMAT_BUILDOUTPUT_DEBUGDIR "${GMAT_BUILDOUTPUT_DIRECTORY}/debug")
MACRO(_SETOUTPUTDIRECTORY TargetName OUTPUT_DIR)
  SET_TARGET_PROPERTIES(${TargetName} PROPERTIES 
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ${GMAT_BUILDOUTPUT_DEBUGDIR}/${OUTPUT_DIR}
    LIBRARY_OUTPUT_DIRECTORY_DEBUG ${GMAT_BUILDOUTPUT_DEBUGDIR}/${OUTPUT_DIR}
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${GMAT_BUILDOUTPUT_DIRECTORY}/${OUTPUT_DIR}
    LIBRARY_OUTPUT_DIRECTORY_RELEASE ${GMAT_BUILDOUTPUT_DIRECTORY}/${OUTPUT_DIR}
    )
ENDMACRO()

# ====================================================================
# Define macro that adds source/header files to appropriate groups
# for IDE-based projects (e.g. VisualStudio, XCode, ...).
MACRO(_ADDSOURCEGROUPS ALL_DIRS)
  # Loop over each directory containing source/header files
  FOREACH(currdir ${ALL_DIRS})
    # Create the source group name from the current directory name
    STRING(REPLACE "/" "\\\\" groupname ${currdir})

    # Add source files (.c, .cpp) to IDE group
    SOURCE_GROUP("Source Files\\${groupname}" REGULAR_EXPRESSION "${currdir}/.*\\.c(pp)?")

    # Add header files (.hpp) to IDE group
    SOURCE_GROUP("Header Files\\${groupname}" REGULAR_EXPRESSION "${currdir}/.*\\.h(pp)?")
  ENDFOREACH(currdir)
ENDMACRO(_ADDSOURCEGROUPS)

# ====================================================================
# Setup GMAT directory structure
INCLUDE(build/install/CMake_INSTALL_Setup.cmake)

# ====================================================================
# Go to src directory and look for CMake instructions there
ADD_SUBDIRECTORY(src)

# ====================================================================
# Setup GMAT install process
ADD_SUBDIRECTORY(build/install)
