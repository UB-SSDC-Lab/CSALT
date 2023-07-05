
MESSAGE("==============================")
MESSAGE("CSALT application with Julia" ${VERSION}) 

SET(TargetName CSALTApplicationWithJulia)

# Find julia
FIND_PACKAGE(Julia 1.7.0 REQUIRED)

# Find jluna
FIND_PACKAGE(jluna 1.0 REQUIRED)

# Add preprocessor definitions
ADD_DEFINITIONS("-D_DYNAMICLINK")

SET(CSALT_APP_JL_DIRS 
    src
    src/util
    src/drivers
    src/pointpath
    ${CSALT_DIRS}
    ${UTIL_DIRS}
)

SET(CSALT_APP_JL_SRCS
    src/csaltApplicationWithJulia.cpp
    src/util/julia_utils.cpp
    src/ConsoleMessageReceiver.cpp
    src/drivers/CsaltDriver.cpp
    src/drivers/DebrisDeorbitDriver.cpp
    src/drivers/AveragedOrbitalElementsDriver.cpp
    src/pointpath/DebrisDeorbitPathObject.cpp 
    src/pointpath/DebrisDeorbitPointObject.cpp 
    src/pointpath/DebrisDeorbitInitialStateConstraint.cpp 
    src/pointpath/DebrisDeorbitFinalStateConstraint.cpp
    src/pointpath/AveragedOrbitalElementsPathObject.cpp 
    src/pointpath/AveragedOrbitalElementsPointObject.cpp
)

# ====================================================================
# Add preprocessor macro for providing Julia code for including 
# Julia source code
SET(JLSRC_INCLUDE_TEMP "${CMAKE_CURRENT_SOURCE_DIR}/jlsrc/include.jl" CACHE FILEPATH "File path to Julia source code include script.")
FILE(TO_CMAKE_PATH "${JLSRC_INCLUDE_TEMP}" JLSRC_INCLUDE)
ADD_DEFINITIONS(-DJLSRC_INCLUDE="${JLSRC_INCLUDE}")

# ====================================================================
# Recursively find all include files, which will be added to IDE-based
# projects (VS, XCode, etc.)
FILE(GLOB_RECURSE CSALT_APP_JL_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)

# ====================================================================

# add the install targets
ADD_EXECUTABLE(${TargetName} ${CSALT_APP_JL_SRCS} ${CSALT_APP_JL_HEADERS})

# Add necessary include directories
TARGET_INCLUDE_DIRECTORIES(${TargetName} PUBLIC ${CSALT_DIRS})
TARGET_INCLUDE_DIRECTORIES(${TargetName} PUBLIC ${CSALT_APP_JL_DIRS})
TARGET_INCLUDE_DIRECTORIES(${TargetName} PUBLIC ${SNOPT_INCLUDE_DIR})
TARGET_INCLUDE_DIRECTORIES(${TargetName} PUBLIC ${Boost_INCLUDE_DIR})
TARGET_INCLUDE_DIRECTORIES(${TargetName} PUBLIC ${jluna_INCLUDE_DIR})

# Link needed libraries
TARGET_LINK_LIBRARIES(${TargetName} PRIVATE ${SNOPT_CPPLIBRARY})
TARGET_LINK_LIBRARIES(${TargetName} PRIVATE ${SNOPT_LIBRARY})
TARGET_LINK_LIBRARIES(${TargetName} PRIVATE CSALT)
TARGET_LINK_LIBRARIES(${TargetName} PRIVATE GmatUtil)
TARGET_LINK_LIBRARIES(${TargetName} PRIVATE ${jluna})
TARGET_LINK_LIBRARIES(${TargetName} PRIVATE "$<BUILD_INTERFACE:Julia::Julia>")
TARGET_LINK_LIBRARIES(${TargetName} PRIVATE "$<BUILD_INTERFACE:Threads::Threads>")

# Set building with C++ STD 20
TARGET_COMPILE_FEATURES(${TargetName} PRIVATE cxx_std_20)

if(UNIX)
  if(APPLE)
    TARGET_LINK_LIBRARIES(${TargetName} PRIVATE ${COREFOUNDATION_LIBRARY} ${CORESERVICES_LIBRARY})
  else()
    TARGET_LINK_LIBRARIES(${TargetName} PRIVATE Threads::Threads)
  endif()
endif()

# Windows-specific link flags
if(WIN32)
  SET_TARGET_PROPERTIES(${TargetName} PROPERTIES LINK_FLAGS "/NODEFAULTLIB:\"libcmt.lib\" /INCREMENTAL:NO")
endif()

# Create build products in bin directory
_SETOUTPUTDIRECTORY(${TargetName} bin)

# Override debug output directory
SET_TARGET_PROPERTIES(${TargetName} PROPERTIES 
  RUNTIME_OUTPUT_DIRECTORY_DEBUG ${GMAT_BUILDOUTPUT_DEBUGDIR}
  LIBRARY_OUTPUT_DIRECTORY_DEBUG ${GMAT_BUILDOUTPUT_DEBUGDIR}
  )

# Specify where to install (make install or VS "INSTALL" project)
if(WIN32)
  INSTALL(TARGETS ${TargetName} RUNTIME DESTINATION bin)
else()
  if(APPLE)
    SET(INSTALL_DIR ${GMAT_MAC_APPBUNDLE_PATH}/Frameworks)
  else()
    SET(INSTALL_DIR bin)
  endif()
  INSTALL(TARGETS ${TargetName} DESTINATION bin)
endif()

# Set RPATH to find shared libraries in default locations on Mac/Linux
if(UNIX)
  if(APPLE)
    SET_TARGET_PROPERTIES(${TargetName} PROPERTIES INSTALL_RPATH
      "@loader_path/;@loader_path/../plugins/;@loader_path/../plugins/proprietary/"
      )
  else()
    SET_TARGET_PROPERTIES(${TargetName} PROPERTIES INSTALL_RPATH
      "\$ORIGIN/;\$ORIGIN/../plugins/;\$ORIGIN/../plugins/proprietary/"
      )
  endif()
endif()
