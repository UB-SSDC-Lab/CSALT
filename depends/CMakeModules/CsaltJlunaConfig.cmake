
IF(NOT GMAT_INCLUDE_CSALT_WITH_JULIA_PROGRAM)
  RETURN()
ENDIF()

# Set jluna paths
SET(jluna_ROOT        ${CMAKE_CURRENT_SOURCE_DIR}/jluna)
SET(jluna_BUILD_DIR   ${GMAT_DEP_BUILD_DIR}/jluna)
SET(jluna_INSTALL_DIR ${jluna_ROOT}/lib)

# If jluna root does not exist, initialize git submodule
IF (NOT EXISTS "${jluna_ROOT}/.src")
  MESSAGE(STATUS "Downloading jluna ...")
  EXECUTE_PROCESS(COMMAND git submodule update --init -- ${CMAKE_CURRENT_SOURCE_DIR}/jluna
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

# Make build directory
FILE(MAKE_DIRECTORY ${jluna_BUILD_DIR})

# If jluna install directory dies not exist, build jluna 
IF (NOT EXISTS ${jluna_INSTALL_DIR})
  MESSAGE(STATUS "Building jluna ...")
  # Configure and generate jluna makefile
  EXECUTE_PROCESS(COMMAND cmake ${jluna_ROOT} 
                                -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                -DCMAKE_INSTALL_PREFIX=${jluna_INSTALL_DIR}
                          WORKING_DIRECTORY ${jluna_BUILD_DIR}
                          OUTPUT_QUIET)
  EXECUTE_PROCESS(COMMAND make install 
                          WORKING_DIRECTORY ${jluna_BUILD_DIR}
                          OUTPUT_QUIET)
ENDIF()

# Find jluna package
SET(jluna_INCLUDE_DIR "${jluna_INSTALL_DIR}/include/jluna" PARENT_SCOPE)