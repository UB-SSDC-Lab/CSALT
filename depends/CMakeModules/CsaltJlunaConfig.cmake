
IF(NOT GMAT_INCLUDE_CSALT_WITH_JULIA_PROGRAM)
  RETURN()
ENDIF()

# Find path to JLUNA root
FIND_PATH(JLUNA_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/jluna)

IF ((NOT JLUNA_ROOT) OR (NOT EXISTS ${JLUNA_ROOT}))
  # Message that we need to build jluna
  MESSAGE("-- Downloading jluna")

  # We have a submodule setup for jluna
ELSE()

ENDIF()