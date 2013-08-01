# - Find MLPack (includes and library)
# This module defines
#  MLPack_INCLUDE_DIR
#  MLPack_LIBRARIES
#  MLPack_FOUND
# also defined, but not for general use are
#  MLPack_LIBRARY, where to find the library.

FIND_PATH(MLPack_INCLUDE_DIR mlpack/core.hpp
HINTS /usr/include/ /usr/local/include/ /usr/include/ /usr/local/include/ )

SET(MLPack_NAMES ${MLPack_NAMES} mlpack)
FIND_LIBRARY(MLPack_LIBRARY
  NAMES ${MLPack_NAMES}
  PATHS /usr/lib64 /usr/local/lib64 /usr/lib /usr/local/lib
  )

IF (MLPack_LIBRARY AND MLPack_INCLUDE_DIR)
    SET(MLPack_LIBRARIES ${MLPack_LIBRARY})
    SET(MLPack_FOUND "YES")
ELSE (MLPack_LIBRARY AND MLPack_INCLUDE_DIR)
  SET(MLPack_FOUND "NO")
ENDIF (MLPack_LIBRARY AND MLPack_INCLUDE_DIR)

IF (MLPack_FOUND)
   IF (NOT MLPack_FIND_QUIETLY)
      MESSAGE(STATUS "Found a MLPack library: ${MLPack_LIBRARIES}")
   ENDIF (NOT MLPack_FIND_QUIETLY)
ELSE (MLPack_FOUND)
   IF (MLPack_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find a MLPack library")
   ENDIF (MLPack_FIND_REQUIRED)
ENDIF (MLPack_FOUND)

MARK_AS_ADVANCED(
  MLPack_LIBRARY
  MLPack_INCLUDE_DIR
  )
