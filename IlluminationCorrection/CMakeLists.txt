SET( GenerateTilesNPairList_SRCS mainTilesNPairs.cpp )
SET( Illumination_SRCS illuminationFull.cpp )
SET( CziToNrrd_SRCS readCziWriteCorrectedTiles.cpp )

INCLUDE_DIRECTORIES( ${Farsight_SOURCE_DIR}/NuclearSegmentation/yousef_core/cell_binarization )

ADD_EXECUTABLE( GenerateTilesNPairList ${GenerateTilesNPairList_SRCS} )
TARGET_LINK_LIBRARIES( GenerateTilesNPairList ${ITK_LIBRARIES} TinyXML ${Boost_LIBRARIES} )

ADD_EXECUTABLE( CziToNrrd ${CziToNrrd_SRCS} )
TARGET_LINK_LIBRARIES( CziToNrrd ${ITK_LIBRARIES} ${Boost_LIBRARIES} )

FIND_PACKAGE( MLPack  REQUIRED )
FIND_PACKAGE( LibXml2 REQUIRED )
INCLUDE_DIRECTORIES( ${MLPack_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${LIBXML2_INCLUDE_DIR} )

ADD_EXECUTABLE( Illumination ${Illumination_SRCS} )
TARGET_LINK_LIBRARIES( Illumination ${ITK_LIBRARIES} ${MLPack_LIBRARY}
				${LIBXML2_LIBRARIES} )

INSTALL( TARGETS GenerateTilesNPairList Illumination RUNTIME DESTINATION ${INSTALL_BIN_DIR} )
