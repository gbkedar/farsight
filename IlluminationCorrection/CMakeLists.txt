SET( GenerateTilesNPairList_SRCS mainTilesNPairs.cpp )
SET( Illumination_SRCS IlluminationFull.cpp )
SET( CziToNrrd_SRCS readCziWriteCorrectedTiles.cpp )
SET( NrrdToDicom_SRCS NrrdToDicom.cpp )

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
TARGET_LINK_LIBRARIES( Illumination ${ITK_LIBRARIES} ${MLPack_LIBRARY} ${LIBXML2_LIBRARIES} )

ADD_EXECUTABLE( NrrdToDicom ${NrrdToDicom_SRCS} )
TARGET_LINK_LIBRARIES( NrrdToDicom ${ITK_LIBRARIES} TinyXML ${Boost_LIBRARIES} )

INSTALL( TARGETS GenerateTilesNPairList Illumination CziToNrrd RUNTIME DESTINATION ${INSTALL_BIN_DIR} )
