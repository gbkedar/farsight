SET( GenerateTilesNPairList_SRCS mainTilesNPairs.cpp )
SET( Illumination_SRCS IlluminationFull.cpp )
SET( CziToNrrd_SRCS readCziWriteCorrectedTiles.cpp )
SET( NrrdToDicom_SRCS NrrdToDicom.cpp )
SET( Resample_SRCS ResampleForRegistration.cpp )

INCLUDE_DIRECTORIES( ${Farsight_SOURCE_DIR}/NuclearSegmentation/yousef_core/cell_binarization )

ADD_EXECUTABLE( GenerateTilesNPairList ${GenerateTilesNPairList_SRCS} )
TARGET_LINK_LIBRARIES( GenerateTilesNPairList ${ITK_LIBRARIES} TinyXML ${Boost_LIBRARIES} )

ADD_EXECUTABLE( CziToNrrd ${CziToNrrd_SRCS} )
TARGET_LINK_LIBRARIES( CziToNrrd ${ITK_LIBRARIES} ${Boost_LIBRARIES} )
FIND_PACKAGE(Armadillo 2.4.2 REQUIRED)
FIND_PACKAGE( MLPack  REQUIRED )
FIND_PACKAGE( LibXml2 REQUIRED )
INCLUDE_DIRECTORIES( ${ARMADILLO_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${MLPack_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${LIBXML2_INCLUDE_DIR} )

message("boost libs ${Boost_LIBRARIES}" )
message("and arma libs ${ARMADILLO_LIBRARIES} with mlpack ${MLPack_LIBRARY}" )

ADD_EXECUTABLE( Illumination ${Illumination_SRCS} )
TARGET_LINK_LIBRARIES( Illumination ${MLPack_LIBRARY} ${ARMADILLO_LIBRARIES} ${Boost_LIBRARIES} ${LIBXML2_LIBRARIES} ${ITK_LIBRARIES} )

ADD_EXECUTABLE( ResampleImage ${Resample_SRCS} )
TARGET_LINK_LIBRARIES( ResampleImage ${ITK_LIBRARIES} )

ADD_EXECUTABLE( NrrdToDicom ${NrrdToDicom_SRCS} )
TARGET_LINK_LIBRARIES( NrrdToDicom ${ITK_LIBRARIES} TinyXML ${Boost_LIBRARIES} )

#INSTALL( TARGETS GenerateTilesNPairList Illumination CziToNrrd RUNTIME DESTINATION ${INSTALL_BIN_DIR} )
