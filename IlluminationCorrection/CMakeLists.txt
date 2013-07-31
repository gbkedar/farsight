SET( GenerateTilesNPairList_SRCS mainTilesNPairs.cpp )
SET( GenetateFG_SRCS mainFG.cpp )

INCLUDE_DIRECTORIES( ${Farsight_SOURCE_DIR}/NuclearSegmentation/yousef_core/cell_binarization )
INCLUDE_DIRECTORIES( ${Farsight_SOURCE_DIR}/InsightJournalFilters/MinErrorThresholding )

#ADD_EXECUTABLE( GenerateTilesNPairList ${GenerateTilesNPairList_SRCS} )
#TARGET_LINK_LIBRARIES( GenerateTilesNPairList ${ITK_LIBRARIES} TinyXML ${Boost_LIBRARIES} )

ADD_EXECUTABLE( GenerateFG ${GenetateFG_SRCS} )
TARGET_LINK_LIBRARIES( GenerateFG ${ITK_LIBRARIES} )

#INSTALL( TARGETS GenerateTilesNPairList GenerateFG RUNTIME DESTINATION ${INSTALL_BIN_DIR} )
