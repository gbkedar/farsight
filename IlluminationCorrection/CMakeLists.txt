SET( GenerateTilesNPairList_SRCS mainTilesNPairs.cpp )
SET( GenetateFG_SRCS mainFG.cpp )

ADD_EXECUTABLE( GenerateTilesNPairList ${GenerateTilesNPairList_SRCS} )
TARGET_LINK_LIBRARIES( GenerateTilesNPairList ${ITK_LIBRARIES} TinyXML ${Boost_LIBRARIES} )

ADD_EXECUTABLE( GenerateFG ${GenetateFG_SRCS} )
TARGET_LINK_LIBRARIES( GenerateFG ${ITK_LIBRARIES} )

INSTALL( TARGETS GenerateTilesNPairList GenerateFG RUNTIME DESTINATION ${INSTALL_BIN_DIR} )
