SET( GenerateTilesNPairList_SRCS mainTilesNPairs.cpp )

ADD_EXECUTABLE( GenerateTilesNPairList ${GenerateTilesNPairList_SRCS} )
TARGET_LINK_LIBRARIES( GenerateTilesNPairList ${ITK_LIBRARIES} TinyXML ${Boost_LIBRARIES} )
INSTALL( TARGETS GenerateTilesNPairList RUNTIME DESTINATION ${INSTALL_BIN_DIR} )
