INCLUDE_DIRECTORIES( ${Farsight_SOURCE_DIR}/InsightJournalFilters/MinErrorThresholding )
ADD_LIBRARY( ftkPreprocess
	ftkPreprocess2.h
	ftkPreprocess2.cpp
	itkBinaryThinningImageFilter3D.h
)
TARGET_LINK_LIBRARIES( ftkPreprocess  TinyXML     ${ITK_LIBRARIES} )

ADD_EXECUTABLE( prep preprocess_main.cpp )
TARGET_LINK_LIBRARIES( prep ftkPreprocess )

INSTALL(TARGETS prep RUNTIME DESTINATION ${INSTALL_BIN_DIR})
