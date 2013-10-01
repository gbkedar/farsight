IF( FARSIGHT_VTK_6 )
	SET( VTKLSM_H )
	SET( VTKLSM_CXX )
ELSE( FARSIGHT_VTK_6 )
	SET( VTKLSM_H vtkLSMReader.h )
	SET( VTKLSM_CXX vtkLSMReader.cxx )
ENDIF( FARSIGHT_VTK_6 )

SET( FTKIMAGE_SRCS    
	${VTKLSM_CXX}
	ftkImage.cpp
)

SET( FTKIMAGE_HDRS
	${VTKLSM_H}
	vtkBXDProcessingWin32Header.h
	ftkImage.h
	ftkImage.txx
)
IF(MSVC)
  SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj" )
ENDIF(MSVC)

ADD_LIBRARY( ftkImage
             ${FTKIMAGE_SRCS} 
             ${FTKIMAGE_HDRS} 
)

TARGET_LINK_LIBRARIES( ftkImage 
  ${VTK_LIBRARIES}
  ${ITK_LIBRARIES} )

#INSTALL( TARGETS ftkImage DESTINATION . )
                           
