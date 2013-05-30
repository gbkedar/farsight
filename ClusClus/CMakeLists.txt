SET( Libraries
    ${VTK_LIBRARIES}
    vcl
    vnl
    ${QT_LIBRARIES})

SET( ClusClus_SRCS 
	Kmeans.cpp
	Tree.cpp
	Biclustering.cpp 
	Dendrogram.cpp
	clusclus.cpp
	clusgap.cpp
	clusclusMainwindow.cpp
	HeatmapWindow.cpp
	Heatmap.cpp
	LocalGeometryRef.cpp)

SET( ClusClus_HDRS 
	Kmeans.h
	Tree.h
	Biclustering.h
	Dendrogram.h
	clusclus.h
	clusgap.h
	clusclusMainwindow.h
	HeatmapWindow.h
	Heatmap.h
	LocalGeometryRef.h)
	
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS})
QT4_WRAP_CPP( ClusClus_MOC_REQD_HDRS HeatmapWindow.h Dendrogram.h clusclusMainwindow.h Heatmap.h
		OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED -DBOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION ) 

ADD_LIBRARY( ClusClus_lib ${ClusClus_SRCS} ${ClusClus_HDRS} ${ClusClus_MOC_REQD_HDRS} )

TARGET_LINK_LIBRARIES( ClusClus_lib  
    ${Libraries}
    ftkGUI
    ftkCommon)           

ADD_EXECUTABLE( ClusClus main.cpp )
  
TARGET_LINK_LIBRARIES( ClusClus ClusClus_lib)

