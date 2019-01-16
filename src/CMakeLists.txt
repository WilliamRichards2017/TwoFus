INCLUDE_DIRECTORIES(
    ${BAMTOOLS_INCLUDE}
    ${ZLIB_INCLUDE}
    ${MINIMAP2_INCLUDE}
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")




add_library(twofus_core STATIC
	clipCoords.cpp
	util.cpp
	input.cpp	
	MEHead.cpp
	polyTail.cpp
	mobileElement.cpp
	insertion.cpp
	contigs.cpp
	vcfWriter.cpp
)

SET(CORE_LIB twofus_core CACHE STRING "")

TARGET_LINK_LIBRARIES(${CORE_LIB}
	${BAMTOOLS_LIB}
	${ZLIB_LIBRARY}
	${MINIMAP2_LIB}
)

add_dependencies(${CORE_LIB} ${TWOFUS_EXTERNAL_PROJECT})
