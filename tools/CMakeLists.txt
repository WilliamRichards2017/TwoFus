

ADD_DEFINITIONS(-DBOOST_FALLTHROUGH)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

INCLUDE_DIRECTORIES(
	${ZLIB_INCLUDE}
	${BAMTOOLS_INCLUDE}
	${PROJECT_SOURCE_DIR}

)

add_executable(runTwofus
	runTwofus.cpp
)

target_link_libraries(runTwofus
	${CORE_LIB}
)
