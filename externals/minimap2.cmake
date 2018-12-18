# =================================
# TwoFus/externals/minimap2.cmake
# =================================

SET_PROPERTY(DIRECTORY PROPERTY "EP_BASE" ${ep_base})

SET(MINIMAP2_PROJECT minimap2_project CACHE INTERNAL "minimap2 project name")
SET(MINIMAP2_DIR ${CMAKE_BINARY_DIR}/externals/minimap2 CACHE INTERNAL "minimap2 project directory")

ExternalProject_Add(${MINIMAP2_PROJECT}
	GIT_REPOSITORY https://github.com/lh3/minimap2.git
	GIT_TAG master
        CONFIGURE_COMMAND ""
	BUILD_COMMAND "make"
	INSTALL_COMMAND ""
	UPDATE_COMMAND ""
	BUILD_IN_SOURCE 1
    PREFIX ${MINIMAP2_DIR}
    CMAKE_CACHE_ARGS
        -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}
)

ExternalProject_Get_Property(${MINIMAP2_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${MINIMAP2_PROJECT} BINARY_DIR)

SET(MINIMAP2_LIB ${BINARY_DIR}/libminimap2.a CACHE INTERNAL "MINIMAP2 Library")
SET(MINIMAP2_INCLUDE ${SOURCE_DIR} CACHE INTERNAL "MINIMAP2 Include")