# =================================
# TwoFus/externals/gtest.cmake
# =================================

SET_PROPERTY(DIRECTORY PROPERTY "EP_BASE" ${ep_base})

SET(GTEST_PROJECT gtest_project CACHE INTERNAL "gtest project name")
SET(GTEST_DIR ${CMAKE_BINARY_DIR}/externals/gtest CACHE INTERNAL "gtest project directory")
SET(GTEST_LIB)
ExternalProject_Add(${GTEST_PROJECT}
	GIT_REPOSITORY https://github.com/google/googletest.git
	GIT_TAG 9ae086a9ebafabdc49b71bb7f3879f551adee09a #lock in the commit id so we don't this doesn't break in the future
	INSTALL_COMMAND ""
	UPDATE_COMMAND ""
	PREFIX ${GTEST_DIR}
    CMAKE_CACHE_ARGS
        -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}
)

ExternalProject_Get_Property(${GTEST_PROJECT} INSTALL_DIR)
ExternalProject_Get_Property(${GTEST_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${GTEST_PROJECT} BINARY_DIR)

SET(GTEST_LIB ${BINARY_DIR}/googlemock/gtest/libgtest.a CACHE INTERNAL "GTEST Lib")
SET(GTEST_SOURCE_DIR ${SOURCE_DIR}/googletest/include/ CACHE INTERNAL "GTEST Include")