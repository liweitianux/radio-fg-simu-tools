#
# CMakeLists.txt for 'fg_simulation'
#
# Aaron LI
# 2015/04/16
#

cmake_minimum_required( VERSION 2.8 )

project( sf_simulation C CXX )

include( CheckCXXCompilerFlag )
CHECK_CXX_COMPILER_FLAG( "-std=c++11" COMPILER_SUPPORTS_CXX11 )
CHECK_CXX_COMPILER_FLAG( "-std=c++0x" COMPILER_SUPPORTS_CXX0X )
if( COMPILER_SUPPORTS_CXX11 )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )
elseif( COMPILER_SUPPORTS_CXX0X )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x" )
else()
    message( STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support." )
endif()

if ( UNIX )
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall" )
endif ( UNIX )

option( USE_GSL_CBLAS "Use gsl cblas library" ON )

set( GSL_LIBS gsl m )
if ( USE_GSL_CBLAS )
    set( GSL_LIBS ${GSL_LIBS} gslcblas )
endif ( USE_GSL_CBLAS )

add_executable( calc_NN0 calc_NN0.c )
target_link_libraries( calc_NN0 ${GSL_LIBS} )
add_executable( calc_Fnu_t calc_Fnu_t.c )
target_link_libraries( calc_Fnu_t ${GSL_LIBS} )

#set( CMAKE_INSTALL_PREFIX ../install )
install( TARGETS calc_NN0 DESTINATION bin )
install( TARGETS calc_Fnu_t DESTINATION bin )

# vim: set ts=8 sw=4 tw=0 fenc=utf-8 ft=cmake: #
