CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

set(CUDA_USE_STATIC_CUDA_RUNTIME OFF)

SET(CMAKE_BUILD_TYPE Release)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -Wall  -O3 -march=native ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall   -O3 -march=native")

# Check C++11 or C++0x support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
   add_definitions(-DCOMPILEDWITHC11)
   message(STATUS "Using flag -std=c++11.")
elseif(COMPILER_SUPPORTS_CXX0X)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
   add_definitions(-DCOMPILEDWITHC0X)
   message(STATUS "Using flag -std=c++0x.")
else()
   message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
   message("OPENMP FOUND")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
   message("OPENMP FOUND")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


#set(SuiteSparse_DIR "/home/jpl/SuiteSparse/")
#set(SuiteSparse_INCLUDE_DIRS "/home/jpl/SuiteSparse/include")
#set(SuiteSparse_LIB_DIRS "/home/jpl/SuiteSparse/lib")

LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

find_package(OpenCV REQUIRED)
#find_package(VTK REQUIRED)
find_package(Eigen3 REQUIRED)
#find_package(Boost REQUIRED)
find_package(dlib REQUIRED)


find_package(SuiteSparse QUIET NO_MODULE)  # 1st: Try to locate the *config.cmake file.


if(NOT SuiteSparse_FOUND)
        #set(SuiteSparse_VERBOSE ON)
        find_package(SuiteSparse REQUIRED) # 2nd: Use FindSuiteSparse.cmake module
        include_directories(${SuiteSparse_INCLUDE_DIRS})
else()
        message(STATUS "Find SuiteSparse : include(${USE_SuiteSparse})")
        include(${USE_SuiteSparse})
endif()
MESSAGE(STATUS "SuiteSparse_LIBS: ${SuiteSparse_LIBRARIES}")


#INCLUDE (${VTK_USE_FILE})

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${OpenCV_INCLUDE_DIRS})
include_directories(${dlib_INCLUDE_DIRS})

set(Trimesh2_DIR "/home/jpl/Programs/trimesh2")
set(Trimesh2_INCLUDE_DIRS "/home/jpl/Programs/trimesh2/include/")
set(Trimesh2_LIB_DIRS "/home/jpl/Programs/trimesh2/lib.Linux64/libtrimesh.a")
include_directories(${Trimesh2_INCLUDE_DIRS})


add_executable(nonrigid_icp main.cpp)
target_link_libraries(nonrigid_icp  ${OpenCV_LIBS} ${SuiteSparse_LIBRARIES} ${dlib_LIBRARIES}  ${Trimesh2_LIB_DIRS})


#add_executable(LaplacianDeformer LaplacianDeformer.cpp)
#target_link_libraries(LaplacianDeformer ${OpenCV_LIBS} ${SuiteSparse_LIBRARIES} ${Trimesh2_LIB_DIRS})
