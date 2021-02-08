cmake_minimum_required(VERSION 3.5)
project(PathPatternMining)

#SET(BOOST_ROOT ~/boost_1_63_0)
#find_package(Boost 1.63.0 REQUIRED COMPONENTS serialization)
#find_library(boost_serialization BOOST_ROOT)
find_package(Boost REQUIRED COMPONENTS serialization)


set(CMAKE_CXX_STANDARD 11)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()


#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")


include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS} )

#output dir
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

set(SRC_List GraphManager.cpp RandSim.cpp RandSim.h
        ExceptionalFactMining.cpp ExceptionalFactMining.h
        QueryInterface.cpp QueryInterface.h
        FMiner.cpp FMiner.h settings.h)

add_executable(PathPatternMining main.cpp ${SRC_List})
target_link_libraries(PathPatternMining ${Boost_LIBRARIES})
