cmake_minimum_required (VERSION 3.1)
project (shannon_cpp)
#set(CMAKE_CXX_COMPILER /usr/bin/gcc-4.9*)
set (CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)

include_directories(/home/bowen/boost_1_72_0/include/)
include_directories(./sparsepp/)

SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
SET(BUILD_SHARED_LIBS OFF)
SET(CMAKE_EXE_LINKER_FLAGS "-static")

#SET(BOOST_ROOT /usr)
#SET(CMAKE_INCLUDE_PATH /usr/include ${CMAKE_INCLUDE_PATH} )
#SET(CMAKE_LIBRARY_PATH /usr/lib ${CMAKE_LIBRARY_PATH} )
#set(BOOST_INCLUDE_DIRS /usr/include/boost)

#set(BOOST_DEBUG)
#find_package(Boost 1.65.1 COMPONENTS program_options system filesystem REQUIRED)
find_library(METIS_LIB_ABS_PATH metis)
find_package(ZLIB)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib_shannon)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib_shannon)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#message( ${Boost_FILESYSTEM_LIBRARY} )
#message( ${Boost_PROGRAM_OPTIONS_LIBRARY} )
#message( ${Boost_FILESYSTEM_LIBRARY} )

SET(boost_filesystem_local /home/bowen/boost_1_72_0/lib/libboost_filesystem.a)
SET(boost_system_local /home/bowen/boost_1_72_0/lib/libboost_system.a)
SET(boost_program_options_local /home/bowen/boost_1_72_0/lib/libboost_program_options.a)
SET(glpk_local /usr/local/lib/libglpk.a)



add_library( primary_lib
	src/log.cpp
	src/basic_log.cpp
	src/local_file_structure.cpp
	src/json_parser.cpp
	src/encoding.cpp
    src/shannon_C_seq_helper.cpp
)
target_link_libraries(primary_lib ${boost_filesystem_local})
target_link_libraries(primary_lib ${boost_program_options_local})
target_link_libraries(primary_lib ${ZLIB_LIBRARIES})

add_library( contig_graph_handler	
	src/Contig_graph_handler.cpp
	src/Kmer_handler.cpp
	src/Contig_handler.cpp	
)
target_link_libraries(contig_graph_handler ${boost_system_local})
target_link_libraries(contig_graph_handler ${METIS_LIB_ABS_PATH})
#set_target_properties(contig_graph_handler PROPERTIES ARCHIVE_OUTPUT_DIRECTORY  ./lib_shannon/ )

add_library( seq_graph_handler 
	src/Sequence_graph_handler.cpp	
)
target_link_libraries(seq_graph_handler pthread)
#set_target_properties(seq_graph_handler PROPERTIES ARCHIVE_OUTPUT_DIRECTORY  ./lib_shannon/ )

add_library( sparse_flow_handler 
	src/Sparse_flow_handler.cpp	
)
target_link_libraries(sparse_flow_handler ${glpk_local})
target_link_libraries(sparse_flow_handler pthread)
target_link_libraries(sparse_flow_handler ${boost_system_local})

#set_target_properties(sparse_flow_handler PROPERTIES ARCHIVE_OUTPUT_DIRECTORY  ./lib_shannon/ )

add_library( multi_graph_handler
	src/Multi_graph_handler.cpp
)
target_link_libraries(multi_graph_handler sparse_flow_handler)
target_link_libraries(multi_graph_handler seq_graph_handler)
#set_target_properties(multi_graph_handler PROPERTIES ARCHIVE_OUTPUT_DIRECTORY  ./lib_shannon/ )

set(SHANNON_DEVELOP_SOURCE 
    src/main.cpp
    src/run_tasks.cpp
)
# add the executable
add_executable(shannon_cpp ${SHANNON_DEVELOP_SOURCE})
target_link_libraries(shannon_cpp multi_graph_handler )
target_link_libraries(shannon_cpp contig_graph_handler )
target_link_libraries(shannon_cpp primary_lib)
