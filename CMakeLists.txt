cmake_minimum_required (VERSION 3.1)
project (shannon_cpp)
set (CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)


include_directories(./sparsepp/)
include_directories(/usr/local/include/)


find_package(Boost 1.65.0 COMPONENTS program_options system filesystem REQUIRED)
find_library(METIS_LIB_ABS_PATH metis)
find_library(local_glpk glpk)


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib_shannon)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib_shannon)



find_package(ZLIB)


add_library( primary_lib
	src/log.cpp
	src/basic_log.cpp
	src/local_file_structure.cpp
	src/json_parser.cpp
	src/encoding.cpp
    src/shannon_C_seq_helper.cpp
)
target_link_libraries(primary_lib ${Boost_FILESYSTEM_LIBRARY})
target_link_libraries(primary_lib ${Boost_PROGRAM_OPTIONS_LIBRARY})
target_link_libraries(primary_lib ${ZLIB_LIBRARIES})

add_library( contig_graph_handler	
	src/Contig_graph_handler.cpp
	src/Kmer_handler.cpp
	src/Contig_handler.cpp	
)
target_link_libraries(contig_graph_handler ${Boost_SYSTEM_LIBRARY})
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
target_link_libraries(sparse_flow_handler ${local_glpk})
target_link_libraries(sparse_flow_handler pthread)
target_link_libraries(sparse_flow_handler ${Boost_SYSTEM_LIBRARY})

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
