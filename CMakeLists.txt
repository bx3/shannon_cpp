cmake_minimum_required (VERSION 3.1)
project (Shannon_RNASeq_Cpp)
#set(CMAKE_CXX_COMPILER /usr/bin/gcc-4.9*)
set (CMAKE_CXX_STANDARD 11)

include_directories(/home/bowen/boost_1_65_1)

SET(CMAKE_INCLUDE_PATH /usr/include ${CMAKE_INCLUDE_PATH} )
SET(CMAKE_LIBRARY_PATH /usr/lib ${CMAKE_LIBRARY_PATH} )
#set(BOOST_ROOT "/usr")
#set(BOOST_INCLUDEDIR "/usr/include")
#set(BOOST_LIBRARYDIR "/usr/lib")
set(BOOST_DEBUG)
find_package(Boost 1.65.1 COMPONENTS system filesystem REQUIRED)

find_library(METIS_LIB_ABS_PATH metis)


set(SHANNON_SOURCE
	src/encoding.cpp
	src/Kmer_handler.cpp
	src/basic_log.cpp
	src/log.cpp
	src/main.cpp
	src/basic_log.cpp
	src/Contig_handler.cpp
	src/Contig_graph_handler.cpp
	src/local_file_structure.cpp
	src/json_parser.cpp
	src/Sequence_graph_handler.cpp
    src/Multibridge_handler.cpp
)
 
# add the executable
add_executable(Shannon_RNASeq_Cpp ${SHANNON_SOURCE})

target_link_libraries(Shannon_RNASeq_Cpp ${Boost_SYSTEM_LIBRARY})
target_link_libraries(Shannon_RNASeq_Cpp ${Boost_FILESYSTEM_LIBRARY})
target_link_libraries(Shannon_RNASeq_Cpp ${METIS_LIB_ABS_PATH})
target_link_libraries(Shannon_RNASeq_Cpp pthread)
