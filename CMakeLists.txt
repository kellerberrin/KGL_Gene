cmake_minimum_required(VERSION 3.8)
project(ReadSamFile)

set(CMAKE_CXX_STANDARD 14)

# add extra include directories
include_directories(/usr/local/lib/pybind11/include)
include_directories(/usr/local/lib/spdlog/include)
include_directories(/usr/include/python2.7)
include_directories(/usr/local/lib/seqan-seqan-v2.3.2/include)
include_directories(/usr/local/lib/boost_1_65_0)

# add boost library directories
set(BOOST_ROOT "/usr/local/lib/boost_1_65_0")
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
include_directories(${Boost_INCLUDE_DIRS})
find_package( Boost COMPONENTS system filesystem timer chrono REQUIRED )

#setup common source directories
set(SOURCE_FILES kgl_logging.h kgl_logging.cc kgl_read_sam.h
        kgl_mt_queue.h kgl_mt_data.h kgl_genome_types.h kgl_consume_sam.h
        kgl_process_sam.h kgl_lock.h kgl_mt_numpy.h kgl_nucleotide.h kgl_parse_sam.h
        kgl_parse_sam.cc kgl_mt_contig.h)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pg -pthread")
#set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wall -pg")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pthread")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wall")


set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../..)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../..)

#generate shared python library.
add_library(read_sam SHARED ${SOURCE_FILES} kgl_python.cc)

#generate kgl_genome executable
add_executable (kgl_genome ${SOURCE_FILES} kgl_main.cc kgl_exec_env.h kgl_exec_env.cc kgl_gff_fasta.h kgl_gff_fasta.cc kgl_main.cc kgl_genome.h kgl_genome.cc)
target_link_libraries(kgl_genome ${Boost_LIBRARIES})
