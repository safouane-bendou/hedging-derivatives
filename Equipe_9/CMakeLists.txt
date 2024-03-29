cmake_minimum_required(VERSION 3.0)
project(mc-pricer CXX)

# Prefix
if (NOT PREFIX)
    set(PREFIX ${CMAKE_CURRENT_BINARY_DIR})
endif (NOT PREFIX)
set(CMAKE_INSTALL_PREFIX ${PREFIX} CACHE STRING "Installation prefix." FORCE)

if (CXX)
    set(CMAKE_CXX_COMPILER ${CXX})
endif (CXX)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Release or Debug
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif ()


# Detect PNL
find_package(Pnl REQUIRED)
set(LIBS ${LIBS} ${PNL_LIBRARIES})
include_directories(${PNL_INCLUDE_DIRS})
if(${CMAKE_BUILD_TYPE} STREQUAL "Release")
    add_definitions(-DPNL_RANGE_CHECK_OFF)
endif()

# Add jlparser
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if (WIN32)
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("/std:c++latest" _cpp_latest_flag_supported)
    if (_cpp_latest_flag_supported)
        add_compile_options("/std:c++latest")
    endif()
endif()
add_subdirectory(src/jlparser)

### Ne rien modifier dans la partie au dessus de cette ligne.

# Test the parser
add_executable(test-parser src/test_parser.cpp)
target_link_libraries(test-parser jl_parser ${LIBS})

add_executable(mc-example src/mc_example.cpp)
target_link_libraries(mc-example jl_parser ${LIBS})

add_executable(TestProject src/Tests/TestProject.cpp src/BlackScholesModel.cpp src/AsianOption.cpp src/BasketOption.cpp src/MonteCarlo.cpp src/PerformanceOption.cpp)
target_link_libraries(TestProject jl_parser ${LIBS})

add_executable(TestAsian src/Tests/TestAsian.cpp src/BlackScholesModel.cpp src/AsianOption.cpp src/BasketOption.cpp src/MonteCarlo.cpp src/PerformanceOption.cpp)
target_link_libraries(TestAsian jl_parser ${LIBS})

add_executable(TestBasket src/Tests/TestBasket.cpp src/BlackScholesModel.cpp src/AsianOption.cpp src/BasketOption.cpp src/MonteCarlo.cpp src/PerformanceOption.cpp)
target_link_libraries(TestBasket jl_parser ${LIBS})

add_executable(TestBasket1 src/Tests/TestBasket1.cpp src/BlackScholesModel.cpp src/AsianOption.cpp src/BasketOption.cpp src/MonteCarlo.cpp src/PerformanceOption.cpp)
target_link_libraries(TestBasket1 jl_parser ${LIBS})

add_executable(TestBasket2 src/Tests/TestBasket2.cpp src/BlackScholesModel.cpp src/AsianOption.cpp src/BasketOption.cpp src/MonteCarlo.cpp src/PerformanceOption.cpp)
target_link_libraries(TestBasket2 jl_parser ${LIBS})

add_executable(TestPerformance src/Tests/TestPerformance.cpp src/BlackScholesModel.cpp src/AsianOption.cpp src/BasketOption.cpp src/MonteCarlo.cpp src/PerformanceOption.cpp)
target_link_libraries(TestPerformance jl_parser ${LIBS})

add_executable(TestBasket1FullTime src/Tests/TestBasket1FullTime.cpp src/BlackScholesModel.cpp src/AsianOption.cpp src/BasketOption.cpp src/MonteCarlo.cpp src/PerformanceOption.cpp)
target_link_libraries(TestBasket1FullTime jl_parser ${LIBS})

add_executable(TestSimulatedDataBasket1 src/Tests/TestSimulatedDataBasket1.cpp src/BlackScholesModel.cpp src/AsianOption.cpp src/BasketOption.cpp src/MonteCarlo.cpp src/PerformanceOption.cpp)
target_link_libraries(TestSimulatedDataBasket1 jl_parser ${LIBS})

add_executable(TestPnl src/Tests/TestPnl.cpp src/BlackScholesModel.cpp src/AsianOption.cpp src/BasketOption.cpp src/MonteCarlo.cpp src/PerformanceOption.cpp)
target_link_libraries(TestPnl jl_parser ${LIBS})

# Print compiler flags
get_directory_property(DirDefs COMPILE_DEFINITIONS)
message(STATUS "COMPILE_DEFINITIONS = ${DirDefs}")
if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_LINKER_EXE_FLAGS}" )
else()
    message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_LINKER_EXE_FLAGS}" )
endif()
