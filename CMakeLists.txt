project("Crawdad")
cmake_minimum_required(VERSION 2.6)

include_directories(msmat)

set(my_sources
# AssemblyInfo.cpp
# ManagedCrawdad.cpp
CrawPeakAnnotator.cpp  
CrawPeak.cpp  
CrawPeakFinder.cpp  
CrawPeakMethod.cpp  
msmat/crawutils.cpp
filters/ChromSmoother.cpp
filters/GaussSmoother.cpp
CrawdadWrapper.cpp  
)

add_library(Crawdad ${my_sources})

add_executable(SimpleTest SimpleTest.C)
target_link_libraries(SimpleTest Crawdad)


configure_file(CrawdadConfig.cmake.in "${PROJECT_BINARY_DIR}/CrawdadConfig.cmake" @ONLY)

export(TARGETS Crawdad FILE CrawdadLibrary.cmake)

