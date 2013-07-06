### example CMakeLists.txt to develop programs using OpenMS
project("Main")
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
)

add_library(Crawdad ${my_sources})

add_executable(SimpleTest SimpleTest.C)
target_link_libraries(SimpleTest Crawdad)

