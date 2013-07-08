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

## platform dependent compiler flags:
include(CheckCXXCompilerFlag)
if (NOT WIN32) # we only want fPIC on non-windows systems (fPIC is implicitly true there)
  CHECK_CXX_COMPILER_FLAG("-fPIC" WITH_FPIC)
  if (WITH_FPIC)
    add_definitions(-fPIC)
  endif()

  # Only add those definitions on non-windows systems
  add_definitions(-Wall -Wextra -Wno-non-virtual-dtor -Wno-long-long -Wno-variadic-macros)
  add_definitions(--pedantic)

  # Careful, this generates several hundred warnings ... 
  add_definitions(-Wconversion)
endif()

add_library(Crawdad ${my_sources})

add_executable(SimpleTest SimpleTest.C)
target_link_libraries(SimpleTest Crawdad)


configure_file(CrawdadConfig.cmake.in "${PROJECT_BINARY_DIR}/CrawdadConfig.cmake" @ONLY)

export(TARGETS Crawdad FILE CrawdadLibrary.cmake)

