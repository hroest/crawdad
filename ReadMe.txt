========================================================================
    DYNAMIC LINK LIBRARY : Crawdad Project Overview
========================================================================

AppWizard has created this Crawdad DLL for you.  

This file contains a summary of what you will find in each of the files that
make up your Crawdad application.

Crawdad.vcproj
    This is the main project file for VC++ projects generated using an Application Wizard. 
    It contains information about the version of Visual C++ that generated the file, and 
    information about the platforms, configurations, and project features selected with the
    Application Wizard.

CMakeLists.txt
    This is the main project file for CMake based projects (See below "Using CMake")

Crawdad.cpp
    This is the main DLL source file.

Crawdad.h
    This file contains a class declaration.

AssemblyInfo.cpp
	Contains custom attributes for modifying assembly metadata.

/////////////////////////////////////////////////////////////////////////////
Other notes:

AppWizard uses "TODO:" to indicate parts of the source code you
should add to or customize.

/////////////////////////////////////////////////////////////////////////////
// Using CMake 
/////////////////////////////////////////////////////////////////////////////

This is a copy of Crawdad from the proteowizard project, taken at revision 4250
using cmake to build.

== Build and Install ==
Simply type

$ cmake .
$ make

== Build against ==

To build against this, use the following code in your CMakeLists

	find_package(Crawdad)
	include_directories(${CRAWDAD_INCLUDE_DIRS})

and then add the library as follows to your executables:

	target_link_libraries(myExecutable Crawdad)

If cmake does not find the library, try to set the path to the place where you
build Crawdad using the Crawdad_DIR parameter, for example:

$ cmake -DCrawdad_DIR=/path/to/crawdad .

In case you have a more complicated setup and would like to identify the actual
path of the libraries, you can use

	find_library(Crawdad_LIBRARY NAMES Crawdad.a Crawdad HINTS ${Crawdad_DIR})

which will put the library name into the "Crawdad_LIBRARY" variable. 

