#include <stdio.h>

#include "CrawPeak.H"
#include "CrawPeakFinder.H"

using namespace crawpeaks;

int main( int argc, const char* argv[] )
{

	SlimCrawPeak* peak = new SlimCrawPeak();
  delete peak;

  StackCrawPeakFinder* finder = new StackCrawPeakFinder();
  delete finder;

  return 0;
}

