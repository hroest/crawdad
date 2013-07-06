#include <stdio.h>
#include <assert.h>

#include "CrawPeak.H"
#include "CrawPeakFinder.H"
#include "CrawdadWrapper.h"

using namespace crawpeaks;

void test_data()
{
  printf( "\nStarting test with real data\n" );

  // this is a simulated SRM experiment where the two traces are not sampled at
  // the exact same time points, thus a resampling is necessary before applying
  // the algorithm.
  static const double rtdata_1[] = {1474.34, 1477.11, 1479.88, 1482.64, 1485.41, 1488.19, 1490.95, 1493.72, 1496.48, 1499.25, 1502.03, 1504.8 , 1507.56, 1510.33, 1513.09, 1515.87, 1518.64, 1521.42};
  static const double rtdata_2[] = {1473.55, 1476.31, 1479.08, 1481.84, 1484.61, 1487.39, 1490.15, 1492.92, 1495.69, 1498.45, 1501.23, 1504   , 1506.76, 1509.53, 1512.29, 1515.07, 1517.84, 1520.62};

  static const double intdata_1[] = {3.26958, 3.74189, 3.31075, 86.1901, 3.47528, 387.864, 13281  , 6375.84, 39852.6, 2.66726, 612.747, 3.34313, 793.12 , 3.29156, 4.00586, 4.1591 , 3.23035, 3.90591};
  static const double intdata_2[] =  {3.44054 , 2142.31 , 3.58763 , 3076.97 , 6663.55 , 45681   , 157694  , 122844  , 86034.7 , 85391.1 , 15992.8 , 2293.94 , 6934.85 , 2735.18 , 459.413 , 3.93863 , 3.36564 , 3.44005};
  int len = 18;


  {
  std::vector<double> intensities;
  intensities.assign(intdata_1, intdata_1 + len);

  std::vector<double> rts;
  rts.assign(rtdata_1, rtdata_1 + len);

  CrawdadWrapper wr;
  wr.SetChromatogram(rts, intensities);
  std::vector<crawpeaks::SlimCrawPeak> res = wr.CalcPeaks();

  assert(res.size() == 1);
  assert(res[0].start_rt_idx == 2);
  assert(res[0].peak_rt_idx == 8);
  assert(res[0].stop_rt_idx == 13);

  for (size_t i = 0; i < res.size(); i++)
  {
    std::cout << " Peak here " << std::endl;

    std::cout << " == found peak at " << res[i].start_rt_idx << " and " << res[i].peak_rt_idx << " to " << res[i].stop_rt_idx << std::endl;
    std::cout << " == which is  " << rts[ res[i].start_rt_idx ] << " and " << rts[ res[i].peak_rt_idx ] << " to " << rts[ res[i].stop_rt_idx ] << std::endl;
  }
  }

  {
  std::vector<double> intensities;
  intensities.assign(intdata_2, intdata_2 + len);

  std::vector<double> rts;
  rts.assign(rtdata_2, rtdata_2 + len);

  CrawdadWrapper wr;
  wr.SetChromatogram(rts, intensities);
  std::vector<crawpeaks::SlimCrawPeak> res = wr.CalcPeaks();

  assert(res.size() == 1);
  assert(res[0].start_rt_idx == 2);
  assert(res[0].peak_rt_idx == 6);
  assert(res[0].stop_rt_idx == 13);

  for (size_t i = 0; i < res.size(); i++)
  {
    std::cout << " == found peak at " << res[i].start_rt_idx << " and " << res[i].peak_rt_idx << " to " << res[i].stop_rt_idx << std::endl;
    std::cout << " == which is  " << rts[ res[i].start_rt_idx ] << " and " << rts[ res[i].peak_rt_idx ] << " to " << rts[ res[i].stop_rt_idx ] << std::endl;
  }
  }

}

void test_creation()
{
  printf( "\nStarting to create objects\n" );

	SlimCrawPeak* peak = new SlimCrawPeak();
  delete peak;

  StackCrawPeakFinder* finder = new StackCrawPeakFinder();
  delete finder;

  CrawdadWrapper* wrapper = new CrawdadWrapper();
  delete wrapper;
}

int main(int /* argc */, const char** /* argv */)
{

  test_creation();
  test_data();

  return 0;
}

