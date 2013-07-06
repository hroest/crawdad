
#include "CrawPeak.H"
#include "CrawPeakFinder.H"
#include "CrawPeakMethod.H"

class CrawdadWrapper
{
  private:
		int _widthDataWings;
    crawpeaks::StackCrawPeakFinder * _pPeakFinder;

  public:

    CrawdadWrapper()
    {
        _pPeakFinder = new crawpeaks::StackCrawPeakFinder();
        _pPeakFinder->slim = true;
        _pPeakFinder->method.peak_location_meth = crawpeaks::MAXIMUM_PEAK;
        _pPeakFinder->method.background_estimation_method = crawpeaks::LOWER_BOUNDARY;
    }

    ~CrawdadWrapper()
    {
    }

    float get_fwhm() { return _pPeakFinder->method.get_fwhm(); }
    void set_fwhm(float fwhm)
    {
      _pPeakFinder->method.set_fwhm(fwhm*3);
      _pPeakFinder->method.min_len = (int)(fwhm/4.0 + 0.5);
    }

    float get_stdev() { return _pPeakFinder->method.get_sd(); }
    void set_stdev(float sd) { _pPeakFinder->method.set_sd(sd); }

    void SetChromatogram(vector<float>& times, vector<float> intensities);

    void SetChromatogram(vector<double>& times, vector<double> intensities);

    std::vector<crawpeaks::SlimCrawPeak> CalcPeaks();

    void getIntensities2d(std::vector<float>& intensities2d);

    void getIntensities1d(std::vector<float>& intensities1d);

    crawpeaks::SlimCrawPeak GetPeak(int startIndex, int endIndex);

  private:
    void SetChromatogram(vector<float>& intensities, int maxIntensityIndex, double baselineIntensity);


};

