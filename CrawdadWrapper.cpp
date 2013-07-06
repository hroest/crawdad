#include "CrawdadWrapper.h"

using namespace crawpeaks;

void CrawdadWrapper::SetChromatogram(vector<double>& times, vector<double> intensities)
{
    // TODO: Check times to make sure they are evenly spaced

    // Marshall intensities to vector for Crawdad
    int len = intensities.size();
    vector<float> intensitiesCrawdad(len);
    double baselineIntensity = numeric_limits<double>::max();
    double maxIntensity = 0;
    int maxIntensityIndex = -1;
    for (int i = 0; i < len; i++)
    {
        float intensity = (float)intensities[i];
        intensitiesCrawdad[i] = intensity;

        // Keep track of where the maximum intensity is
        if (intensity > maxIntensity)
        {
            maxIntensity = intensity;
            maxIntensityIndex = i;
        }
	if (intensity < baselineIntensity)
		baselineIntensity = intensity;
    }
if (baselineIntensity == numeric_limits<double>::max())
	baselineIntensity = 0;

    SetChromatogram(intensitiesCrawdad, maxIntensityIndex, baselineIntensity);
}

void CrawdadWrapper::SetChromatogram(vector<float>& times, vector<float> intensities)
{
    // TODO: Check times to make sure they are evenly spaced

    // Marshall intensities to vector for Crawdad
    int len = intensities.size();
    vector<float> intensitiesCrawdad(len);
    double baselineIntensity = numeric_limits<double>::max();
    double maxIntensity = 0;
    int maxIntensityIndex = -1;
    for (int i = 0; i < len; i++)
    {
        float intensity = intensities[i];
        intensitiesCrawdad[i] = intensity;

        // Keep track of where the maximum intensity is
        if (intensity > maxIntensity)
        {
            maxIntensity = intensity;
            maxIntensityIndex = i;
        }
	if (intensity < baselineIntensity)
		baselineIntensity = intensity;
    }
if (baselineIntensity == numeric_limits<double>::max())
	baselineIntensity = 0;

    SetChromatogram(intensitiesCrawdad, maxIntensityIndex, baselineIntensity);
}

void CrawdadWrapper::SetChromatogram(vector<float>& intensities, int maxIntensityIndex, double baselineIntensity)
{
    // Find the peak width of the maximum intensity point at
    // half its height.
    int fwhm = 6;
    if (maxIntensityIndex != -1)
    {
        double halfHeight = (intensities[maxIntensityIndex] - baselineIntensity)/2 + baselineIntensity;
        int iStart = 0;
        for (int i = maxIntensityIndex - 1; i >= 0; i--)
        {
            if (intensities[i] < halfHeight)
            {
                iStart = i;
                break;
            }
        }
        int len = (int)intensities.size();
        int iEnd = len - 1;
        for (int i = maxIntensityIndex + 1; i < len; i++)
        {
            if (intensities[i] < halfHeight)
            {
                iEnd = i;
                break;
            }
        }
        fwhm = max(fwhm, iEnd - iStart);
    }
    set_fwhm(fwhm);

_widthDataWings = (int)(get_fwhm()*2);

if (_widthDataWings > 0)
{
  intensities.insert(intensities.begin(), _widthDataWings, (float)baselineIntensity);
  intensities.insert(intensities.end(), _widthDataWings, (float)baselineIntensity);
}

_pPeakFinder->clear();
    _pPeakFinder->set_chrom(intensities, 0);
}

std::vector<SlimCrawPeak> CrawdadWrapper::CalcPeaks()
{
    // Find peaks
    _pPeakFinder->call_peaks();

    // Marshall found peaks to managed list
    std::vector<SlimCrawPeak> result;
    result.reserve((int)_pPeakFinder->sps.size());

    vector<SlimCrawPeak>::iterator itPeak = _pPeakFinder->sps.begin();
    vector<SlimCrawPeak>::iterator itPeakEnd = _pPeakFinder->sps.end();
    double totalArea = 0;
int stop_rt = (int)_pPeakFinder->chrom.size() - _widthDataWings - 1;
int adjust_stop_rt = stop_rt - _widthDataWings;
    while (itPeak != itPeakEnd)
    {
	if (itPeak->start_rt_idx < stop_rt && itPeak->stop_rt_idx > _widthDataWings)
	{
		double rheight = itPeak->peak_height / itPeak->raw_height;
		double rarea = itPeak->peak_area / itPeak->raw_area;

		if (rheight > 0.02 && rarea > 0.02)
		{
      itPeak->start_rt_idx = std::max(_widthDataWings, itPeak->start_rt_idx);
      itPeak->start_rt_idx -= _widthDataWings;
      itPeak->peak_rt_idx = std::max(_widthDataWings, std::min(stop_rt, itPeak->peak_rt_idx));
      itPeak->peak_rt_idx -= _widthDataWings;
      itPeak->stop_rt_idx = std::max(_widthDataWings, std::min(stop_rt, itPeak->stop_rt_idx));
      itPeak->stop_rt_idx -= _widthDataWings;

      result.push_back(*itPeak);

			totalArea += itPeak->peak_area;
		}
	}
        itPeak++;
    }
    return result;
}

void CrawdadWrapper::getIntensities2d(std::vector<float>& intensities2d)
{
    // Make sure the 2d chromatogram is populated
    if (_pPeakFinder->chrom_2d.size() != _pPeakFinder->chrom.size())
    {
        _pPeakFinder->chrom_2d.resize(_pPeakFinder->chrom.size());
        _pPeakFinder->get_2d_chrom(_pPeakFinder->chrom, _pPeakFinder->chrom_2d);
    }

    vector<float>::iterator it = _pPeakFinder->chrom_2d.begin() + _widthDataWings;
    vector<float>::iterator itEnd = _pPeakFinder->chrom_2d.end() - _widthDataWings;
    while (it < itEnd)
    {
        intensities2d.push_back(*it);
        it++;
    }
}

void CrawdadWrapper::getIntensities1d(std::vector<float>& intensities1d)
{
    // Make sure the 2d chromatogram is populated
    if (_pPeakFinder->chrom_1d.size() != _pPeakFinder->chrom.size())
    {
        _pPeakFinder->chrom_1d.resize(_pPeakFinder->chrom.size());
        _pPeakFinder->get_1d_chrom(_pPeakFinder->chrom, _pPeakFinder->chrom_1d);
    }

    vector<float>::iterator it = _pPeakFinder->chrom_1d.begin() + _widthDataWings;
    vector<float>::iterator itEnd = _pPeakFinder->chrom_1d.end() - _widthDataWings;
    while (it < itEnd)
    {
        intensities1d.push_back(*it);
        it++;
    }
}

crawpeaks::SlimCrawPeak CrawdadWrapper::GetPeak(int startIndex, int endIndex)
{
startIndex += _widthDataWings;
endIndex += _widthDataWings;
    
crawpeaks::SlimCrawPeak peak;
    _pPeakFinder->annotator.reannotate_peak(peak, startIndex, endIndex);
_pPeakFinder->annotator.calc_fwhm(peak);

peak.start_rt_idx -= _widthDataWings;
peak.stop_rt_idx -= _widthDataWings;
peak.peak_rt_idx -= _widthDataWings;

return peak;
}

