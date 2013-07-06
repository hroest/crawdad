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

std::vector<SlimCrawPeak> CrawdadWrapper::CalcPeaks(int max, std::vector<int> idIndices)
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

    // If max is not -1, then return the max most intense peaks, plus any
    // peaks that have been identified with MS/MS peptide search results
    if (max != -1)
    {
        // Shorten the list before performing the slow sort by intensity.
        // The sort shows up as bottleneck in a profiler.
        int lenResult = result.size();
        float intensityCutoff = 0;

        // TODO 
        throw "Not implemented";

#if 0
        FindIntensityCutoff(result, 0, (float)(totalArea/lenResult)*2, max, 1, intensityCutoff, lenResult);

	    List<KeyValuePair<CrawdadPeak^, bool>>^ resultFiltered =
            gcnew List<KeyValuePair<CrawdadPeak^, bool>>(lenResult);
        for (int i = 0, lenOrig = result->Count; i < lenOrig; i++)
        {
            CrawdadPeak^ peak = result[i];
            bool isIdentified = peak->IsIdentified(idIndices);
            if (isIdentified || peak->Area >= intensityCutoff || intensityCutoff == 0)
                resultFiltered->Add(KeyValuePair<CrawdadPeak^, bool>(peak, isIdentified));
        }

        resultFiltered->Sort(gcnew Comparison<KeyValuePair<CrawdadPeak^, bool>>(OrderIdAreaDesc));
        if (max < resultFiltered->Count)
            resultFiltered->RemoveRange(max, resultFiltered->Count - max);

        result = gcnew List<CrawdadPeak^>(resultFiltered->Count);
        for each (KeyValuePair<CrawdadPeak^, bool> peakId in resultFiltered)
            result->Add(peakId.Key);
#endif
    }

    return result;
}
