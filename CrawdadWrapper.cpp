#include "CrawdadWrapper.h"

void CrawdadPeakFinder::SetChromatogram(IList<double>^ times, IList<double>^ intensities)
{
    // TODO: Check times to make sure they are evenly spaced

    // Marshall intensities to vector for Crawdad
    int len = intensities->Count;
    vector<float> intensitiesCrawdad(len);
double baselineIntensity = Double::MaxValue;
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
if (baselineIntensity == Double::MaxValue)
  baselineIntensity = 0;

    SetChromatogram(intensitiesCrawdad, maxIntensityIndex, baselineIntensity);
}

void CrawdadPeakFinder::SetChromatogram(IList<float>^ times, IList<float>^ intensities)
{
    // TODO: Check times to make sure they are evenly spaced

    // Marshall intensities to vector for Crawdad
    int len = intensities->Count;
    vector<float> intensitiesCrawdad(len);
double baselineIntensity = Double::MaxValue;
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
if (baselineIntensity == Double::MaxValue)
  baselineIntensity = 0;

    SetChromatogram(intensitiesCrawdad, maxIntensityIndex, baselineIntensity);
}

void CrawdadPeakFinder::SetChromatogram(vector<float>& intensities, int maxIntensityIndex, double baselineIntensity)
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
    FullWidthHalfMax = (float) fwhm;

_widthDataWings = (int)(FullWidthHalfMax*2);

if (_widthDataWings > 0)
{
  intensities.insert(intensities.begin(), _widthDataWings, (float)baselineIntensity);
  intensities.insert(intensities.end(), _widthDataWings, (float)baselineIntensity);
}

_pPeakFinder->clear();
    _pPeakFinder->set_chrom(intensities, 0);
}
