/*
 * Original author: Hannes Roest <roest .at. imsb.biol.ethz.ch>,
 *                  Aebersold Lab, IMSB, ETH Zurich
 *
 * Copyright 2013 ETH Zurich
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "CrawPeak.H"
#include "CrawPeakFinder.H"
#include "CrawPeakMethod.H"

class CrawdadWrapper
{
  private:
		int _widthDataWings;
    crawpeaks::StackCrawPeakFinder * _pPeakFinder;

  public:

    /** Constructor
     *
     * Initializes the PeakFinder and sets a few defaults.
    */
    CrawdadWrapper()
    {
        _pPeakFinder = new crawpeaks::StackCrawPeakFinder();
        _pPeakFinder->slim = true;
        _pPeakFinder->method.peak_location_meth = crawpeaks::MAXIMUM_PEAK;
        _pPeakFinder->method.background_estimation_method = crawpeaks::LOWER_BOUNDARY;

        /*
         * The following defaults are set in the Crawdad algorithm by default
         *
          this->set_sd(4.0f);
          saved_weights = false;
          mean_cutoff   = false;
          minimum_level = 0.0f;
          min_len       = 3;
          switch_len    = 2;
          bound_meth = GAUSS_2D_BOUND;
          peak_location_meth = MAXIMUM_PEAK;
          chrom_smooth_method = SAVITZKY_GOLAY_SMOOTHER;
          extend_peak_to_lower_bound = false;
          extend_to_zero_crossing = false;
          extend_allowed_asymmetry = 1.0f;
          fraction_to_valley = 0.0f;
          background_estimation_method = PEAK_BOUNDARY_ESTIMATE;
          this->one_peak_slope_merge_constraint = 0.0f;
          this->mean_slope_merge_constraint = 0.0f;
          this->exclude_extension_overlaps_by_peakrt = false;
          this->ratchet_back_to_frac_maxval = -1.0f;
          this->extend_from_peak_rt = false;
          merge_peaks_list_based = false;
          extend_peak_set = false;
        */

    }

    /** Destructor
    */
    ~CrawdadWrapper()
    {
    }

    /** Getter / Setter for FWHM
     *
     * calls the underlying fwhm methods of the PeakFinder
    */
    float get_fwhm() { return _pPeakFinder->method.get_fwhm(); }
    void set_fwhm(float fwhm)
    {
      _pPeakFinder->method.set_fwhm(fwhm*3);
      _pPeakFinder->method.min_len = (int)(fwhm/4.0 + 0.5);
    }

    /** Getter / Setter for stdev
     *
     * calls the underlying stdev methods of the PeakFinder
    */
    float get_stdev() { return _pPeakFinder->method.get_sd(); }
    void set_stdev(float sd) { _pPeakFinder->method.set_sd(sd); }

    /** Initialize with a float-chromatogram
     *
     * Determines the maximal and minimal intensities and then calls SetChromatogram
    */
    void SetChromatogram(std::vector<float>& times, std::vector<float> intensities);

    /** Initialize with a double-chromatogram
     *
     * Determines the maximal and minimal intensities and then calls SetChromatogram
    */
    void SetChromatogram(std::vector<double>& times, std::vector<double> intensities);

    /** Execute the algorithm
     *
     * Calls call_peaks of the PeakFinder, then applies some Skyline-specific
     * post-processing on the result before returning it.
    */
    std::vector<crawpeaks::SlimCrawPeak> CalcPeaks();

    /** Get the intensities of the 2nd derivative
     *
    */
    void getIntensities2d(std::vector<float>& intensities2d);

    /** Get the intensities of the 1st derivative
     *
    */
    void getIntensities1d(std::vector<float>& intensities1d);

    /** Get a peak between the start and end indices provided
     *
    */
    crawpeaks::SlimCrawPeak GetPeak(int startIndex, int endIndex);

  private:
    void SetChromatogram(std::vector<float>& intensities, int maxIntensityIndex, double baselineIntensity);


};

