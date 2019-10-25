# ReCom
ReCom is a modified version of the [Comet](http://comet-ms.sourceforge.net/) project for the thorough analysis of post-translational modifications, developed at the [Cardiovascular Proteomics laboratory at CNIC (Spain)](https://www.cnic.es/).

The `.params` file for ReCom includes several new parameters:
* `use_delta_closest`: 0=no (default), 1=yes. Setting this to 1 will enable ReCom, setting this to 0 will perform a normal Comet-PTM search.
* `deltamass_peaks`: This is the path to the txt file containing the list of theoretical deltamasses.
* `deltamass_peaks_tolerance`: This is the tolerance in Da for matching experimental deltamasses to the list (default 1, will search for theoretical deltamasses within ±1 Da of the experimental deltamass.
* `cand_number`: This is the number of candidates to rescore for each identification (default 1)
* `use_xcorr_corr`: 0=no (default), 1=yes. Setting this to 1 will calculate the Corrected Xcorr value and use it to sort candidates and choose between experimental or theoretical. Regular Xcorr will still be calculated but not used for sorting.
