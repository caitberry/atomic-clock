Paper re-organized as of 3/6/25 to not include details of data preparation 
Need to communicate to scientists:
    1) How to exclude data marked as bad quality from analysis 

Having obtained clock shift files that include the three variables $t_{MJD}$, $f_{shift}$, and $IS_{GOOD}$, a comb data file containing $t_{MJD}$ and any other variables necessary to compute the optical frequency for each clock, and ensuring the values of the variables are high precision decimal type, the first step in processing the data for analysis is to exclude any clock observations where $IS_{GOOD}=0$. 

For example, suppose the $Al+$ clock data is stored in a data frame called *shift_data_Al*. In order to analyze only the data that is marked as reliable by the clock scientist, we can use the following Python code to create a subset containing only the good data, named *shift_data_Al_good*. 

```{python, eval=FALSE, echo=TRUE}
good_condition_al = shift_data_Al["IS_GOOD"] == 1
shift_data_Al_good = shift_data_Al[good_condition_al].reset_index(drop=True)
```

    2) Each time series (comb-clock comparison values, shift values) are irregularly sampled 

    Due to the nature of the machinery involved in generating clock and comb data, any time series analyzed will be an irregularly observed one. That is, the interval between consecutive measurement times $t_{i}$ and $t_{i+1}$, where $t_{i}$ is a particular MDJ value, is non-constant. This is true even for comb data, though in a much less drastic way. For example, most intervals for data from the ErYb comb are around $t_{i+1} - t_{i} = 0.000012$ however some observational gaps can be almost double that size with, say, $t_{j+1} - t_{j} = 0.000022$. In any case, these discrepancies in the time at which the frequency data is collected result in an irregularly sampled time series which poses its own challenges in analysis. Typically, such irregularities are corrected by interpolation which is the topic of a later section.  

    3) Distinction between missing vs irregular sampling remains ambiguous, especially when taking into consideration filtering and deglitching processes


    