# Critical_Slowing_Epilepsy

Sample code from Maturana et al. 2020 - Critical slowing down as a biomarker for seizure susceptibility. 
Nature Communications

Three files are included that will generate some of the figures from the paper, and other figures that 
are similar but may vary due plots being made on a subset of data. 
Code was written in Matlab 2017b

Demo.m - will produce several results figures 
Demo_methods.m - will produce some of the methods figures
Demo_stats.m - will produce some of the population results

Access to the trial data is necessary to plot the remaining figures. Some data, including all seizures are
available on EpilepsyEcosystem.org. Additional data used in this study can be made available under a 
collaborative agreement and upon reasonable request to corresponding author

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Additional code has been included to demonstrate how features from iEEG data was computed. To run the code
you will need access to iEEG portal. 

CSD_compute_features.m - Downloads data in chunks of ~30 mins. For each chunk, the data is divided into
2-minute segments. In each segment, the variance and autocorrelation width is computed.

CSD_Forecasting.m - Uses the data obtained in CSD_compute_features.m. The code will compute filter the 
data into fast and slow components. Then the hilbert transform is used to estimate the signal phase. Seizure
times are used to determine the signal phase at the time of each seizure. Three predictions are computed:
1. Forecasting using the entire dataset, where two optimal thresholds are computed and data is thresholded 
based on these optimal thresholds
2. Forecasting in an iterative manner - thresholds are computed in an iterative manner. After each seizure, the
optimal threshold are recomputed and used until the next seizure. 
3. Random markov model - the transition probabilities (i.e from low to medium risk, medium to high etc) are 
computed. Using the transition probabilities, the risk is randomly computed at each time. Forecasts in methods
1 and 2 are compared to this random forecast



