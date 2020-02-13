# Critical_Slowing_Epilepsy

The code provided shows how to compute the autocorrelation and variance from iEEG data. Code assumes you
have access to iEEG portal and are able to access data from there. 

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

For more details see Maturana et al. 2020 - Critical slowing down as a biomarker for seizure susceptibility. 
Nature Communications