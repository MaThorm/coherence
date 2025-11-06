# Tools to analyse V1/V4 data

This datapipeline is used to preprocess, analyse and visualize LFP data from V1 and V4. 


Parameters for analyses are created using create_params.m. Processing is done using processing_pip_full.m. Analyses are stored in the analyses folder. 

The general process following pre-processing is the following: 

<img width="2000" height="1939" alt="process_overview" src="https://github.com/user-attachments/assets/caf1657f-076c-43ad-9c81-e0bd65af3d7f" />



## Some more example plots: 

### Instantaneous frequency over time 
<img width="1600" height="536" alt="grand_average_inst_freq" src="https://github.com/user-attachments/assets/fd670027-22e4-4071-afa9-0d1a2c47bcf9" />

### Instantaneous frequency over time, gray areas indicate significant differences between V1a and V1n
<img width="1000" height="536" alt="V1dif_per_T_env_bpfilt" src="https://github.com/user-attachments/assets/a5df8957-4179-4d80-bff3-71e349743fd0" />

### Average IF per recording site for different SSD bounds
<img width="840" height="642" alt="mean_freq_per_filtwidth_SSD" src="https://github.com/user-attachments/assets/684aedfe-6c12-497a-8e57-8f02a165e730" />

### Comparison of wavelet and hilbert approach
<img width="4000" height="856" alt="wave_inst_plot_in_t3" src="https://github.com/user-attachments/assets/b63fc197-0991-4400-9eae-bc1e80b4fe3b" />

### Performance of smoothing filters
<img width="2000" height="1039" alt="filt_compare_sess5_trial2" src="https://github.com/user-attachments/assets/e2c68905-f0b9-47fc-b72e-c382591a9938" />
