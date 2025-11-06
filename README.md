# Tools to analyse V1/V4 data

This datapipeline is used to preprocess, analyse and visualize LFP data from V1 and V4. 


Parameters for analyses are created using create_params.m. Processing is done using processing_pip_full.m. Analyses are stored in the analyses folder. 

The general process following pre-processing is the following: 

<img width="2000" height="1939" alt="process_overview" src="https://github.com/user-attachments/assets/caf1657f-076c-43ad-9c81-e0bd65af3d7f" />

Some more example plots: 

Instantaneous frequency over time 
<img width="1600" height="536" alt="grand_average_inst_freq" src="https://github.com/user-attachments/assets/fd670027-22e4-4071-afa9-0d1a2c47bcf9" />

Phase locking over time, grey areas indicate significant differences
<img width="1000" height="536" alt="V1dif_per_T_env_bpfilt" src="https://github.com/user-attachments/assets/a5df8957-4179-4d80-bff3-71e349743fd0" />
