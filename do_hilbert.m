function [phase_angles_wr] = do_hilbert(trials)
%DO_HILBERT Summary of this function goes here
% Does the Hilbert transform on a fieldtrip file. Calculates the phase
% angles and the unwrapped phase angles 
% Inputs should be: 
% Trials: result of pre_processing_pip_trials


% Hilberting
cfg = [];
cfg.channel = 'all'
cfg.hilbert = 'angle';
phase_angles_wr = ft_preprocessing(cfg,trials);



