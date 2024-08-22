function [results] = do_fft(full_trials,params)
%DO_FOOF Calculates FFT
for ii = 1:length(full_trials)
    trials = full_trials(ii);
    cfg               = [];
    cfg.foilim        = [1 150];
    cfg.channel       = trials.label{1,1};
    cfg.pad           = 6;
    cfg.tapsmofrq     = 2;
    cfg.method        = 'mtmfft';
    cfg.output        = 'pow';
    cfg.channel       = 'all';
    results(ii) = ft_freqanalysis(cfg, trials);
end 
end 