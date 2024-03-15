function [fractal, original,oscillatory,oscillatory_alt] = do_foof(trials)
%DO_FOOF Summary of this function goes here
%   Detailed explanation goes here
cfg               = [];
cfg.foilim        = [1 150]; 
cfg.pad           = 6;
cfg.tapsmofrq     = 2;
cfg.method        = 'mtmfft';
cfg.output        = 'fooof_aperiodic';
fractal = ft_freqanalysis(cfg, trials);
cfg.output        = 'pow';
original = ft_freqanalysis(cfg, trials);

% subtract the fractal component from the power spectrum
cfg               = [];
cfg.parameter     = 'powspctrm';
cfg.operation     = 'x2-x1';
oscillatory = ft_math(cfg, fractal, original);

% original implementation by Donoghue et al. 2020
cfg               = [];
cfg.parameter     = 'powspctrm';
cfg.operation     = 'x2./x1';  % equivalent to 10^(log10(x2)-log10(x1))
oscillatory_alt = ft_math(cfg, fractal, original);
end

