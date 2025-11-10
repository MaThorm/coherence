function [trials] = do_bpfilt(trials,params)
%Does a bpfilt to all structures 
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = params.bpwidth;
cfg.padding = 4;
cfg.bpfilttype = 'fir'; 
cfg.bpfiltord = params.order;
for i_s = 1:length(trials)
    trials(i_s) = ft_preprocessing(cfg,trials(i_s));
end

end 
