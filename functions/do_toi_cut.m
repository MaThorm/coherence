function [new_trials] = do_toi_cut(trials,toilim)
%Cuts the ft_struct into the toi of interest

for i_s = 1:length(trials)
    cfg = [];
    cfg.toilim = toilim;
    new_trials(i_s) = ft_redefinetrial(cfg, trials(i_s));
end

end 

