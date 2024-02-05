function [hilbertData,diffHilbertData,medfiltHilbert,instchange] = do_hilbert(trials,medianfiltord)
%DO_HILBERT Summary of this function goes here
%  


% Hilberting
cfg = [];
cfg.channel = 'all';
cfg.hilbert = 'angle';
hilbertData = ft_preprocessing(cfg,trials);
hilbertData.trial = cellfun(@unwrap,hilbertData.trial,'UniformOutput',false);

% Ableiting 
cfg = [];
cfg.channel = 'all';
cfg.absdiff = 'yes';
diffHilbertData = ft_preprocessing(cfg,hilbertData);
diffHilbertData.trial = cellfun(@(x) x*1000/(pi*2),diffHilbertData.trial,'UniformOutput',false);
% Old Median filtering 
%{
medfiltHilbert = diffHilbertData;
for ii = 1:length(diffHilbertData.trial)
    medfiltHilbert.trial = medfilt1(diffHilbertData.trial{ii},medfiltwidth);
end
end 
%}

%Median Filtering
cfg = [];
cfg.channel = 'all';
cfg.medianfilter = 'yes';
cfg.medianfiltord = medianfiltord;
medfiltHilbert = ft_preprocessing(cfg,diffHilbertData);

%{ 
%old inst change
cfg = [];
cfg.channel = 'all';
cfg.derivative = 'yes'
instchange = ft_preprocessing(cfg,medfiltHilbert);
%}
%
% five-point numerical derivation !!testing 9 right now!!
for ii = 1: length(medfiltHilbert.trial)
    instchange.diff{ii} = cent_diff_n(medfiltHilbert.trial{ii},1,5);
end 
%


