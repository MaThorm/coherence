function [hilbertData,diffHilbertData,medfiltHilbert,instchange] = do_hilbert(trials,filttype,framelen,filtord)
%DO_HILBERT Summary of this function goes here
% Does the Hilbert transform on a fieldtrip file. Using the resulting angle
% the instantaneous frequency is calcualted by taking the derivative. The
% resulting frequency is then filtered. 
% Inputs should be: 
% Trials: result of pre_processing_pip_trials
% filttype: either medfilt ord sgolay
% framelen: length of filter frame
% filtord: order of filter, only relevant for sgolay 


% Hilberting
cfg = [];
cfg.channel = 'all'
cfg.hilbert = 'angle';
hilbertData = ft_preprocessing(cfg,trials);
%hilbertData.trial =cellfun(@unwrap,hilbertData.trial,'UniformOutput',false); old cellfun
hilbertData.trial = cellfun(@(x) unwrap(x,[],2),hilbertData.trial,'UniformOutput',false);



% Ableiting 
cfg = [];
cfg.channel = trials.label{1}
cfg.absdiff = 'yes';
diffHilbertData = ft_preprocessing(cfg,hilbertData);
diffHilbertData.trial = cellfun(@(x) x*1000/(pi*2),diffHilbertData.trial,'UniformOutput',false);


%Filtering somehow
if filttype == 'medfilt'
    cfg = [];
    cfg.channel = trials.label{1}
    cfg.medianfilter = 'yes';
    cfg.medianfiltord = framelen;
    medfiltHilbert = ft_preprocessing(cfg,diffHilbertData);
elseif filttype == 'sgolay'
    medfiltHilbert = diffHilbertData;
    for ii = 1:length(diffHilbertData.trial)
        medfiltHilbert.trial{1,ii} = sgolayfilt(diffHilbertData.trial{1,ii},filtord,framelen);
    end 
end  

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


