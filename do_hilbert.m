function [hilbertData_wr,hilbertData,diffHilbertData,filtHilbert,instchange] = do_hilbert(trials,filttype,framelen,filtord)
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
hilbertData_wr = ft_preprocessing(cfg,trials);
%hilbertData.trial =cellfun(@unwrap,hilbertData.trial,'UniformOutput',false); old cellfun

% Unwrapping
hilbertData = hilbertData_wr;
hilbertData.trial = cellfun(@(x) unwrap(x,[],2),hilbertData.trial,'UniformOutput',false);

% if filtertype == sgolay 
if filttype == 'sgolay'
     filtHilbert = hilbertData;
    for ii = 1:length(hilbertData.trial)
        for i_t = 1:size(filtHilbert.trial{1,ii},1)
            filtHilbert.trial{1,ii}(i_t,:) = sgolayfilt(hilbertData.trial{1,ii}(i_t,:),filtord,framelen);
        end
    end 
 

    % Ableiting 
    cfg = [];
    cfg.channel = 'all'
    
    %cfg.absdiff = 'yes';
    cfg.derivative = 'yes';
    diffHilbertData = ft_preprocessing(cfg,filtHilbert);
    diffHilbertData.trial = cellfun(@(x) x*1000/(pi*2),diffHilbertData.trial,'UniformOutput',false);

    % five-point numerical derivation !!testing 9 right now!!
    for ii = 1: length(diffHilbertData.trial)
        instchange.diff{ii} = cent_diff_n(diffHilbertData.trial{ii},1,5);
    end 
elseif filttype == 'medfilt'
    % Ableiting, in the case of medfilt first
    cfg = [];
    cfg.channel = 'all'
    %cfg.absdiff = 'yes';
    cfg.derivative = 'yes';
    diffHilbertData = ft_preprocessing(cfg,hilbertData);
    diffHilbertData.trial = cellfun(@(x) x*1000/(pi*2),diffHilbertData.trial,'UniformOutput',false);

    %Filtering somehow
    cfg = [];
    cfg.channel = 'all'
    cfg.medianfilter = 'yes';
    cfg.medianfiltord = framelen;
    filtHilbert = ft_preprocessing(cfg,diffHilbertData);
    
    % five-point numerical derivation !!testing 9 right now!!
    for ii = 1: length(filtHilbert.trial)
        instchange.diff{ii} = cent_diff_n(filtHilbert.trial{ii},1,5);
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

%


