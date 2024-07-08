function [phase_angles,inst_freq,filt_data] = do_filtinst(phase_angles_wr,filttype,framelen,filtord)
%DO_HILBERT Summary of this function goes here
% Does the Hilbert Transform on 
% Inputs should be: 
% Trials: result of pre_processing_pip_trials
% filttype: either medfilt ord sgolay
% framelen: length of filter frame
% filtord: order of filter, only relevant for sgolay

% Unwrapping
phase_angles = phase_angles_wr
phase_angles.trial = cellfun(@(x) unwrap(x,[],2),phase_angles.trial,'UniformOutput',false);


if filttype == "sgolay"
    filt_data = phase_angles;
    for ii = 1:length(phase_angles.trial)
        for i_t = 1:size(filt_data.trial{1,ii},1)
            filt_data.trial{1,ii}(i_t,:) = sgolayfilt(filt_data.trial{1,ii}(i_t,:),filtord,framelen);
        end
    end 
 

    % Ableiting 
    cfg = [];
    cfg.channel = 'all'
    
    %cfg.absdiff = 'yes';
    cfg.derivative = 'yes';
    inst_freq = ft_preprocessing(cfg,filt_data);
    inst_freq.trial = cellfun(@(x) x*1000/(pi*2),inst_freq.trial,'UniformOutput',false);

    % five-point numerical derivation !!testing 9 right now!!
    for ii = 1: length(inst_freq.trial)
        instchange.diff{ii} = cent_diff_n(inst_freq.trial{ii},1,5);
    end 
elseif filttype == "medfilt"
    % Ableiting, in the case of medfilt first
    cfg = [];
    cfg.channel = 'all'
    %cfg.absdiff = 'yes';
    cfg.derivative = 'yes';
    inst_freq = ft_preprocessing(cfg,phase_angles);
    inst_freq.trial = cellfun(@(x) x*1000/(pi*2),inst_freq.trial,'UniformOutput',false);

    % fieldtrip way of median filtering 
    cfg = [];
    cfg.channel = 'all'
    cfg.medianfilter = 'yes';
    cfg.medianfiltord = framelen;
    cfg.padding = 2.5;
    cfg.padtype = 'mirror';
    filt_data = ft_preprocessing(cfg,inst_freq);
%   New way of median filtering
%     filt_data = inst_freq;
%     for ii = 1:length(inst_freq.trial)
%         for i_t = 1:size(filt_data.trial{1,ii},1)
%             filt_data.trial{1,ii}(i_t,:) = medfilt1(filt_data.trial{1,ii}(i_t,:),framelen,[],2,"includenan",'truncate');
%         end
%     end 
    % five-point numerical derivation !!testing 9 right now!!
for ii = 1: length(filt_data.trial)
    instchange.diff{ii} = cent_diff_n(filt_data.trial{ii},1,5);
end 
end  
end

