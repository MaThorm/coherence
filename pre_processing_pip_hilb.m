function [hilbert_angles_wr,hilbert_angles,inst_freq,filt_inst_freq,instchange] = pre_processing_pip_hilb(trials,filttype,framelen,filtord)
%PRE_PROCESSING_PIP_HILB Gives back multiple structures relating to the
%instantaneous frequency. Angles: Result of Hilbert Transform, inst_freq:
%derivative of angles, filt_inst_freq: filtered inst_freq, instchange:
%instantaneous change, summary: means & SDs
% Input should be
% trials: a trialselected fieldtrip struct
% filttype: either 'medfilt' or 'sgolay'
% framelen: framelength of the filter
% filtord: order of the filter applied (only for sgolay)

parfor ii = 1:length(trials)
    [hilbert_angles_wr(ii),hilbert_angles(ii), inst_freq(ii), filt_inst_freq(ii), instchange(ii)] = do_hilbert(trials(ii),filttype,framelen,filtord); 
end 

% mean over trials 
for ii = 1:length(filt_inst_freq)
    aMat = cell2matnan(filt_inst_freq(ii).trial,1);
    filt_inst_freq(ii).sess_mean = squeeze(mean(aMat,1,'omitnan'));
    filt_inst_freq(ii).sess_sd = squeeze(std(aMat,1,'omitnan'));
end 

%{
mean over sites  
site_cell = struct2cell(filt_inst_freq);
site_cell = squeeze(site_cell(9,:,:));
aMat = cell2matnan(site_cell',1);
summary.mean = mean(aMat,1,'omitnan');
summary.SEM = std(aMat,1,'omitnan')/sqrt(length(filt_inst_freq));
summary.std = std(aMat,1,'omitnan');
%}
end