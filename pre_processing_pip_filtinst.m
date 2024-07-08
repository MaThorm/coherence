function [phase_angles,inst_freq,filt_data] = pre_processing_pip_hilb(phase_angles_wr,filttype,framelen,filtord)
%PRE_PROCESSING_PIP_HILB Feeds fieldtrip structs to do_filtinst function
%which filters the data and takes the instantaneous frequency. 
% Input should be
% angles: a fieldtrip struct containing wrapped phase angles 
% filttype: either 'medfilt' or 'sgolay'
% framelen: framelength of the filter
% filtord: order of the filter applied (only for sgolay)

for ii = 1:length(phase_angles_wr)
    [phase_angles(ii), inst_freq(ii), filt_data(ii)] = do_filtinst(phase_angles_wr(ii),filttype,framelen,filtord); 
end 

% % mean over trials 
% for ii = 1:length(filt_data)
%     aMat = cell2matnan(filt_data(ii).trial,1);
%     filt_data(ii).sess_mean = squeeze(mean(aMat,1,'omitnan'));
%     filt_data(ii).sess_sd = squeeze(std(aMat,1,'omitnan'));
% end 


end