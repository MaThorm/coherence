function [grand_struct,angles,inst_freq] = pre_processing_pip_hilb(in_trials,out_trials,V4_trials,medianfiltord)
%PRE_PROCESSING_PIP_HILB Gives back three structures: the Nx1 Hilbert angle
%structure with the result of the hilbert transformation using fieldtrip
%for each session (N). Same for the inst. freq structure which includes the
%derivative of the hilbert transformation. Finally the grand struct
%includes the median filtered inst. freq as well as some summary
%descriptives. 

% Input should be
% Result of pre_processing_pip_trials for 
% Attended 
% Unattended
% V4
% medianfiltord: order of the median filter applied

% Hilberting, Deriving, median filtering, inst change taking, taking mean 
% over trials and recording sites 

parfor ii = 1:length(in_trials)
    [in_hilbertData(ii), in_diffHilbertData(ii), in_medfiltHilbert(ii),in_instchange(ii)] = do_hilbert(in_trials(ii),medianfiltord); 
end 

parfor ii = 1:length(out_trials)
    [out_hilbertData(ii), out_diffHilbertData(ii), out_medfiltHilbert(ii),out_instchange(ii)] = do_hilbert(out_trials(ii),medianfiltord); 
end 

parfor ii = 1:length(V4_trials)
    [V4_hilbertData(ii), V4_diffHilbertData(ii), V4_medfiltHilbert(ii), V4_instchange(ii)] = do_hilbert(V4_trials(ii),medianfiltord); 
end 

% mean over trials in 
for ii = 1:length(in_medfiltHilbert)
    aMat = cell2matnan(in_medfiltHilbert(ii).trial);
    in_medfiltHilbert(ii).sess_mean = squeeze(mean(aMat,1,'omitnan'));
    in_medfiltHilbert(ii).sess_sd = squeeze(std(aMat,1,'omitnan'));
end 
% mean over trials out
for ii = 1:length(out_medfiltHilbert)
    aMat = cell2matnan(out_medfiltHilbert(ii).trial);
    out_medfiltHilbert(ii).sess_mean = squeeze(mean(aMat,1,'omitnan'));
    out_medfiltHilbert(ii).sess_sd = squeeze(std(aMat,1,'omitnan'));
end 

% mean over trials V4
for ii = 1:length(V4_medfiltHilbert)
    aMat = cell2matnan(V4_medfiltHilbert(ii).trial);
    V4_medfiltHilbert(ii).sess_mean = squeeze(mean(aMat,1,'omitnan'));
    V4_medfiltHilbert(ii).sess_sd = squeeze(std(aMat,1,'omitnan'));
end 

grand_struct.out_medfiltHilbert = out_medfiltHilbert;
grand_struct.in_medfiltHilbert = in_medfiltHilbert;
grand_struct.V4_medfiltHilbert = V4_medfiltHilbert; 

% Storing processing steps in seperate files
angles.in_hilbert = in_hilbertData;
angles.out_hilbert = out_hilbertData;
angles.V4_hilbert = V4_hilbertData;
inst_freq.in_inst = in_diffHilbertData; % named it wrong and too lazy to change
inst_freq.out_inst = out_diffHilbertData; % should be inst_freq
inst_freq.V4_inst = V4_diffHilbertData;


% mean over sites in 
site_cell = struct2cell(out_medfiltHilbert);
site_cell = squeeze(site_cell(9,:,:));
aMat = cell2matnan(site_cell');
grand_struct.outsummary.mean = mean(aMat,1,'omitnan');
grand_struct.outsummary.SEM = std(aMat,1,'omitnan')/sqrt(length(out_medfiltHilbert));
grand_struct.outsummary.std = std(aMat,1,'omitnan');

% mean over sites out 
site_cell = struct2cell(in_medfiltHilbert);
site_cell = squeeze(site_cell(9,:,:));
aMat = cell2matnan(site_cell');
grand_struct.insummary.mean = mean(aMat,1,'omitnan');
grand_struct.insummary.SEM = std(aMat,1,'omitnan')/sqrt(length(in_medfiltHilbert));
grand_struct.insummary.std = std(aMat,1,'omitnan');

% mean over sites V4
site_cell = struct2cell(V4_medfiltHilbert);
site_cell = squeeze(site_cell(9,:,:));
aMat = cell2matnan(site_cell');
grand_struct.V4summary.mean = mean(aMat,1,'omitnan');
grand_struct.V4summary.SEM = std(aMat,1,'omitnan')/sqrt(length(V4_medfiltHilbert));
grand_struct.V4summary.std = std(aMat,1,'omitnan');
end

