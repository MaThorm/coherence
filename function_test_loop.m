clc 
clear all
load('attout_dataset.mat')
load('attin_dataset.mat')
load("V4_dataset.mat")
matpath = '/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files';
blow1 = 30;
blow2 = 55;
bhigh1 = 70;
bhigh2 = 100;
bstep = 5;
for b1= blow1:bstep:blow2
    for b2 = bhigh1:bstep:bhigh2
    bpwidth = [b1 b2]   
for ii = 1:length(attout_dataset)
    in_trials(ii) = do_trialselection(attin_dataset(ii).path,attin_dataset(ii).file,attin_dataset(ii).chan,attin_dataset(ii).stimno,bpwidth);
end 

% Attout trials
for ii = 1:length(attout_dataset)
    out_trials(ii) = do_trialselection(attout_dataset(ii).path,attout_dataset(ii).file,attout_dataset(ii).chan,attout_dataset(ii).stimno,bpwidth);
end 

% AttV4 trials
for ii = 1:length(V4_dataset)
    V4_trials(ii) = do_trialselection(V4_dataset(ii).path,V4_dataset(ii).file,V4_dataset(ii).chan,V4_dataset(ii).stimno,bpwidth);
end 
%

save(fullfile(matpath,sprintf('attin%d%d.trials',bpwidth(1),bpwidth(2))),'in_trials');
save(fullfile(matpath,sprintf('attout%d%d.trials',bpwidth(1),bpwidth(2))),'out_trials');
save(fullfile(matpath,sprintf('V4%d%d.trials',bpwidth(1),bpwidth(2))),'V4_trials');

% Hilberting, Deriving, median filtering, inst change taking, taking mean 
% over trials and recording sites 
medianfiltord = 30;
for ii = 1:length(in_trials)
    [in_hilbertData(ii), in_diffHilbertData(ii), in_medfiltHilbert(ii),in_instchange(ii)] = do_hilbert(in_trials(ii),medianfiltord); 
end 

for ii = 1:length(out_trials)
    [out_hilbertData(ii), out_diffHilbertData(ii), out_medfiltHilbert(ii),out_instchange(ii)] = do_hilbert(out_trials(ii),medianfiltord); 
end 

for ii = 1:length(V4_trials)
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

% Mean over trials inst change
% Change: should put this as a function 
for ii = 1:length(in_instchange)
    aMat = cell2matnan(in_instchange(ii).diff);
    in_instchange(ii).sess_mean = squeeze(mean(aMat,1,'omitnan'));
    in_instchange(ii).sess_sd = squeeze(std(aMat,1,'omitnan'));
end 

% mean over trials out
for ii = 1:length(out_instchange)
    aMat = cell2matnan(out_instchange(ii).diff);
    out_instchange(ii).sess_mean = squeeze(mean(aMat,1,'omitnan'));
    out_instchange(ii).sess_sd = squeeze(std(aMat,1,'omitnan'));
end 

% mean over trials V4
for ii = 1:length(V4_instchange)
    aMat = cell2matnan(V4_instchange(ii).diff);
    V4_instchange(ii).sess_mean = squeeze(mean(aMat,1,'omitnan'));
    V4_instchange(ii).sess_sd = squeeze(std(aMat,1,'omitnan'));
end 

grand_struct.in_instchange = in_instchange;
grand_struct.out_instchange = out_instchange;
grand_struct.V4_instchange = V4_instchange;
%
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

% mean over sites inst change
% mean over sites in 
site_cell = struct2cell(out_instchange);
site_cell = squeeze(site_cell(2,:,:));
aMat = cell2matnan(site_cell');
grand_struct.outsummary.change_mean = mean(aMat,1,'omitnan');
grand_struct.outsummary.change_SEM = std(aMat,1,'omitnan')/sqrt(length(out_instchange));
grand_struct.outsummary.change_std = std(aMat,1,'omitnan');

% mean over sites out 
site_cell = struct2cell(in_instchange);
site_cell = squeeze(site_cell(2,:,:));
aMat = cell2matnan(site_cell');
grand_struct.insummary.change_mean = mean(aMat,1,'omitnan');
grand_struct.insummary.change_SEM = std(aMat,1,'omitnan')/sqrt(length(in_instchange));
grand_struct.insummary.change_std = std(aMat,1,'omitnan');

% mean over sites V4
site_cell = struct2cell(V4_instchange);
site_cell = squeeze(site_cell(2,:,:));
aMat = cell2matnan(site_cell');
grand_struct.V4summary.change_mean = mean(aMat,1,'omitnan');
grand_struct.V4summary.change_SEM = std(aMat,1,'omitnan')/sqrt(length(V4_instchange));
grand_struct.V4summary.change_std = std(aMat,1,'omitnan');

% saving
save(fullfile(matpath,sprintf('grand_structbp%d%dmed%d.mat',bpwidth(1),bpwidth(2),medianfiltord)),'grand_struct')
clearvars -except blow1 blow2 bhigh1 bhigh2 b1 b2 attout_dataset attin_dataset V4_dataset matpath
    end
end 