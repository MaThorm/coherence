clc 
clearvars -except med_i medianorders bpwidth
load('attout_dataset.mat')
load('attin_dataset.mat')
load("V4_dataset.mat")
matpath = '/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files'
medianorders = [20];
bpwidth = [55 95];
toi = [0 4.3];
store_hilbert = true;
%% Creating trials and bp filtering 
% Attin trials
parfor ii = 1:length(attout_dataset)
    in_trials(ii) = do_trialselection(attin_dataset(ii).path,attin_dataset(ii).file,[attin_dataset(ii).chan attin_dataset(ii).V4chan{1}],attin_dataset(ii).stimno,bpwidth,toi);
end 
% Attout trials
parfor ii = 1:length(attout_dataset)
    out_trials(ii) = do_trialselection(attout_dataset(ii).path,attout_dataset(ii).file,[attout_dataset(ii).chan attout_dataset(ii).V4chan{1}],attout_dataset(ii).stimno,bpwidth,toi);
end 
% AttV4 trials
parfor ii = 1:length(V4_dataset)
    V4_trials(ii) = do_trialselection(V4_dataset(ii).path,V4_dataset(ii).file,V4_dataset(ii).chan,V4_dataset(ii).stimno,bpwidth,toi);
end 
%
%%
save(fullfile(matpath,sprintf('attin_trials%d%dtoi%f-%f.mat',bpwidth(1),bpwidth(2),toi(1), toi(2))),'in_trials');
save(fullfile(matpath,sprintf('attout_trials%d%dtoi%f-%f.mat',bpwidth(1),bpwidth(2),toi(1), toi(2))),'out_trials');
save(fullfile(matpath,sprintf('V4_trials%d%dtoi%f-%f.mat.mat',bpwidth(1),bpwidth(2),toi(1), toi(2))),'V4_trials');

%% Hilberting, Deriving, median filtering, inst change taking, taking mean 
% over trials and recording sites 

for med_i = 1:length(medianorders)
clearvars -except med_i medianorders bpwidth store_hilbert toi
matpath = '/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files'

load(fullfile(matpath,sprintf('attin_trials%d%dtoi%f-%f.mat',bpwidth(1),bpwidth(2),toi(1), toi(2))));
load(fullfile(matpath,sprintf('attout_trials%d%dtoi%f-%f.mat',bpwidth(1),bpwidth(2),toi(1), toi(2))));
load(fullfile(matpath,sprintf('V4_trials%d%dtoi%f-%f.mat.mat',bpwidth(1),bpwidth(2),toi(1), toi(2))));


medianfiltord = medianorders(med_i);
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

if store_hilbert == true
    grand_struct.add.in_hilbert = in_hilbertData;
    grand_struct.add.out_hilbert = out_hilbertData;
    grand_struct.add.V4_hilbert = V4_hilbertData;
    %grand_struct.add.inst_freq = in_diffHilbertData;
    %grand_struct.add.outst_freq = out_diffHilbertData;
    %grand_struct.add.V4st_freq = V4_diffHilbertData;
end 

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

% saving
save(fullfile(matpath,sprintf('grand_structbp%d%dmed%dtoi%f-%f.mat',bpwidth(1),bpwidth(2),medianfiltord, toi(1), toi(2))),'grand_struct')
end


%% Loading and combining different structures
clc 
clearvars -except med_i medianorders bpwidth
s_in = struct();
path = '/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files'
for ii = 1:length(medianorders)
    s(ii) = load(fullfile(path,sprintf('grand_structbp45100med%d',medianorders(ii))));
    s_in(ii).input = s(ii).grand_struct.in_medfiltHilbert;
    s_out(ii).input = s(ii).grand_struct.out_medfiltHilbert;
    s_V4(ii).input = s(ii).grand_struct.V4_medfiltHilbert;
end 
[in_cell, in_sdcell] = comb_avg(s_in);
[out_cell, out_sdcell] = comb_avg(s_out);
[V4_cell, V4_sdcell] = comb_avg(s_V4);
summary_struct.in_sites = in_cell;
summary_struct.out_sites = out_cell;
summary_struct.V4_sites = V4_cell;
summary_struct.in_sdsites = in_sdcell;
summary_struct.out_sdsites = out_sdcell;
summary_struct.V4_sdsites = V4_sdcell;

summary_struct.in_mean = mean(in_cell,1,'omitnan');
summary_struct.in_sdmean = mean(in_sdcell,1,'omitnan');
summary_struct.out_mean = mean(out_cell,1,'omitnan');
summary_struct.out_sdmean = mean(out_sdcell,1,'omitnan');
summary_struct.V4_mean = mean(V4_cell,1,'omitnan');
summary_struct.V4_sdmean = mean(V4_sdcell,1,'omitnan');

save(fullfile(path,sprintf('summary_structbp%d%dmed%d_%d_%d_%d_%d.mat',bpwidth(1),bpwidth(2),medianorders(1),medianorders(2),medianorders(3),medianorders(4),medianorders(5))),'summary_struct')