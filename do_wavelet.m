clc 
clear all
load('attout_dataset.mat')
load('attin_dataset.mat')
load("V4_dataset.mat")
% Creating trials and bp filtering 
% Attin trials
bpwidth = [1 150];

parfor ii = 1:length(attin_dataset)
    in_trials(ii) = do_trialselection(attin_dataset(ii).path,attin_dataset(ii).file,attin_dataset(ii).chan,attin_dataset(ii).stimno,bpwidth);
end 
parfor ii = 1:length(attout_dataset)
    out_trials(ii) = do_trialselection(attout_dataset(ii).path,attout_dataset(ii).file,attout_dataset(ii).chan,attout_dataset(ii).stimno,bpwidth);
end 
parfor ii = 1:length(V4_dataset)
    V4_trials(ii) = do_trialselection(V4_dataset(ii).path,V4_dataset(ii).file,V4_dataset(ii).chan,V4_dataset(ii).stimno,bpwidth);
end 

%%
% calculating wavelet
cfg = [];
cfg.width = 7;
cfg.method = 'wavelet'
cfg.output = 'pow';
cfg.foi = 1:4:100;
cfg.toi = -1.3:0.001:4;
for ii = 1:length(in_trials)
    in_TFRwave(ii) = ft_freqanalysis(cfg,in_trials(ii));
end 
for ii = 1:length(out_trials)
    out_TFRwave(ii) = ft_freqanalysis(cfg,out_trials(ii));
end 
for ii = 1:length(V4_trials)
    V4_TFRwave(ii) = ft_freqanalysis(cfg,V4_trials(ii));
end 
cfg = [];
cfg.baseline = [-1.3 -0.5];
cfg.baselinetype = 'relchange';
for ii = 1:length(in_TFRwave)
    in_TFRwave_bl(ii) = ft_freqbaseline(cfg,in_TFRwave(ii));
end 
for ii = 1:length(out_TFRwave)
    out_TFRwave_bl(ii) = ft_freqbaseline(cfg,out_TFRwave(ii));
end 
for ii = 1:length(V4_TFRwave)
    V4_TFRwave_bl(ii) = ft_freqbaseline(cfg,V4_TFRwave(ii));
end 

attin.site_cell_bl = struct2cell(in_TFRwave_bl);
attin.site_cell_bl = squeeze(attin.site_cell_bl(5,:,:));
attin.site_mat_bl = cell2mat(attin.site_cell_bl);
attin.site_mat_bl = squeeze(mean(attin.site_mat_bl,1,'omitnan'));
attout.site_cell_bl = struct2cell(out_TFRwave_bl);
attout.site_cell_bl = squeeze(attout.site_cell_bl(5,:,:));
attout.site_mat_bl = cell2mat(attout.site_cell_bl);
attout.site_mat_bl = squeeze(mean(attout.site_mat_bl,1,'omitnan'));
V4.site_cell_bl = struct2cell(V4_TFRwave_bl);
V4.site_cell_bl = squeeze(V4.site_cell_bl(5,:,:));
V4.site_mat_bl = cell2mat(V4.site_cell_bl);
V4.site_mat_bl = squeeze(mean(V4.site_mat_bl,1,'omitnan'));
%
attin.freq_site_bl = attin.site_mat_bl(:,2300:4300);
attin.freq_site_bl = mean(attin.freq_site_bl,2,'omitnan');
attout.freq_site_bl = attout.site_mat_bl(:,2300:4300);
attout.freq_site_bl = mean(attout.freq_site_bl,2,'omitnan');
V4.freq_site_bl = V4.site_mat_bl(:,2300:4300);
V4.freq_site_bl = mean(V4.freq_site_bl,2,'omitnan');

%
attin.site_cell = struct2cell(in_TFRwave);
attin.site_cell = squeeze(attin.site_cell(5,:,:));
attin.site_mat = cell2mat(attin.site_cell);
attin.site_mat = squeeze(mean(attin.site_mat,1,'omitnan'));
attout.site_cell = struct2cell(out_TFRwave);
attout.site_cell = squeeze(attout.site_cell(5,:,:));
attout.site_mat = cell2mat(attout.site_cell);
attout.site_mat = squeeze(mean(attout.site_mat,1,'omitnan'));
V4.site_cell = struct2cell(V4_TFRwave);
V4.site_cell = squeeze(V4.site_cell(5,:,:));
V4.site_mat = cell2mat(V4.site_cell);
V4.site_mat = squeeze(mean(V4.site_mat,1,'omitnan'));
%
attin.freq_site = attin.site_mat(:,2300:4300);
attin.freq_site = mean(attin.freq_site,2,'omitnan');
attout.freq_site = attout.site_mat(:,2300:4300);
attout.freq_site = mean(attout.freq_site,2,'omitnan');
V4.freq_site = V4.site_mat(:,2300:4300);
V4.freq_site = mean(V4.freq_site,2,'omitnan');

%% Saving
matpath = "/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files"
wavelet_struct.in.trials = in_trials;
wavelet_struct.out.trials = out_trials;
wavelet_struct.V4.trials = V4_trials; 
wavelet_struct.in.TFR = in_TFRwave;
wavelet_struct.out.TFR = out_TFRwave;
wavelet_struct.V4.TFR = V4_TFRwave;
wavelet_struct.in.TFRbl = in_TFRwave_bl;
wavelet_struct.out.TFRbl = out_TFRwave_bl;
wavelet_struct.V4.TFRbl = V4_TFRwave_bl;
wavelet_struct.in.freqsite = attin.freq_site;
wavelet_struct.out.freqsite = attout.freq_site;
wavelet_struct.V4.freqsite = V4.freq_site;
wavelet_struct.in.site_mat = attin.site_mat;
wavelet_struct.out.site_mat = attout.site_mat;
wavelet_struct.V4.site_mat = V4.site_mat;
wavelet_struct.in.freqsite_bl = attin.freq_site_bl;
wavelet_struct.out.freqsite_bl = attout.freq_site_bl;
wavelet_struct.V4.freqsite_bl = V4.freq_site_bl;
wavelet_struct.in.site_mat_bl = attin.site_mat_bl;
wavelet_struct.out.site_mat_bl = attout.site_mat_bl;
wavelet_struct.V4.site_mat_bl = V4.site_mat_bl;
save(fullfile(matpath,'wavelet_struct.mat'),"wavelet_struct");
%%
cfg = [];
cfg.maskstyle = ['saturation'];
cfg.colorbar = 'yes';
figure
ft_singleplotTFR(cfg,wavelet_struct.V4.TFR(6))
clim([0 5e-8])

%% 
figure 
imagesc(flip(wavelet_struct.in.site_mat))
colorbar()
ax = gca;
ax.YTick = 1:25;
ax.YTickLabels = round(flip(wavelet_struct.in.TFRbl(1).freq),2);

%% plotting bl correct
semilogx(wavelet_struct.in.TFRbl(1).freq,wavelet_struct.in.freqsite_bl)
hold on 
semilogx(wavelet_struct.out.TFRbl(1).freq,wavelet_struct.out.freqsite_bl)
semilogx(wavelet_struct.V4.TFRbl(1).freq,wavelet_struct.V4.freqsite_bl)
legend('AttIn','AttOut','V4')

%% plotting nonbl correct
semilogx(wavelet_struct.in.TFRbl(1).freq,wavelet_struct.in.freqsite)
hold on  
semilogx(wavelet_struct.out.TFRbl(1).freq,wavelet_struct.out.freqsite)
hold off
legend('AttIn','AttOut')
%% 
semilogx(wavelet_struct.V4.TFRbl(1).freq,wavelet_struct.V4.freqsite)