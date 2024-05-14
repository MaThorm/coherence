clc
clear all

% Trial preprocessing parameters
bpfilt = false;
bpwidth = [1 150];
toi = [2.3 4.3]; % Time region of interest [2.3 4.3] is MC 2&3 full trial with - 100 + 100 lol 


% SSD parameters
fs = 1000; % Sampling frequency of the signal
th = 0.01; % residual variance threshhold   
lower = 52;
upper = 90;

matpath = '/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files';


%% Loading structs 
LogicalStr = {'false', 'true'};
load(fullfile(matpath,'trials',sprintf('in_trials_bpfilt%s_toi%.1f-%.1f.mat',LogicalStr{bpfilt+1},toi(1),toi(2))))
load(fullfile(matpath,'trials',sprintf('out_trials_bpfilt%s_toi%.1f-%.1f.mat',LogicalStr{bpfilt+1},toi(1),toi(2))))
load(fullfile(matpath,'trials',sprintf('V4_trials_bpfilt%s_toi%.1f-%.1f.mat',LogicalStr{bpfilt+1},toi(1),toi(2))))
in_SSD_trials = load(fullfile(matpath,'SSD',sprintf('in_ssd_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))));
in_SSD_trials = in_SSD_trials.in_trials;
out_SSD_trials = load(fullfile(matpath,'SSD',sprintf('out_ssd_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))));
out_SSD_trials = out_SSD_trials.out_trials;
V4_SSD_trials = load(fullfile(matpath,'SSD',sprintf('V4_ssd_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))));
V4_SSD_trials = V4_SSD_trials.V4_trials;

%%
% calculating wavelet
clear in_TFRwave out_TFRwave V4_TFRwave
cur_in = in_SSD_trials; % Needs to be swapped out in case SSD should be calculated
cur_out = out_SSD_trials;
cur_V4 = V4_SSD_trials; 
cfg = [];
cfg.width = 7;
cfg.method = 'wavelet'
cfg.output = 'pow';
cfg.foi = exp(linspace(log(5),log(160),35));
cfg.toi = -1.3:0.001:4;
cfg.keeptrials = 'no';
for ii = 1:length(cur_in)
    cfg.channel = cur_in(ii).label{1,1};
    in_TFRwave(ii) = ft_freqanalysis(cfg,cur_in(ii));
end 
for ii = 1:length(cur_out)
    cfg.channel = cur_out(ii).label{1,1};
    out_TFRwave(ii) = ft_freqanalysis(cfg,cur_out(ii));
end 
for ii = 1:length(cur_V4)
    cfg.channel = cur_V4(ii).label{1,1};
    V4_TFRwave(ii) = ft_freqanalysis(cfg,cur_V4(ii));
end 
% Performing baseline correction 
cfg = [];
cfg.baseline = [-1.3 -0.5];
cfg.baselinetype = 'relchange';
clear in_TFRwave_bl out_TFRwave_bl V4_TFRwave_bl
for ii = 1:length(in_TFRwave)
    in_TFRwave_bl(ii) = ft_freqbaseline(cfg,in_TFRwave(ii));
end 
for ii = 1:length(out_TFRwave)
    out_TFRwave_bl(ii) = ft_freqbaseline(cfg,out_TFRwave(ii));
end 
for ii = 1:length(V4_TFRwave)
    V4_TFRwave_bl(ii) = ft_freqbaseline(cfg,V4_TFRwave(ii));
end 

%% Saving structs non SSD 
mt = [];
filename = 'wavelet_red.mat'
savepath = "/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files/wavelet/";
save(fullfile(savepath,filename),'mt','-v7.3')
m = matfile(fullfile(savepath,filename),'Writable',true)
m.in = in_TFRwave_bl;
m.out = out_TFRwave_bl;
m.V4 = V4_TFRwave_bl; 

%% Saving the structs SSD
mt = [];
filename = sprintf('SSDwavelet%d-%d_red.mat',lower,upper);
savepath = "/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files/wavelet/"
save(fullfile(savepath,filename),'mt','-v7.3')
m = matfile(fullfile(savepath,filename),'Writable',true)
m.in = in_TFRwave;
m.out = out_TFRwave;
m.V4 = V4_TFRwave;



%% Plotting an example trial & session 
%sess1mean = mean(in_TFRwave_bl(1).powspctrm,1,'omitnan')
trial = 3
sess = 1;
cur_cond = in_TFRwave;
%for trial = 1:length(cur_cond(sess).trialinfo)
    cur_mean = cur_cond(sess).powspctrm(trial,:,:,:);
    cur_mean = reshape(cur_mean,1,length(cur_cond(sess).freq),length(cur_mean));
    test_struct = cur_cond(sess);
    test_struct.powspctrm = cur_mean;
    cfg = [];
    cfg.maskstyle = ['saturation'];
    cfg.colorbar = 'yes';
    f = figure;
    f.Units = 'normalized';
    f.Position = [0 0 0.5 0.5];
    ft_singleplotTFR(cfg,test_struct)
    xlim([1 2])
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
    yticks = get(gca, 'YTick');
    a(trial) = sum(sum(isnan(cur_mean)));
    f = figure;
    cur_mean_r = reshape(cur_mean,size(cur_mean,2),size(cur_mean,3))
    ft_plot_matrix(in_TFRwave(1).time,in_TFRwave(1).freq,cur_mean_r,'tag', 'cip')
    axis xy
    xlim([1 2])
    colorbar()
%end 
%histogram(a,10);


% ~~~~~~~~~~~~~~ Analysis focused on averages over trials ~~~~~~~~~~~%%%
% Only works when cfg.keeptrials is set to no
%% Plotting avg. of sessions 
figpath = "/data/projects/V1V4coherence/03_results_max/wavelets/SSD/th52-90/";
cur_var = out_TFRwave;
var_name = getVarName(out_TFRwave);
a = plot_avg(figpath,cur_var,var_name)



%% Plotting grand average In
out = in_TFRwave_bl(1);
lfax = log(out.freq');
nYtick = 17;
dyt = (lfax(end)-lfax(1))/(nYtick-1);
lyt = (0:nYtick-1)*dyt+lfax(1);
yt = exp(lyt);
lfax2 =  1:35;
dyt2   = (lfax2(end)-lfax2(1))/(nYtick-1);
lyt2   = (0:nYtick-1)*dyt2+lfax2(1);
yt2    = (lyt2);
figpath = "/data/projects/V1V4coherence/03_results_max/wavelets/blcor/";

cur_var = in_TFRwave_bl;
mean_array = zeros(1,size(cur_var(1).powspctrm,2),size(cur_var(1).powspctrm,3))
for ii = 1:length(cur_var);
    mean_array = mean_array + cur_var(ii).powspctrm;
end 
mean_array = mean_array / length(cur_var)
f = figure;
f.Units = 'normalized';
f.Position = [0 0 0.5 0.5];
imagesc(out.time,1:35,squeeze(mean_array))
title('Attended Wavelet Summary')
ax = gca;
set(ax,...
'color','none',...
'TickDir','out',...
'YDir','normal',...
'YTick',yt2,...
'YTickLabel',num2str(yt','%3.2f'));
colormap(jet(256));
colorbar();
xlabel('Time [s]')
ylabel('Frequency [Hz]')
clim([-0.5 2.5])
saveas(gcf,fullfile(figpath,'AttIn_Wavelet_Summary.fig'));
saveas(gcf,fullfile(figpath,'AttIn_Wavelet_Summary.jpg'));

% Saving out
cur_var = out_TFRwave_bl;
mean_array = zeros(1,size(cur_var(1).powspctrm,2),size(cur_var(1).powspctrm,3))
for ii = 1:length(cur_var);
    mean_array = mean_array + cur_var(ii).powspctrm;
end 
mean_array = mean_array / length(cur_var)
cur_var(1).powspctrm = mean_array

f = figure;
f.Units = 'normalized';
f.Position = [0 0 0.5 0.5];
imagesc(out.time,1:35,squeeze(mean_array))
title('Non-Attended Wavelet Summary')
ax = gca;
set(ax,...
'color','none',...
'TickDir','out',...
'YDir','normal',...
'YTick',yt2,...
'YTickLabel',num2str(yt','%3.2f'));
colormap(jet(256));
colorbar();
xlabel('Time [s]')
ylabel('Frequency [Hz]')
clim([-0.5 2.5])
saveas(gcf,fullfile(figpath,'AttOut_Wavelet_Summary.fig'));
saveas(gcf,fullfile(figpath,'AttOut_Wavelet_Summary.jpg'));

% Saving V4
cur_var = V4_TFRwave_bl;
mean_array = zeros(1,size(cur_var(1).powspctrm,2),size(cur_var(1).powspctrm,3))
for ii = 1:length(cur_var);
    mean_array = mean_array + cur_var(ii).powspctrm;
end 
mean_array = mean_array / length(cur_var)
cur_var(1).powspctrm = mean_array
f = figure;
f.Units = 'normalized';
f.Position = [0 0 0.5 0.5];
imagesc(out.time,1:35,squeeze(mean_array))
colorbar()
title('V4 Wavelet Summary')
ax = gca;
set(ax,...
'color','none',...
'TickDir','out',...
'YDir','normal',...
'YTick',yt2,...
'YTickLabel',num2str(yt','%3.2f'));
xlabel('Time [s]')
ylabel('Frequency [Hz]')
saveas(gcf,fullfile(figpath,'V4_Wavelet_Summary.fig'));
saveas(gcf,fullfile(figpath,'V4_Wavelet_Summary.jpg'));

%%
function [a] = plot_avg(figpath,cur_var,var_name)
    cfg = [];
    cfg.maskstyle = ['saturation'];
    cfg.colorbar = 'yes';
    for i_s = 1:length(cur_var)
        out = cur_var(i_s);
        lfax = log(out.freq');
        nYtick = 17;
        dyt = (lfax(end)-lfax(1))/(nYtick-1);
        lyt = (0:nYtick-1)*dyt+lfax(1);
        yt = exp(lyt);
        lfax2 =  1:35;
        dyt2   = (lfax2(end)-lfax2(1))/(nYtick-1);
        lyt2   = (0:nYtick-1)*dyt2+lfax2(1);
        yt2    = (lyt2);
        f = figure;
        f.Units = 'normalized';
        f.Position = [0 0 0.5 0.5];
        imagesc(out.time(2300:4300),1:35,squeeze(out.powspctrm(:,:,2300:4300)))
        title(sprintf("%s session %i",var_name,i_s))
        colormap(jet(256));
        colorbar();
        ax = gca;
        set(ax,...
        'color','none',...
        'TickDir','out',...
        'YDir','normal',...
        'YTick',yt2,...
        'YTickLabel',num2str(yt','%3.2f'));
        xlabel('Time [s]')
        ylabel('Frequency [Hz]')  
        %xlim([1 3])
        %clim([-0.5 6])
        foldername = fullfile(figpath,sprintf('session_%d',i_s));
        if ~exist(foldername,'dir')
            mkdir(foldername)
        end 
        saveas(gcf,fullfile(foldername,sprintf('%s_sess%d.fig',var_name,i_s)));
        saveas(gcf,fullfile(foldername,sprintf('%s_sess%d.jpg',var_name,i_s)));
        close all 
        a = 1;
    end 
end 


