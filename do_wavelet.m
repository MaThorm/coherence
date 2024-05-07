clc
clear all

% Trial preprocessing parameters
bpfilt = false;
bpwidth = [1 150];
toi = [-100 100]; % Time region of interest [2.3 4.3] is MC 2&3 full trial with - 100 + 100 lol 


% SSD parameters
fs = 1000; % Sampling frequency of the signal
th = 0.01; % residual variance threshhold   
lower = 52;
upper = 110;

matpath = '/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files';


%% Loading structs 
LogicalStr = {'false', 'true'};
load(fullfile(matpath,'trials',sprintf('in_trials_bpfilt%s_toi%.1f-%.1f.mat',LogicalStr{bpfilt+1},toi(1),toi(2))))
load(fullfile(matpath,'trials',sprintf('out_trials_bpfilt%s_toi%.1f-%.1f.mat',LogicalStr{bpfilt+1},toi(1),toi(2))))
load(fullfile(matpath,'trials',sprintf('V4_trials_bpfilt%s_toi%.1f-%.1f.mat',LogicalStr{bpfilt+1},toi(1),toi(2))))
%in_SSD_trials = load(fullfile(matpath,'SSD',sprintf('in_ssd_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))));
%in_SSD_trials = in_SSD_trials.in_trials;
%out_SSD_trials = load(fullfile(matpath,'SSD',sprintf('out_ssd_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))));
%out_SSD_trials = out_SSD_trials.out_trials;
%V4_SSD_trials = load(fullfile(matpath,'SSD',sprintf('V4_ssd_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))));
%V4_SSD_trials = V4_SSD_trials.V4_trials;

%%
% calculating wavelet
clear in_TFRwave out_TFRwave V4_TFRwave
cur_in = in_trials;
cur_out = out_trials;
cur_V4 = V4_trials; 
cfg = [];
cfg.width = 7;
cfg.method = 'wavelet'
cfg.output = 'pow';
cfg.foi = exp(linspace(log(5),log(160),35));
cfg.toi = -1.3:0.002:4;
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



%% PLotting and saving single trials only MC 2&3
% Attended 
figpath = '/data/projects/V1V4coherence/03_results_max/wavelets/blcor_MC23'
cur_TFR = in_TFRwave_bl;
for i_s = 1:length(cur_TFR)
    for i_t = 1:length(cur_TFR(i_s).trialinfo)
        cur_mean = cur_TFR(i_s).powspctrm(i_t,:,:,:);
        cur_mean = reshape(cur_mean,1,length(cur_TFR(i_s).freq),length(cur_mean));
        test_struct = cur_TFR(i_s);
        test_struct.powspctrm = cur_mean;
        cfg = [];
        cfg.maskstyle = ['saturation'];
        cfg.colorbar = 'yes';
        cfg.title = sprintf('V1 attended Power Spectrum MC2/3: Session %d, Trial %d, Electrode %s',i_s,i_t,cur_TFR(i_s).label{1});
        f = figure;
        f.Units = 'normalized';
        f.Position = [0 0 0.5 0.5];
        ft_singleplotTFR(cfg,test_struct)
        if length(in_trials(i_s).trial{1,i_t}) <= 4000
            xlim([1 2])
        else 
            xlim([1 3])
        end
        xlabel('Time [s]')
        ylabel('Frequency [Hz]')
        set(gcf, 'Visible','off')
        foldername = fullfile(figpath,sprintf('session_%d/V1/attend_in',i_s));
        if ~exist(foldername,'dir')
            mkdir(foldername)
        end 
        saveas(gcf,fullfile(foldername,sprintf('attendin_wavelet_MC23_sess%d_trial%d.fig',i_s,i_t)));
        saveas(gcf,fullfile(foldername,sprintf('attendin_wavelet_MC23_sess%d_trial%d.jpg',i_s,i_t)));
        close all 
    end 
end 

% Attend out
cur_TFR = out_TFRwave_bl;
for i_s = 1:length(cur_TFR)
    for i_t = 1:length(cur_TFR(i_s).trialinfo)
        cur_mean = cur_TFR(i_s).powspctrm(i_t,:,:,:);
        cur_mean = reshape(cur_mean,1,length(cur_TFR(i_s).freq),length(cur_mean));
        test_struct = cur_TFR(i_s);
        test_struct.powspctrm = cur_mean;
        cfg = [];
        cfg.maskstyle = ['saturation'];
        cfg.colorbar = 'yes';  
        cfg.title = sprintf('V1 non-attended Power Spectrum MC 2/3: Session %d, Trial %d, Electrode %s',i_s,i_t,cur_TFR(i_s).label{1});
        f = figure;
        f.Units = 'normalized';
        f.Position = [0 0 0.5 0.5];
        ft_singleplotTFR(cfg,test_struct)
        if length(out_trials(i_s).trial{1,i_t}) <= 4000
            xlim([1 2])
        else 
            xlim([1 3])
        end
        xlabel('Time [s]')
        ylabel('Frequency [Hz]')
        set(gcf, 'Visible','off')
        foldername = fullfile(figpath,sprintf('session_%d/V1/attend_out',i_s));
        if ~exist(foldername,'dir')
            mkdir(foldername)
        end 
        saveas(gcf,fullfile(foldername,sprintf('attendout_wavelet_MC23_sess%d_trial%d.fig',i_s,i_t)));
        saveas(gcf,fullfile(foldername,sprintf('attendout_wavelet_MC23_sess%d_trial%d.jpg',i_s,i_t)));
        close all 
    end 
end 

% V4
cur_TFR = V4_TFRwave_bl;
for i_s = 1:length(cur_TFR)
    for i_t = 1:length(cur_TFR(i_s).trialinfo)
        cur_mean = cur_TFR(i_s).powspctrm(i_t,:,:,:);
        cur_mean = reshape(cur_mean,1,length(cur_TFR(i_s).freq),length(cur_mean));
        test_struct = cur_TFR(i_s);
        test_struct.powspctrm = cur_mean;
        cfg = [];
        cfg.maskstyle = ['saturation'];
        cfg.colorbar = 'yes';  
        cfg.title = sprintf('V1 attended Power Spectrum MC 2/3: Session %d, Trial %d, Electrode %s',i_s,i_t,cur_TFR(i_s).label{1});
        f = figure;
        f.Units = 'normalized';
        f.Position = [0 0 0.5 0.5];
        ft_singleplotTFR(cfg,test_struct)
        if length(V4_trials(i_s).trial{1,i_t}) <= 4000
            xlim([1 2])
        else 
            xlim([1 3])
        end
        xlabel('Time [s]')
        ylabel('Frequency [Hz]')    
        set(gcf, 'Visible','off')
        foldername = fullfile(figpath,sprintf('session_%d/V4',i_s));
        if ~exist(foldername,'dir')
            mkdir(foldername)
        end 
        saveas(gcf,fullfile(foldername,sprintf('V4_wavelet_MC23_sess%d_trial%d.fig',i_s,i_t)));
        saveas(gcf,fullfile(foldername,sprintf('V4_wavelet_MC23_sess%d_trial%d.jpg',i_s,i_t)));
        close all 
    end 
end 
%% Saving structs non SSD 
mt = [];
savepath = "/data/projects/V1V4coherence/03_results_max/wavelets/";
save(fullfile(savepath,'Wavelet.mat'),'mt','-v7.3')
m = matfile(fullfile(savepath,'Wavelet.mat'),'Writable',true)
m.in = in_TFRwave_bl;
m.out = out_TFRwave_bl;
m.V4 = V4_TFRwave_bl;  
%% PLotting and saving SSD components 
% Attended 
figpath = sprintf('/data/projects/V1V4coherence/03_results_max/wavelets/SSD/th%d-%d',lower,upper)
cur_TFR = in_TFRwave;
for i_s = 1:length(cur_TFR)
    for i_t = 1:length(cur_TFR(i_s).trialinfo)
        cur_mean = cur_TFR(i_s).powspctrm(i_t,:,:,:);
        cur_mean = reshape(cur_mean,1,length(cur_TFR(i_s).freq),length(cur_mean));
        test_struct = cur_TFR(i_s);
        test_struct.powspctrm = cur_mean;
        cfg = [];
        cfg.maskstyle = ['saturation'];
        cfg.colorbar = 'yes';
        cfg.title = sprintf('V1 attended Power Spectrum MC2/3: Session %d, Trial %d, Electrode %s',i_s,i_t,cur_TFR(i_s).label{1});
        f = figure;
        f.Units = 'normalized';
        f.Position = [0 0 0.5 0.5];
        ft_singleplotTFR(cfg,test_struct)
        if length(in_trials(i_s).trial{1,i_t}) <= 1500
            xlim([1 2])
        else 
            xlim([1 3])
        end
        xlabel('Time [s]')
        ylabel('Frequency [Hz]')
        set(gcf, 'Visible','off')
        foldername = fullfile(figpath,sprintf('session_%d/V1/attend_in',i_s));
        if ~exist(foldername,'dir')
            mkdir(foldername)
        end 
        saveas(gcf,fullfile(foldername,sprintf('attendin_wavelet_SSD_sess%d_trial%d.fig',i_s,i_t)));
        saveas(gcf,fullfile(foldername,sprintf('attendin_wavelet_SSD_sess%d_trial%d.jpg',i_s,i_t)));
        close all 
    end 
end 

% Attend out
cur_TFR = out_TFRwave;
for i_s = 1:length(cur_TFR)
    for i_t = 1:length(cur_TFR(i_s).trialinfo)
        cur_mean = cur_TFR(i_s).powspctrm(i_t,:,:,:);
        cur_mean = reshape(cur_mean,1,length(cur_TFR(i_s).freq),length(cur_mean));
        test_struct = cur_TFR(i_s);
        test_struct.powspctrm = cur_mean;
        cfg = [];
        cfg.maskstyle = ['saturation'];
        cfg.colorbar = 'yes';  
        cfg.title = sprintf('V1 non-attended Power Spectrum MC 2/3: Session %d, Trial %d, Electrode %s',i_s,i_t,cur_TFR(i_s).label{1});
        f = figure;
        f.Units = 'normalized';
        f.Position = [0 0 0.5 0.5];
        ft_singleplotTFR(cfg,test_struct)
        if length(out_trials(i_s).trial{1,i_t}) <= 1500
            xlim([1 2])
        else 
            xlim([1 3])
        end
        xlabel('Time [s]')
        ylabel('Frequency [Hz]')
        set(gcf, 'Visible','off')
        foldername = fullfile(figpath,sprintf('session_%d/V1/attend_out',i_s));
        if ~exist(foldername,'dir')
            mkdir(foldername)
        end 
        saveas(gcf,fullfile(foldername,sprintf('attendout_wavelet_SSD_sess%d_trial%d.fig',i_s,i_t)));
        saveas(gcf,fullfile(foldername,sprintf('attendout_wavelet_SSD_sess%d_trial%d.jpg',i_s,i_t)));
        close all 
    end 
end 

% V4
cur_TFR = V4_TFRwave;
for i_s = 1:length(cur_TFR)
    for i_t = 1:length(cur_TFR(i_s).trialinfo)
        cur_mean = cur_TFR(i_s).powspctrm(i_t,:,:,:);
        cur_mean = reshape(cur_mean,1,length(cur_TFR(i_s).freq),length(cur_mean));
        test_struct = cur_TFR(i_s);
        test_struct.powspctrm = cur_mean;
        cfg = [];
        cfg.maskstyle = ['saturation'];
        cfg.colorbar = 'yes';  
        cfg.title = sprintf('V1 attended Power Spectrum MC 2/3: Session %d, Trial %d, Electrode %s',i_s,i_t,cur_TFR(i_s).label{1});
        f = figure;
        f.Units = 'normalized';
        f.Position = [0 0 0.5 0.5];
        ft_singleplotTFR(cfg,test_struct)
        if length(V4_trials(i_s).trial{1,i_t}) <= 1500
            xlim([1 2])
        else 
            xlim([1 3])
        end
        xlabel('Time [s]')
        ylabel('Frequency [Hz]')    
        set(gcf, 'Visible','off')
        foldername = fullfile(figpath,sprintf('session_%d/V4',i_s));
        if ~exist(foldername,'dir')
            mkdir(foldername)
        end 
        saveas(gcf,fullfile(foldername,sprintf('V4_wavelet_SSD_sess%d_trial%d.fig',i_s,i_t)));
        saveas(gcf,fullfile(foldername,sprintf('V4_wavelet_SSD_sess%d_trial%d.jpg',i_s,i_t)));
        close all 
    end 
end 

%% Saving the structs SSD
mt = [];
savepath = "/data/projects/V1V4coherence/03_results_max/wavelets/SSD/"
save(fullfile(savepath,sprintf('SSD_wavelet%d-%d.mat',lower,upper)),'mt','-v7.3')
m = matfile(fullfile(savepath,sprintf('SSD_wavelet%d-%d.mat',lower,upper)),'Writable',true)
m.in = in_TFRwave;
m.out = out_TFRwave;
m.V4 = V4_TFRwave;




% ~~~~~~~~~~~~~~ Analysis focused on averages over trials ~~~~~~~~~~~%%%
% Only works when cfg.keeptrials is set to no
%% Plotting avg. of sessions 
figpath = "/data/projects/V1V4coherence/03_results_max/wavelets/blcor_MC23/";
cur_var = V4_TFRwave_bl;
var_name = getVarName(V4_TFRwave_bl);
a = plot_avg(figpath,cur_var,var_name)

%% Plotting grand average In
figpath = "/data/projects/V1V4coherence/03_results_max/wavelets/blcor_MC23/";
cur_var = in_TFRwave_bl;
mean_array = zeros(1,size(cur_var(1).powspctrm,2),size(cur_var(1).powspctrm,3))
for ii = 1:length(cur_var);
    mean_array = mean_array + cur_var(ii).powspctrm;
end 
mean_array = mean_array / length(cur_var)
cur_var(1).powspctrm = mean_array
cfg = [];
cfg.maskstyle = ['saturation'];
cfg.colorbar = 'yes';
f = figure;
f.Units = 'normalized';
f.Position = [0 0 0.5 0.5];
ft_singleplotTFR(cfg,cur_var(1))
xlabel('Time [s]')
xlim([1 3])
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
cfg = [];
cfg.maskstyle = ['saturation'];
cfg.colorbar = 'yes';
f = figure;
f.Units = 'normalized';
f.Position = [0 0 0.5 0.5];
ft_singleplotTFR(cfg,cur_var(1))
xlim([1 3])
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
cfg = [];
cfg.maskstyle = ['saturation'];
cfg.colorbar = 'yes';
f = figure;
f.Units = 'normalized';
f.Position = [0 0 0.5 0.5];
ft_singleplotTFR(cfg,cur_var(1))
xlim([1 3])
xlabel('Time [s]')
ylabel('Frequency [Hz]')
saveas(gcf,fullfile(figpath,'V4_Wavelet_Summary.fig'));
saveas(gcf,fullfile(figpath,'V4_Wavelet_Summary.jpg'));

%% Functions 
function [a] = plot_avg(figpath,cur_var,var_name)
    cfg = [];
    cfg.maskstyle = ['saturation'];
    cfg.colorbar = 'yes';
    for i_s = 1:length(cur_var)
        f = figure;
        f.Units = 'normalized';
        f.Position = [0 0 0.5 0.5];
        ft_singleplotTFR(cfg,cur_var(i_s))
        xlabel('Time [s]')
        ylabel('Frequency [Hz]')  
        xlim([1 3])
        %clim([-0.5 3])
        set(gcf, 'Visible','off')
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


