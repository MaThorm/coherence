% Script to test filter specifications 
clear 
% Trial selection 
% Loading
load("/data/projects/V1V4coherence/02_analysis_max/git_repos/params.mat")
load('attin_dataset.mat')

[in_trials] = pre_processing_pip_trials(attin_dataset,params.bpfilt,params.bpwidth,params.order)


%%
Fs = 1000;
T = 2;
f1 = 55;
A = 10; 
t = 0:1/Fs:T-1/Fs;
trials = in_trials(1);
% generating a fake signal 
noise_mods = 1:0.5:50;
for ii_n = 1:length(noise_mods)
noise_mod = noise_mods(ii_n);
Fs = 1000;
T = 2;
f1 = 55;
% f2 = 65;
A = 10; 
t = 0:1/Fs:T-1/Fs;
signal_a = A * sin(2*pi*f1*t);
% signal_b = A * sin(2*pi*f2*t);
signal_b = 0;
signal = signal_a + signal_b;
noise = randn(size(t))*noise_mod;
n_sig = signal + noise;
% Exchangin things 
for ii = 1:length(trials.trial)
    trials.trial{1,ii}= n_sig;
    trials.trial{1,ii}(2,:)= n_sig;
    trials.time{1,ii} = t;
end 


%
% filtering the data
filt_trials = trials;
Fbp = [50 100];
order = 1000;
type = 'fir';
dir = 'twopass';
for i_t = 1:length(trials.trial)
    for i_c = 1:size(trials.trial{1,i_t},1)
        cur_trial = trials.trial{1,i_t}(i_c,:);
        filt_trial = ft_preproc_bandpassfilter(cur_trial,Fs,Fbp,order,type,dir);
        filt_trials.trial{1,i_t}(i_c,:) = filt_trial;
    end 
end 
hilb_trial = pre_processing_pip_hilb(filt_trials);
[hilb_angles,inst_freq,filt_data]  = pre_processing_pip_filtinst(hilb_trial,params.filttype,params.framelen,params.filtord,params.toi);
filt_data.label{1,1} = 'V1';
cfg = [];
cfg.channel = 'V1';
s1avg = ft_timelockanalysis(cfg,filt_data);

%
% f = figure;
% plot(s1avg.avg);
b(ii_n) = mean(s1avg.avg);
end 
%%
plot(b)
%% Pip with mannual filter 
clc
clear 
% Trial selection 
% Loading
load("/data/projects/V1V4coherence/02_analysis_max/git_repos/params.mat")
load('attin_dataset.mat')
load('V4_dataset.mat')

% Trialselection 
[in_trials] = pre_processing_pip_trials(attin_dataset,params.bpfilt,params.bpwidth,params.order);
in_trials_nbp = pre_processing_pip_trials(attin_dataset,false,params.bpwidth,params.order);
[V4_trials] = pre_processing_pip_trials(V4_dataset,params.bpfilt,params.bpwidth,params.order);
V4_trials_nbp = pre_processing_pip_trials(V4_dataset,false,params.bpwidth,params.order);

% fft and hilbert 
in_fft = do_fft(in_trials,params)
V4_fft = do_fft(V4_trials,params)
in_fft_nbp = do_fft(in_trials_nbp,params)
V4_fft_nbp = do_fft(V4_trials_nbp,params)

[hilb_angles.wrapped.in]  = pre_processing_pip_hilb(in_trials);
[hilb_angles.wrapped.V4]  = pre_processing_pip_hilb(V4_trials);

[hilb_angles.in,inst_freq.in,filt_data.in]  = pre_processing_pip_filtinst(hilb_angles.wrapped.in,params.filttype,params.framelen,params.filtord,params.toi);
[hilb_angles.V4,inst_freq.V4,filt_data.V4]  = pre_processing_pip_filtinst(hilb_angles.wrapped.V4,params.filttype,params.framelen,params.filtord,params.toi);

filt_inst_freq.in = filt_data.in;
filt_inst_freq.V4 = filt_data.V4;

% Averaging
[inst_avg.in inst_gavg.in] = do_grand_avg(filt_inst_freq.in);
[inst_avg.V4 inst_gavg.V4] = do_grand_avg(filt_inst_freq.V4);

for ii = 1:length(inst_avg.in)
    s_avg.in(ii) = mean(inst_avg.in(ii).avg);
    s_avg.V4(ii) = mean(inst_avg.V4(ii).avg);
end 
fft_gavg.in = do_mean_pow(in_fft)
fft_gavg.V4 = do_mean_pow(V4_fft)
fft_gavg_nbp.in = do_mean_pow(in_fft_nbp)
fft_gavg_nbp.V4 = do_mean_pow(V4_fft_nbp)

b = [s_avg.in;s_avg.V4];
%%
f = figure;
f.Units = 'normalized';
f.Position = [0.25 0.25 0.5 0.7]
sgtitle(sprintf('Inst frequencies and fft using a BP-filter with order %i and bounds %i - %i Hz',params.order,params.lower,params.upper));
subplot(2,1,1)
boxplot(b','Labels',{'V1a','V4'})
ylabel('Frequency [Hz]')

subplot(2,2,3)
plot(in_fft(1).freq,fft_gavg.in)
hold on 
plot(in_fft_nbp(1).freq,fft_gavg_nbp.in)
hold off
title('Power spectrum V1a')
xlabel('Frequency [Hz]')
ylabel('Power')
ylim([0 max(fft_gavg.in)*3])
subplot(2,2,4)
plot(V4_fft(1).freq,fft_gavg.V4);
hold on 
plot(V4_fft_nbp(1).freq,fft_gavg_nbp.V4)
ylim([0 max(fft_gavg.V4)*3])
title('Power spectrum V4')
xlabel('Frequency [Hz]')
ylabel('Power')
foldername = fullfile(params.figpath,'Filter_testing')
if ~exist(foldername,'dir')
    mkdir(foldername)
end
saveas(f,fullfile(foldername,sprintf('filter_inst_fft_order%i_bounds%i_%i.fig',params.order,params.lower,params.upper)))
saveas(f,fullfile(foldername,sprintf('filter_inst_fft_order%i_bounds%i_%i.jpg',params.order,params.lower,params.upper)))

%% Pip with mannual filter 
clc
clear 
% Trial selection 
% Loading
load("/data/projects/V1V4coherence/02_analysis_max/git_repos/params.mat")
load('V4_dataset.mat')
% Trialselection 
[V4_trials] = pre_processing_pip_trials(V4_dataset,params.bpfilt,params.bpwidth)
[hilb_angles.wrapped.V4]  = pre_processing_pip_hilb(V4_trials);
[hilb_angles.V4,inst_freq.V4,filt_data.V4]  = pre_processing_pip_filtinst(hilb_angles.wrapped.V4,params.filttype,params.framelen,params.filtord,params.toi);
filt_inst_freq.V4 = filt_data.V4;
[V4_avg V4_g_avg] = do_grand_avg(filt_inst_freq.V4);
for ii = 1:length(V4_avg)
    s_avg(ii) = mean(V4_avg(ii).avg);
end 
f = figure;
boxplot(s_avg)


%%
function mean_pow = do_mean_pow(fft_struct)
for i_s = 1:length(fft_struct)
    f_pow(i_s,:) = fft_struct(i_s).powspctrm(1,:);
end 
mean_pow = mean(f_pow);
end 