% Script to test the envelope information of the hilbert transform
clc
clear 
% Trial selection 
% Loading
load("/data/projects/V1V4coherence/02_analysis_max/git_repos/params.mat")
load('attin_dataset.mat')
load('V4_dataset.mat')

[in_trials] = pre_processing_pip_trials(attin_dataset,params.bpfilt,params.bpwidth,params.order);
[V4_trials] = pre_processing_pip_trials(V4_dataset,params.bpfilt,params.bpwidth,params.order);

%% generating fake data
amp_range = [0.1 100];
noise_mod = 0.1;
freqs = [100]
f_in_trials = do_create_data(in_trials,freqs,amp_range,noise_mod)
f_V4_trials = do_create_data(V4_trials,freqs,amp_range,noise_mod)
% re bpfiltering
in_trials = do_bpfilt(f_in_trials,params);
V4_trials = do_bpfilt(f_in_trials,params);

%
[hilb_angles.wrapped.in,hilb_env.in]  = pre_processing_pip_hilb(in_trials);
[hilb_angles.wrapped.V4,hilb_env.V4]  = pre_processing_pip_hilb(V4_trials);

[hilb_angles.in,inst_freq.in,filt_data.in]  = pre_processing_pip_filtinst(hilb_angles.wrapped.in,params.filttype,params.framelen,params.filtord,params.toi);
[hilb_angles.V4,inst_freq.V4,filt_data.V4]  = pre_processing_pip_filtinst(hilb_angles.wrapped.V4,params.filttype,params.framelen,params.filtord,params.toi);

hilb_env.in = do_toi_cut(hilb_env.in,params.toi);
hilb_env.V4 = do_toi_cut(hilb_env.V4,params.toi);
hilb_angles.wrapped.in = do_toi_cut(hilb_angles.wrapped.in,params.toi);
hilb_angles.wrapped.V4 = do_toi_cut(hilb_angles.wrapped.V4,params.toi);

[mean_env.in.avg mean_env.in.g_avg] = do_grand_avg(hilb_env.in);
[mean_env.V4.avg mean_env.V4.g_avg] = do_grand_avg(hilb_env.V4);

[mean_inst.in.avg mean_inst.in.g_avg] = do_grand_avg(filt_data.in);
[mean_inst.V4.avg mean_inst.V4.g_avg] = do_grand_avg(filt_data.V4);


% Calculating for each point in time the respective frequency and amplitude
inst_trials.in = cat_all_trials(filt_data.in);
inst_trials.V4 = cat_all_trials(filt_data.V4);
env_trials.in = cat_all_trials(hilb_env.in);
env_trials.V4 = cat_all_trials(hilb_env.V4);
bins = [100,100];
%

% binning
cur_inst = inst_trials.in;
cur_env = env_trials.in;
[N,x1,y1] = histcounts2(cur_inst,cur_env,bins);
% N_norm = N ./ sum(N,1);
N_norm = N ./ max(N);
N_norm(isnan(N_norm)) = 0;


%% PLotting and calcuating heat map of power vs inst freq: V1a: Normalized
cur_inst = inst_trials.in;
cur_env = env_trials.in;
[N,x1,y1] = histcounts2(cur_inst,cur_env,bins);
% N_norm = N ./ sum(N,1);
N_norm = N ./ max(N);
N_norm(isnan(N_norm)) = 0;

f = figure;
f.Units = 'normalized';
f.Position = [0.25 0.25 0.7 0.7]
imagesc(y1(1:end-1),x1(1:end-1),N_norm)
title(sprintf('Heatmap of V1a instantaneous power vs instantaneous frequency, normalized over power: BP %i - %i Hz',params.lower,params.upper))
xlabel('Power')
ylabel('Frequency [Hz]')
xlim([0 4])
axis xy
% hold on 
% plot(1:60)
colorbar()
foldername = fullfile(params.figpath,'envelope_testing');
if ~exist(foldername,'dir')
    mkdir(foldername);
end 
% saveas(f,fullfile(foldername,sprintf('pow_vs_inst_V1a_norm_%i_%i.fig',params.lower,params.upper)));
% saveas(f,fullfile(foldername,sprintf('pow_vs_inst_V1a_norm_%i_%i.jpg',params.lower,params.upper)));


%% PLotting and calcuating heat map of power vs inst freq: V1a
cur_inst = inst_trials.in;
cur_env = env_trials.in;
[N,x1,y1] = histcounts2(cur_inst,cur_env,bins);
% N_norm = N ./ sum(N,1);
N_norm = N ./ max(N);
N_norm(isnan(N_norm)) = 0;

f = figure;
f.Units = 'normalized';
f.Position = [0.25 0.25 0.7 0.7]
imagesc(y1(1:end-1),x1(1:end-1),N)
xlim([0 2])
title(sprintf('Heatmap of V1a instantaneous power vs instantaneous frequency: BP %i - %i Hz',params.lower,params.upper))
xlabel('Power')
ylabel('Frequency [Hz]')
axis xy
colorbar()
foldername = fullfile(params.figpath,'envelope_testing');
if ~exist(foldername,'dir')
    mkdir(foldername);
end 
% saveas(f,fullfile(foldername,sprintf('pow_vs_inst_V1a_%i_%i.fig',params.lower,params.upper)));
% saveas(f,fullfile(foldername,sprintf('pow_vs_inst_V1a_%i_%i.jpg',params.lower,params.upper)));

%% PLotting and calcuating heat map of power vs inst freq: V1a: Normalized
cur_inst = inst_trials.V4;
cur_env = env_trials.V4;
[N,x1,y1] = histcounts2(cur_inst,cur_env,bins);
% N_norm = N ./ sum(N,1);
N_norm = N ./ max(N);
N_norm(isnan(N_norm)) = 0;

f = figure;
f.Units = 'normalized';
f.Position = [0.25 0.25 0.7 0.7]
imagesc(y1(1:end-1),x1(1:end-1),N_norm)
title(sprintf('Heatmap of V4 instantaneous power vs instantaneous frequency, normalized over power: BP %i - %i Hz',params.lower,params.upper))
xlabel('Power')
ylabel('Frequency [Hz]')
axis xy
colorbar()
foldername = fullfile(params.figpath,'envelope_testing');
if ~exist(foldername,'dir')
    mkdir(foldername);
end 
% saveas(f,fullfile(foldername,sprintf('pow_vs_inst_V4_norm_%i_%i.fig',params.lower,params.upper)));
% saveas(f,fullfile(foldername,sprintf('pow_vs_inst_V4_norm_%i_%i.jpg',params.lower,params.upper)));


%% PLotting and calcuating heat map of power vs inst freq: V4
cur_inst = inst_trials.V4;
cur_env = env_trials.V4;
[N,x1,y1] = histcounts2(cur_inst,cur_env,bins);
% N_norm = N ./ sum(N,1);
N_norm = N ./ max(N);
N_norm(isnan(N_norm)) = 0;

f = figure;
f.Units = 'normalized';
f.Position = [0.25 0.25 0.7 0.7]
imagesc(y1(1:end-1),x1(1:end-1),N)
title(sprintf('Heatmap of V4 instantaneous power vs instantaneous frequency: BP %i - %i Hz',params.lower,params.upper))
xlabel('Power')
ylabel('Frequency [Hz]')
axis xy
colorbar()
foldername = fullfile(params.figpath,'envelope_testing');
if ~exist(foldername,'dir')
    mkdir(foldername);
end 
saveas(f,fullfile(foldername,sprintf('pow_vs_inst_V4_%i_%i.fig',params.lower,params.upper)));
saveas(f,fullfile(foldername,sprintf('pow_vs_inst_V4_%i_%i.jpg',params.lower,params.upper)));
%% recalculating signal from abs and ang
% reconsturicting cos of the phase 
rec_signal = hilb_env.in;
for i_s = 1:length(hilb_angles.wrapped.in)
    for i_t = 1:length(hilb_angles.wrapped.in(i_s).trial)
        for i_c = 1:size(hilb_angles.wrapped.in(i_s).trial{1,i_t},1)
            cur_pha = hilb_angles.wrapped.in(i_s).trial{1,i_t}(i_c,:);
            cur_env = hilb_env.in(i_s).trial{1,i_t}(i_c,:);
            rec_sig = cur_env .* cos(cur_pha);
            rec_signal(i_s).trial{1,i_t}(i_c,:) = rec_sig;
        end 
    end 
end 

% plotting side by side the original and reconstructed signal 
f = figure;
sess = 1;
f.Units = 'normalized';
f.Position = [0 0 1 1]

for i_s = 1:length(in_trials(sess).trial)
    subplot(2,1,1)
    plot(in_trials(sess).trial{1,i_s}(1,:))
    subplot(2,1,2)
    plot(rec_signal(sess).trial{1,i_s}(1,:))    
    w = waitforbuttonpress;
end 
%%
subplot(2,1,1)
plot(mean_env.g_avg.avg)
subplot(2,1,2)
plot(mean_inst.g_avg.avg)

%% functions
function full_data = cat_all_trials(trials)
trials = [trials.trial];
trials = cellfun(@(x) x(1,:),trials,'UniformOutput',false);
full_data = cell2mat(trials);
end 

 