clc
clear all
load('attout_dataset.mat')
load('attin_dataset.mat')
load("V4_dataset.mat")

% Trial preprocessing parameters
bpfilt = true;
bpwidth = [50 90];

% SSD parameters
toi = [2.3 4.3]; % Time region of interest [2.3 4.3] is MC 2&3
fs = 1000; % Sampling frequency of the signal
th = 0.01; % residual variance threshhold   

% Hilbert parameters
filttype = "sgolay"; %either medfilt or sgolay
framelen = 31;
filtord = 1;

matpath = '/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files';

% Trialselection 
[in_trials] = pre_processing_pip_trials(attin_dataset,bpfilt,bpwidth,toi)
[out_trials] = pre_processing_pip_trials(attout_dataset,bpfilt,bpwidth,toi)
[V4_trials] = pre_processing_pip_trials(V4_dataset,bpfilt,bpwidth,toi)

% Performing SSD
[in_trials] = do_SSD(in_trials,fs,th)
[out_trials] = do_SSD(out_trials,fs,th)
[V4_trials] = do_SSD(V4_trials,fs,th)

% Hilbert Angles, instantaneous frequency, filtered instantaneous
% frequency, instantaneous change, summary struct
[hilb_angles.wrapped.in,hilb_angles.in,inst_freq.in,filt_inst_freq.in,inst_change.in,summary.in]  = pre_processing_pip_hilb(in_trials,filttype,framelen,filtord);
[hilb_angles.wrapped.out,hilb_angles.out,inst_freq.out,filt_inst_freq.out,inst_change.out,summary.out]  = pre_processing_pip_hilb(out_trials,filttype,framelen,filtord);
[hilb_angles.wrapped.V4,hilb_angles.V4,inst_freq.V4,filt_inst_freq.V4,inst_change.V4,summary.V4]  = pre_processing_pip_hilb(V4_trials,filttype,framelen,filtord);
timebar = -1.3:0.001:5;


%% 
in_plv = nan(23,252); %creating nan array of all sessions and max trials
out_plv = nan(23,231); % differs from in because different max 
count = 1;
for i_sess = 1:length(hilb_angles.wrapped.out)
    out_phase_diff = cell(1,length(hilb_angles.wrapped.out(i_sess).trial));
    for i_chan = 1:length(hilb_angles.wrapped.out(i_sess).label)-1
        for i_trial = 1:length(hilb_angles.wrapped.out(i_sess).trial)
            signal_1 = hilb_angles.wrapped.out(i_sess).trial{1,i_trial}(1,:);
            signal_2 = hilb_angles.wrapped.out(i_sess).trial{1,i_trial}(i_chan+1,:);
            out_phase_diff{:,i_trial} = signal_1-signal_2;
            out_plv(count,i_trial) = abs(mean(exp(1j * out_phase_diff{:,i_trial})));
        end 
        count = count + 1;
    end 
end 
count = 1;
for i_sess = 1:length(hilb_angles.wrapped.in)
    in_phase_diff = cell(1,length(hilb_angles.wrapped.in(i_sess).trial));
    for i_chan = 1:length(hilb_angles.wrapped.in(i_sess).label)-1
        for i_trial = 1:length(hilb_angles.wrapped.in(i_sess).trial)
            signal_1 = hilb_angles.wrapped.in(i_sess).trial{1,i_trial}(1,:);
            signal_2 = hilb_angles.wrapped.in(i_sess).trial{1,i_trial}(i_chan+1,:);
            in_phase_diff{:,i_trial} = signal_1-signal_2;
            in_plv(count,i_trial) = abs(mean(exp(1j * in_phase_diff{:,i_trial})));
        end 
        count = count + 1;
    end 
end 
in_plv_m = mean(in_plv,2,'omitnan');
out_plv_m = mean(out_plv,2,'omitnan');
in_plv_fin = mean(in_plv_m);
out_plv_fin = mean(out_plv_m);


%% 
new = cat(2,in_plv_m,out_plv_m);
figure;
boxplot(new,'Labels',{'PLV attended','PLV unattended'})
ylabel('PLV')
title('PLV in the frequency range from 52 to 78 Hz using conjugate method')
ylim([0.1 0.25])
[p, value] = signrank(in_plv_m, out_plv_m)

%% Calculating PC the Iris way
comb_data = hilb_angles.wrapped.in;
counter = 1;
for i_sess = 1:length(comb_data)
    var = comb_data(i_sess);
    for i_chan = 1:length(var.label)-1
        contrib = zeros(1,2001);
        dCos = zeros(1,2001);
        dSin = zeros(1,2001);
        for i_trial = 1:length(var.trial)
            temp_array = NaN(1,2001);
            temp_array(:,1:length(var.trial{i_trial})) = 1;
            phase_diff = var.trial{i_trial}(1,:) - var.trial{i_trial}(i_chan+1,:);
            contrib = contrib + ~isnan(temp_array);
            dCos(:,1:length(var.trial{i_trial})) = dCos(:,1:length(var.trial{i_trial})) + cos(phase_diff);
            dSin(:,1:length(var.trial{i_trial})) = dSin(:,1:length(var.trial{i_trial})) + sin(phase_diff);
        end 
        dCos = dCos./contrib;
        dSin = dSin./contrib;
        VL = sqrt((dCos).^2+(dSin).^2);
        VLEXP = sqrt(pi)./(2.*sqrt(contrib));
        VL = VL - VLEXP;
        meanAngle = acos((dCos)./VL);
        in_meanVL(counter) = mean(VL);
        counter = counter + 1;
    end 
end 
in_mean_meanVL = mean(in_meanVL)

comb_data = hilb_angles.wrapped.out;
counter = 1;
for i_sess = 1:length(comb_data)
    var = comb_data(i_sess);
    for i_chan = 1:length(var.label)-1
        contrib = zeros(1,2001);
        dCos = zeros(1,2001);
        dSin = zeros(1,2001);
        for i_trial = 1:length(var.trial)
            temp_array = NaN(1,2001);
            temp_array(:,1:length(var.trial{i_trial})) = 1;
            phase_diff = var.trial{i_trial}(1,:) - var.trial{i_trial}(i_chan+1,:);
            contrib = contrib + ~isnan(temp_array);
            dCos(:,1:length(var.trial{i_trial})) = dCos(:,1:length(var.trial{i_trial})) + cos(phase_diff);
            dSin(:,1:length(var.trial{i_trial})) = dSin(:,1:length(var.trial{i_trial})) + sin(phase_diff);
        end 
        dCos = dCos./contrib;
        dSin = dSin./contrib;
        VL = sqrt((dCos).^2+(dSin).^2);
        VLEXP = sqrt(pi)./(2.*sqrt(contrib));
        VL = VL - VLEXP;
        meanAngle = acos((dCos)./VL);
        out_meanVL(counter) = mean(VL);
        counter = counter + 1;
    end 
end 
out_mean_meanVL = mean(out_meanVL)

%%
new = cat(1,in_meanVL,out_meanVL);
[p, value] = signrank(in_meanVL', out_meanVL')
figure;
boxplot(new','Labels',{'PLV attended','PLV unattended'})
ylabel('PLV')
ylim([0 0.21])
title('PLVs per trial, BP-filter 52 - 78 Hz')

