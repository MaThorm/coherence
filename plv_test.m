clc
clear all
dir = "/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files";
toi = [2.3 4.3];
bp = [52 78];
medfiltord = 20;
load('attout_dataset.mat')
load('attin_dataset.mat')
load("V4_dataset.mat")
[in_trials,out_trials,V4_trials] = pre_processing_pip_trials(attin_dataset,attout_dataset,V4_dataset,bp,toi);
[grand_struct,angles,inst_freq] = pre_processing_pip_hilb(in_trials,out_trials,V4_trials,medfiltord)


%% Creating trials and bp filtering
in_plv = nan(23,252); %creating nan array of all sessions and max trials
out_plv = nan(23,231); % differs from in because different max 
count = 1;
for i_sess = 1:length(angles.out_hilbert)
    out_phase_diff = cell(1,length(angles.out_hilbert(i_sess).trial));
    for i_chan = 1:length(angles.out_hilbert(i_sess).label)-1
        for i_trial = 1:length(angles.out_hilbert(i_sess).trial)
            signal_1 = angles.out_hilbert(i_sess).trial{1,i_trial}(1,:);
            signal_2 = angles.out_hilbert(i_sess).trial{1,i_trial}(i_chan+1,:);
            out_phase_diff{:,i_trial} = signal_1-signal_2;
            out_plv(count,i_trial) = abs(mean(exp(1j * out_phase_diff{:,i_trial})));
        end 
        count = count + 1;
    end 
end 
count = 1;
for i_sess = 1:length(angles.in_hilbert)
    in_phase_diff = cell(1,length(angles.in_hilbert(i_sess).trial));
    for i_chan = 1:length(angles.in_hilbert(i_sess).label)-1
        for i_trial = 1:length(angles.in_hilbert(i_sess).trial)
            signal_1 = angles.in_hilbert(i_sess).trial{1,i_trial}(1,:);
            signal_2 = angles.in_hilbert(i_sess).trial{1,i_trial}(i_chan+1,:);
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
title('PLV in the frequency range from 52 to 78 Hz')
[p, value] = signrank(in_plv_m, out_plv_m)


%% Calculating PC the Iris way
for i_sess = 1:length(angles.in_hilbert)
    contrib = zeros(1,2001);
    dCos = zeros(1,2001);
    dSin = zeros(1,2001);
    var = angles.in_hilbert(i_sess);
    for i_trial = 1:length(var.trial)
        temp_array = NaN(1,2001);
        temp_array(:,1:length(var.trial{i_trial})) = 1;
        phase_diff = var.trial{i_trial}(1,:) - var.trial{i_trial}(2,:);
        contrib = contrib + ~isnan(temp_array);
        dCos(:,1:length(var.trial{i_trial})) = dCos(:,1:length(var.trial{i_trial})) + cos(phase_diff);
        dSin(:,1:length(var.trial{i_trial})) = dSin(:,1:length(var.trial{i_trial})) + cos(phase_diff);
    end 
    dCos = dCos./contrib;
    dSin = dSin./contrib;
    VL = sqrt((dCos).^2+(dSin).^2);
    VLEXP = sqrt(pi)./(2.*sqrt(contrib));
    VL = VL - VLEXP;
    meanAngle = acos((dCos)./VL);
    meanVL(i_sess) = mean(VL);
end 
mean_meanVL = mean(meanVL)