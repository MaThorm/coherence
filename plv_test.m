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


%% Creating trials and bp filtering !!!Important!!! I do not use all V4/V1 combinations yet
in_plv = nan(23,252); %creating nan array of all sessions and max trials
out_plv = nan(23,231); % differs from in because different max 
count = 1;
for i_sess = 1:length(angles.out_hilbert)
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
%%
in_plv_m = mean(in_plv,2,'omitnan');
out_plv_m = mean(out_plv,2,'omitnan');
in_plv_fin = mean(in_plv_m);
out_plv_fin = mean(out_plv_m);


%%
new = cat(2,in_plv_m,out_plv_m)
boxplot(new,'Labels',{'PLV attended','PLV unattended'})
ylabel('PLV')
title('PLV in the frequency range from 52 to 78 Hz')

%violin(new)
