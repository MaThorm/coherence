clc
clear all
%%


saving = true;

% Trial preprocessing parameters
bpfilt = false;
bpwidth = [1 150];
toi = [2.3 4.3]; % Time region of interest [2.3 4.3] is MC 2&3


% SSD parameters
fs = 1000; % Sampling frequency of the signal
th = 0.01; % residual variance threshhold   
lower = 50;
upper = 80;

% Hilbert parameters
filttype = "sgolay"; %either medfilt or sgolay
framelen = 31;
filtord = 1;
timebar = -1.3:0.001:5;


matpath = '/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files';
figpath = '/data/projects/V1V4coherence/03_results_max';
%
% loading inst freq struct 
filename = sprintf("inst_freq_thr%.2f_comp%d_lower%d_upper%d.mat",th,10,lower,upper);
m = matfile(fullfile(matpath,'Inst_freq',filename),'Writable',true);
angles = m.angle;



%% Calculating PC the Iris way
comb_data = angles.in;
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

comb_data = angles.out;
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

%
new = cat(1,in_meanVL,out_meanVL);
[p, value] = signrank(in_meanVL', out_meanVL')
figure;
boxplot(new','Labels',{'PLV attended','PLV unattended'})
ylabel('PLV')
ylim([0 0.21])
title(sprintf('PLVs per trial, SSD Boundaries: %d, %d',lower,upper))
saveas(gcf,fullfile(figpath,sprintf('PLV/PLV_SSDbounds%d-%d.fig',lower,upper)))
saveas(gcf,fullfile(figpath,sprintf('PLV/PLV_SSDbounds%d-%d.png',lower,upper)))

%%

in_plv = nan(23,252); %creating nan array of all sessions and max trials
out_plv = nan(23,231); % differs from in because different max 
count = 1;
for i_sess = 1:length(angles.out)
    out_phase_diff = cell(1,length(angles.out(i_sess).trial));
    for i_chan = 1:length(angles.out(i_sess).label)-1
        for i_trial = 1:length(angles.out(i_sess).trial)
            signal_1 = angles.out(i_sess).trial{1,i_trial}(1,:);
            signal_2 = angles.out(i_sess).trial{1,i_trial}(i_chan+1,:);
            out_phase_diff{:,i_trial} = signal_1-signal_2;
            out_plv(count,i_trial) = abs(mean(exp(1j * out_phase_diff{:,i_trial})));
        end 
        count = count + 1;
    end 
end 
count = 1;
for i_sess = 1:length(angles.in)
    in_phase_diff = cell(1,length(angles.in(i_sess).trial));
    for i_chan = 1:length(angles.in(i_sess).label)-1
        for i_trial = 1:length(angles.in(i_sess).trial)
            signal_1 = angles.in(i_sess).trial{1,i_trial}(1,:);
            signal_2 = angles.in(i_sess).trial{1,i_trial}(i_chan+1,:);
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
