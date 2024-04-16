clc
clear all
load('attout_dataset.mat')
load('attin_dataset.mat')
load("V4_dataset.mat")

% Trial preprocessing parameters
bpfilt = true;
bpwidth = [30 100];

% SSD parameters
toi = [3.3 4.3]; % Time region of interest [2.3 4.3] is MC 2&3
fs = 1000; % Sampling frequency of the signal
th = 0.01; % residual variance threshhold   

% Hilbert parameters
filttype = "sgolay"; %either medfilt or sgolay
framelen = 31;
filtord = 1;

matpath = '/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files';

%
[in_trials] = pre_processing_pip_trials(attin_dataset,bpfilt,bpwidth,toi)
[out_trials] = pre_processing_pip_trials(attout_dataset,bpfilt,bpwidth,toi)
[V4_trials] = pre_processing_pip_trials(V4_dataset,bpfilt,bpwidth,toi)

% Performing SSD
[in_trials] = do_SSD(in_trials,fs,th)
[out_trials] = do_SSD(out_trials,fs,th)
[V4_trials] = do_SSD(V4_trials,fs,th)
% Hilbert Angles, instantaneous frequency, filtered instantaneous
% frequency, instantaneous change, summary struct
[hilb_angles.in,inst_freq.in,filt_inst_freq.in,inst_change.in,summary.in]  = pre_processing_pip_hilb(in_trials,filttype,framelen,filtord);
[hilb_angles.out,inst_freq.out,filt_inst_freq.out,inst_change.out,summary.out]  = pre_processing_pip_hilb(out_trials,filttype,framelen,filtord);
[hilb_angles.V4,inst_freq.V4,filt_inst_freq.V4,inst_change.V4,summary.V4]  = pre_processing_pip_hilb(V4_trials,filttype,framelen,filtord);
timebar = -1.3:0.001:5;

%%
attin_inst = grand_struct.in_medfiltHilbert;
attout_inst = grand_struct.out_medfiltHilbert;
insummary = grand_struct.insummary;
outsummary = grand_struct.outsummary;
V4_inst = grand_struct.V4_medfiltHilbert;
V4summary = grand_struct.V4summary;

%% Plotting all attin sessions with error bars 
input = filt_inst_freq.V4; %Change last parameter: in,out,V4
figure;
sgtitle('Instantaneous frequencies per recording site')
sel = 1:length(input(1).sess_mean);
x = timebar(sel);
for ii = 1:length(input)
    subplot(4,4,ii)
    y = input(ii).sess_mean(sel);
    sd = input(ii).sess_sd(sel);
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
    hold on
    plot(x,y,'r')
    hold off
end 

%% Plotting AttIn Summary
figure; 
sel = 1:length(summary.in.mean);
x = timebar(sel)
y = summary.in.mean(sel);
sd = summary.in.std(sel);
patch([x fliplr(x)], [y-sd  fliplr(y+sd)], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
title('AttIn Summary')
xlabel('Time [s]')
ylabel('Freqeuncy [Hz]')
hold on 
label = {'Static', 'MS1','MS2','MS3','MS4'}
xl = xline([-0.5 0 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
plot(x,y,'r');
hold off
%% Plotting attout Summary
sel = 1:length(summary.out.mean);
x = timebar(sel)
figure
y = summary.out.mean(sel);
sd = summary.out.std(sel);
patch([x fliplr(x)], [y-sd  fliplr(y+sd)], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
label = {'Static', 'MS1','MS2','MS3','MS4'}
xl = xline([-0.5 0 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
xlabel('Time [s]')
ylabel('Freqeuncy [Hz]')
title('AttOut Summary')
hold on 
plot(x,y,'r');
hold off

%% plotting V4 summary 
figure
sel = 1:length(summary.V4.mean);
x = timebar(sel)
y = summary.V4.mean(sel);
sd = summary.V4.std(sel);
patch([x fliplr(x)], [y-sd  fliplr(y+sd)], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
label = {'Static', 'MS1','MS2','MS3','MS4'}
xl = xline([-0.5 0 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
xlabel('Time [s]')
ylabel('Freqeuncy [Hz]')%dividing the remaining variance of the largest component by the remaining variance of all other 
title('V4 Summary')
hold on 
plot(x,y,'r');
hold off

%% plotting all together 
sel = 1:length(summary.out.mean);
x = timebar(sel);
figure
plot(x,summary.in.mean(sel),'r');
hold on 
plot(x,summary.out.mean(sel),'b');
plot(x,summary.V4.mean(sel));
title("Summary of all conditions")
legend('AttIn','AttOut','V4','AutoUpdate','off')
label = {'Static', 'MS1','MS2','MS3','MS4'};
%xl = xline([-0.5 0 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
%ylim([60 80])
hold off
xlabel('Time [s]')
ylabel('Frequency [Hz]')

%% plotting all together, but only MS 2 & 3 
sel = 1:length(summary.in.mean);
x = timebar(sel);
figure;
plot(x,summary.in.mean,'r');
hold on 
plot(x,summary.out.mean,'b');
plot(x,summary.V4.mean);
title("Summary of all conditions")
legend('AttIn','AttOut','V4','AutoUpdate','off')
label = {'MS2','MS3','MS4'};
%xl = xline([ 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
hold off
xlabel('Time [s]')
ylabel('Frequency [Hz]')

%% Quick derivative check
a = diff(summary.in.mean(sel))
plot(a)
%% Plotting derivative of derivative over all recording sites
sel = 1:5000;
x = timebar(sel);
subplot(3,1,1)
ylim([-2 2])
plot(x,insummary.change_mean(sel),'r')
hold on 
title('AttIn')
subplot(3,1,2)
plot(x,outsummary.change_mean(sel),'b')
title('AttOut')
subplot(3,1,3)
plot(x,V4summary.change_mean(sel))
title('V4')

%% more derivative testing
for ii = 1:length(V4_inst)
    cfg = [];
    cfg.channel = 'all';
    cfg.derivative = 'yes';
    instchange(ii) = ft_preprocessing(cfg,V4_inst(ii));
end 
%% 
sel = 1:5000;
for ii = 1:length(V4_instchange)
    subplot(4,4,ii)
    plot(timebar(sel),V4_instchange(ii).sess_mean(sel));
end 
%% even more derivative testing
for ii = 1:length(V4_instchange)
    mean_dv{ii} = cent_diff_n(V4_inst(ii).sess_mean,1,5);
end 
maxNumCol = max(cellfun(@(c) size(c,2), mean_dv));
aMat = cell2mat(cellfun(@(c){[c nan(1,maxNumCol-numel(c))]}, mean_dv)');
colMeans = mean(aMat,1,'omitnan');
new_instchange.V4 = colMeans;

for ii = 1:length(attin_inst)
    mean_dv{ii} = cent_diff_n(attin_inst(ii).sess_mean,1,5);
end 
maxNumCol = max(cellfun(@(c) size(c,2), mean_dv));
aMat = cell2mat(cellfun(@(c){[c nan(1,maxNumCol-numel(c))]}, mean_dv)');
colMeans = mean(aMat,1,'omitnan');
new_instchange.in = colMeans;

for ii = 1:length(attout_inst)
    mean_dv{ii} = cent_diff_n(attout_inst(ii).sess_mean,1,5);
end 
maxNumCol = max(cellfun(@(c) size(c,2), mean_dv));
aMat = cell2mat(cellfun(@(c){[c nan(1,maxNumCol-numel(c))]}, mean_dv)');
colMeans = mean(aMat,1,'omitnan');
new_instchange.out = colMeans;
%
figure
subplot(3,1,1)
title('AttIn')
sel = 1:5000;
yyaxis left
plot(timebar(sel),new_instchange.in(sel));
ylabel('Change in Frequency')
hold on 
ylim([-0.15 0.15])
yyaxis right
plot(timebar(sel),insummary.mean(sel));
ylabel('Frequency [Hz]')
% 
subplot(3,1,2)
title('AttOut')

yyaxis left
plot(timebar(sel),new_instchange.out(sel));
ylabel('Change in Frequency')
hold on 
ylim([-0.15 0.15])
yyaxis right
plot(timebar(sel),outsummary.mean(sel));
ylabel('Frequency [Hz]')
%
subplot(3,1,3)
title('V4')

yyaxis left
plot(timebar(sel),new_instchange.V4(sel));
ylabel('Change in Frequency')
hold on 
ylim([-0.15 0.15])
yyaxis right
plot(timebar(sel),V4summary.mean(sel));
ylabel('Frequency [Hz]')
 
%%