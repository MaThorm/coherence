clc
clear all
load('attout_dataset.mat')
load('attin_dataset.mat')
load("V4_dataset.mat")
bpwidth = [55 95];
toi = [-100 100];
medianfiltord = 20;
matpath = '/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files'

[in_trials,out_trials,V4_trials] = pre_processing_pip_trials(attin_dataset,attout_dataset,V4_dataset,bpwidth,toi)
[grand_struct,angles,inst_freq] = pre_processing_pip_hilb(in_trials,out_trials,V4_trials,medianfiltord)
attin_inst = grand_struct.in_medfiltHilbert;
attout_inst = grand_struct.out_medfiltHilbert;
insummary = grand_struct.insummary;
outsummary = grand_struct.outsummary;
V4_inst = grand_struct.V4_medfiltHilbert;
V4summary = grand_struct.V4summary;
timebar = -1.3:0.001:5;

%% Plotting all attin sessions with error bars 
figure;
sgtitle('Instantaneous frequencies per recording site AttIn')
sel = 1:5000;
x = timebar(sel);
for ii = 1:length(attout_inst)
    subplot(4,4,ii)
    y = attin_inst(ii).sess_mean(sel);
    sd = attin_inst(ii).sess_sd(sel);
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
    hold on
    plot(x,y,'r')
    hold off
end 

%% Plotting all attout sessions with error bars 
figure;
sgtitle('Instantaneous frequencies per recording site AttOut')
sel = 1:5000;
x = timebar(sel);
for ii = 1:length(attout_inst)
    subplot(4,4,ii)
    y = attout_inst(ii).sess_mean(sel);
    sd = attout_inst(ii).sess_sd(sel);
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
    hold on
    plot(x,y,'r')
    hold off
end 

%% Plotting AttIn Summary
figure; 
sel = 1:5000;
x = timebar(sel)
y = insummary.mean(sel);
sd = insummary.std(sel);
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
sel = 1:5000;
x = timebar(sel)
figure
y = outsummary.mean(sel);
sd = outsummary.std(sel);
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
sel = 1:5000;
x = timebar(sel)
y = V4summary.mean(sel);
sd = V4summary.std(sel);
patch([x fliplr(x)], [y-sd  fliplr(y+sd)], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
label = {'Static', 'MS1','MS2','MS3','MS4'}
xl = xline([-0.5 0 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
xlabel('Time [s]')
ylabel('Freqeuncy [Hz]')
title('V4 Summary')
hold on 
plot(x,y,'r');
hold off

%% plotting all together 
sel = 1:5000;
x = timebar(sel);
figure
plot(x,insummary.mean(sel),'r');
hold on 
plot(x,outsummary.mean(sel),'b');
plot(x,V4summary.mean(sel));
title("Summary of all conditions")
legend('AttIn','AttOut','V4','AutoUpdate','off')
label = {'Static', 'MS1','MS2','MS3','MS4'};
xl = xline([-0.5 0 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
ylim([60 80])
hold off
xlabel('Time [s]')
ylabel('Frequency [Hz]')

%% plotting all together, but only MS 2 & 3 
sel = 2000:4500;
x = timebar(sel);
figure;
plot(x,insummary.mean(sel),'r');
hold on 
plot(x,outsummary.mean(sel),'b');
plot(x,V4summary.mean(sel));
title("Summary of all conditions")
xlim([0.7 3.2])
legend('AttIn','AttOut','V4','AutoUpdate','off')
label = {'MS2','MS3','MS4'};x = timebar(sel);

xl = xline([ 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
hold off
xlabel('Time [s]')
ylabel('Frequency [Hz]')

%% plotting all together with combined median filter orders
sel = 1:5000;
x = timebar(sel);
figure;
plot(x,summary_struct.in_mean(sel),'r');
hold on 
plot(x,summary_struct.out_mean(sel),'b');
plot(x,summary_struct.V4_mean(sel));
title("Summary of all conditio0s,multiple median filters")
legend('AttIn','AttOut','V4','AutoUpdate','off')
label = {'Static', 'MS1','MS2','MS3','MS4'};
xl = xline([-0.5 0 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
hold off
xlabel('Time [s]')
ylabel('Frequency [Hz]')


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

%% plotting some shit
for ii = 1:length(attin_instchange)
    subplot(4,4,ii)
    plot(attin_instchange(ii).sess_mean)
end 
%% 
for ii = 1:length(grand_struct.grand_struct.in_medfiltHilbert)
    subplot(4,4,ii)
    plot(grand_struct.grand_struct.in_medfiltHilbert(ii).sess_mean)
end 

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