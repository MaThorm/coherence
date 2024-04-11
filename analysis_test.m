clc
clear all
load('attout_dataset.mat')
load('attin_dataset.mat')
load("V4_dataset.mat")

% Trial preprocessing parameters
bpfilt = true;
bpwidth = [30 100];

% SSD parameters
toi = [2.3 4.3];
fs = 1000;
th = 0.01;

% Hilbert parameters
filttype = "sgolay"; %either medfilt or sgolay
framelen = 21;
filtord = 1;
ord = [1 3 5];
len = [11 21 31]
matpath = '/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files';
paths = '/home/mthormann@brain.uni-bremen.de/V1V4coherence/03_results_max/inst_freq_loop/sgolay_filts/'
for i_ord = 1:length(ord)
    filtord = ord(i_ord)
    for i_len = 1:length(len)
        framelen = len(i_len)
        [in_trials] = pre_processing_pip_trials(attin_dataset,bpfilt,bpwidth,toi)
        [out_trials] = pre_processing_pip_trials(attout_dataset,bpfilt,bpwidth,toi)
        [V4_trials] = pre_processing_pip_trials(V4_dataset,bpfilt,bpwidth,toi)
        [in_trials] = do_SSD(in_trials,fs,th)
        [out_trials] = do_SSD(out_trials,fs,th)
        [V4_trials] = do_SSD(V4_trials,fs,th)
        [grand_struct,angles,inst_freq] = pre_processing_pip_hilb(in_trials,out_trials,V4_trials,filttype,framelen,filtord)
        
        attin_inst = grand_struct.in_medfiltHilbert;
        attout_inst = grand_struct.out_medfiltHilbert;
        insummary = grand_struct.insummary;
        outsummary = grand_struct.outsummary;
        V4_inst = grand_struct.V4_medfiltHilbert;
        V4summary = grand_struct.V4summary;
        timebar = -1.3:0.001:5;
        
        %plotting all together, but only MS 2 & 3 
        sel = 1:length(outsummary.mean);
        x = timebar(sel);
        figure('units','normalized','outerposition',[0 0 1 1],'visible','off')
        plot(x,insummary.mean(sel),'r');
        hold on 
        plot(x,outsummary.mean(sel),'b');
        plot(x,V4summary.mean(sel));
        title(sprintf("Summary: sgolay_ord%d_len%d",ord(i_ord),len(i_len)))
        legend('AttIn','AttOut','V4','AutoUpdate','off')
        label = {'MS2','MS3','MS4'};
        ylim([60 75])
        %xl = xline([ 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
        hold off
        xlabel('Time [s]')
        ylabel('Frequency [Hz]')
        saveas(gcf,fullfile(paths,sprintf("sgolay_ord%d_len%d.jpg",ord(i_ord),len(i_len))))
    end 
end 
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
sel = 1:length(outsummary.mean);
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
sel = 1:length(insummary.mean);
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
sel = 1:length(outsummary.mean);
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
sel = 1:length(V4summary.mean);
x = timebar(sel)
y = V4summary.mean(sel);
sd = V4summary.std(sel);
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
sel = 1:length(outsummary.mean);
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



%% Quick derivative check
a = diff(insummary.mean(sel))
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
 
%% Testing script of ssd 
% EXAMPLE:
% Creation of some simple timeserie
y  = sin(2*pi*5*(0:999)/1000);
y2 = 0.1*sin(2*pi*15*(0:999)/1000);
y3 = y+y2;
y3(500:999) = y3(500:999)+sin(2*pi*75*(500:999)/1000);

x1 = sin(2*pi*5*(0:999)/1000);
x2 = [zeros(1,500) sin(2*pi*75*(501:1000)/1000)];
x3 = 0.1*sin(2*pi*15*(0:999)/1000);

v  = y3;
% Sampling frequency 1000 and threshold of 0.005
SSDcomponents = SSD(v,1000,0.01);