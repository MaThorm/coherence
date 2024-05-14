clc
clear all
%%
load('attout_dataset.mat')
load('attin_dataset.mat')
load("V4_dataset.mat")

saving = true;

% Trial preprocessing parameters
bpfilt = false;
bpwidth = [1 150];
toi = [2.3 4.3]; % Time region of interest [2.3 4.3] is MC 2&3


% SSD parameters
fs = 1000; % Sampling frequency of the signal
th = 0.01; % residual variance threshhold   
lower = 52;
upper = 90;

% Hilbert parameters
filttype = "sgolay"; %either medfilt or sgolay
framelen = 31;
filtord = 1;

matpath = '/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files';

%%
% Trialselection 
[in_trials] = pre_processing_pip_trials(attin_dataset,bpfilt,bpwidth,toi)
[out_trials] = pre_processing_pip_trials(attout_dataset,bpfilt,bpwidth,toi)
[V4_trials] = pre_processing_pip_trials(V4_dataset,bpfilt,bpwidth,toi)

if saving == true
    LogicalStr = {'false', 'true'};
    save(fullfile(matpath,'trials',sprintf('in_trials_bpfilt%s_toi%.1f-%.1f.mat',LogicalStr{bpfilt+1},toi(1),toi(2))),"in_trials")
    save(fullfile(matpath,'trials',sprintf('out_trials_bpfilt%s_toi%.1f-%.1f.mat',LogicalStr{bpfilt+1},toi(1),toi(2))),"out_trials")
    save(fullfile(matpath,'trials',sprintf('V4_trials_bpfilt%s_toi%.1f-%.1f.mat',LogicalStr{bpfilt+1},toi(1),toi(2))),"V4_trials")
end 

%% Performing SSD
[in_SSD_trials,inc_in,testing_in] = do_SSD(in_trials,fs,th,lower,upper)
[out_SSD_trials,inc_out,testing_out] = do_SSD(out_trials,fs,th,lower,upper)
[V4_SSD_trials,inc_V4,testing0_V4] = do_SSD(V4_trials,fs,th,lower,upper)
%
if saving == true
    save(fullfile(matpath,'SSD',sprintf('in_ssd_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))),"in_SSD_trials")
    save(fullfile(matpath,'SSD',sprintf('out_ssd_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))),"out_SSD_trials")
    save(fullfile(matpath,'SSD',sprintf('V4_ssd_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))),"V4_SSD_trials")
    save(fullfile(matpath,'SSD',sprintf('inc_in_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))),"inc_in")
    save(fullfile(matpath,'SSD',sprintf('inc_out_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))),"inc_out")
    save(fullfile(matpath,'SSD',sprintf('inc_V4_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))),"inc_V4")
    save(fullfile(matpath,'SSD',sprintf('inc_in_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))),"inc_in")
    save(fullfile(matpath,'SSD',sprintf('inc_out_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))),"inc_out")
    save(fullfile(matpath,'SSD',sprintf('inc_V4_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))),"inc_V4")
    m = matfile(fullfile(matpath,'SSD_testing',sprintf('testing_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))),'Writable',true);
    m.in = testing_in;
    m.out = testing_out;
    m.V4 = testing_V4;
end 

%% Hilbert Angles, instantaneous frequency, filtered instantaneous
% frequency, instantaneous change

loading = true;
if loading == true
    load(fullfile(matpath,'SSD',sprintf('in_ssd_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))));
    load(fullfile(matpath,'SSD',sprintf('out_ssd_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))));
    load(fullfile(matpath,'SSD',sprintf('V4_ssd_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))));
    load(fullfile(matpath,'SSD',sprintf('inc_in_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))))
    load(fullfile(matpath,'SSD',sprintf('inc_out_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))))
    load(fullfile(matpath,'SSD',sprintf('inc_V4_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))))
end 
% While SSD_trials are still called in_trials
in_SSD_trials = in_trials;
out_SSD_trials = out_trials;
V4_SSD_trials = V4_trials;
% removning non-included trials 
removing = false;
if removing == true
    in_SSD_trials = remover(in_SSD_trials,inc_in)
    out_SSD_trials = remover(out_SSD_trials,inc_out)
    V4_SSD_trials = remover(V4_SSD_trials,inc_V4)
end 
    
% Hilberting: VERY IMPORTANT: Depending on the filter, the filtered instantaneous changes is either inst_freq (in the case of sgolay) or filt_inst_freq (in the case of medfilt)    
[hilb_angles.wrapped.in,hilb_angles.in,inst_freq.in,filt_inst_freq.in,inst_change.in]  = pre_processing_pip_hilb(in_SSD_trials,filttype,framelen,filtord);
[hilb_angles.wrapped.out,hilb_angles.out,inst_freq.out,filt_inst_freq.out,inst_change.out]  = pre_processing_pip_hilb(out_SSD_trials,filttype,framelen,filtord);
[hilb_angles.wrapped.V4,hilb_angles.V4,inst_freq.V4,filt_inst_freq.V4,inst_change.V4]  = pre_processing_pip_hilb(V4_SSD_trials,filttype,framelen,filtord);
timebar = -1.3:0.001:5;

if filttype == 'sgolay'     
    cur_in = inst_freq.in;
    cur_out = inst_freq.out;
    cur_V4 = inst_freq.V4
elseif filttype == 'medfilt'
    cur_in = filt_inst_freq.in;
    cur_out = filt_inst_freq.out;
    cur_V4 = filt_inst_freq.V4;
end 

%% Saving the file 
mt = [];
filename = sprintf("inst_freq_thr%.2f_comp%d_lower%d_upper%d.mat",th,10,lower,upper);
save(fullfile(matpath,'Inst_freq',filename),'mt','-v7.3')
m = matfile(fullfile(matpath,'Inst_freq',filename),'Writable',true);
m.inst_in = cur_in;
m.inst_out = cur_out;
m.inst_V4 = cur_V4;
m.angle_in = hilb_angles.wrapped.in;
m.angle_out = hilb_angles.wrapped.out;
m.angle_V4 = hilb_angles.wrapped.V4;

%%


% Calculating trial arrays
inst_in = (struct2matnan(cur_in,1));
inst_out = (struct2matnan(cur_out,1));
inst_V4 = (struct2matnan(cur_V4,1));
inst_inV4 = struct2matnan(cur_in,2); %Leaving come channels out since there ere often multiples
inst_outV4 = (struct2matnan(cur_out,2));
inst_inV4dif = inst_in - inst_inV4;
inst_outV4dif = inst_out - inst_outV4;
inst_inV4dif_mean = squeeze(mean(inst_inV4dif,2,'omitnan')); % So this should be x23 not x16
inst_outV4dif_mean = squeeze(mean(inst_outV4dif,2,'omitnan'));
inst_in_mean = squeeze(mean(inst_in,2,'omitnan'));
inst_out_mean = squeeze(mean(inst_out,2,'omitnan'));
inst_V4_mean = squeeze(mean(inst_V4,2,'omitnan'));
inst_in_meanmean = mean(inst_in_mean,1,'omitnan');
inst_out_meanmean = mean(inst_out_mean,1,'omitnan');
inst_V4_meanmean = mean(inst_V4_mean,1,'omitnan');
clear inst_in inst_out inst_V4 inst_inV4 inst_outV4 inst_inV4dif inst_outV4dif


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
x = timebar(toi(1)*1000:toi(2)*1000);
figure;
plot(x,inst_in_meanmean,'r');
hold on 
plot(x,inst_out_meanmean,'b');
plot(x,inst_V4_meanmean);
title([sprintf("Summary: bp %d - %d, filter",bpwidth(1),bpwidth(2)) filttype sprintf("length: %d, order: %d",framelen,filtord)])
legend('AttIn','AttOut','V4','AutoUpdate','off')
label = {'MS2','MS3','MS4'};
xl = xline([ 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
xlim([toi(1)-1.3 toi(2)-1.3])
%ylim([62 74])
hold off
xlabel('Time [s]')
ylabel('Frequency [Hz]')


%% plotting all together: per session 
figure;
for ii = 1:16
    subplot(2,1,1)
    x = timebar(toi(1)*1000:toi(2)*1000);
    plot(x,inst_in_mean(ii,:),'r');
    grid on 
    hold on 
    plot(x,inst_out_mean(ii,:),'b');
    plot(x,inst_V4_mean(ii,:));
    title([sprintf("Session %d: bp %d - %d, filter",ii, bpwidth(1),bpwidth(2)) filttype sprintf("length: %d, order: %d",framelen,filtord)])
    legend('AttIn','AttOut','V4','AutoUpdate','off')
    label = {'MS2','MS3','MS4'};
    xl = xline([ 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
    xlim([toi(1)-1.3 toi(2)-1.3])
    %ylim([62 74])
    hold off
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
    subplot(2,1,2)
    plot(x,inst_inV4dif_mean(ii,:),'r')
    grid on 
    hold on 
    plot(x,inst_outV4dif_mean(ii,:),'b')
    yline(mean(inst_outV4dif_mean(ii,:)),'b')    
    yline(mean(inst_inV4dif_mean(ii,:)),'r')
    title(sprintf('Session: %d',ii))
    hold off
    legend('In','Out')
    w = waitforbuttonpress;
    clf 
end 




%% Plotting inst_freq difference 
for ii = 1:size(inst_dif_mat,1)
    plot(inst_dif_mat(ii,:));
    w = waitforbuttonpress;
    clf
end 
%% Plotting mean of inst freq differences
inst_dif_mean = mean(inst_dif_mat,1)
figure
plot(inst_dif_mean)

%%
figure
sess = 1;
for ii = 1:1000
    subplot(2,1,1)
    plot(hilb_angles.wrapped.in(sess).trial{1, ii}(1,:));
    title(sprintf('Session %d, Trial %d',sess,ii))
    subplot(2,1,2)
    plot(inst_freq.in(sess).trial{1, ii}(1,:));
    w = waitforbuttonpress;
    clf;
end 

%% Quick derivative check
a = diff(summary.in.mean)
figure
plot(a)
ylim([-0.12 0.17])

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