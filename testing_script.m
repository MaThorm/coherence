clc
clear all
load('attout_dataset.mat')
load('attin_dataset.mat')
load("V4_dataset.mat")

% Trial preprocessing parameters
bpfilt = true;
bpwidth = [30 100];

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
[in1_trials] = pre_processing_pip_trials(attin_dataset,bpfilt,bpwidth,toi)
[out1_trials] = pre_processing_pip_trials(attout_dataset,bpfilt,bpwidth,toi)
[V41_trials] = pre_processing_pip_trials(V4_dataset,bpfilt,bpwidth,toi)

%Performing SSD
[in1ssd_trials] = do_SSD(in1_trials,fs,th)
[out1ssd_trials] = do_SSD(out1_trials,fs,th)
[V41ssd_trials] = do_SSD(V41_trials,fs,th)

bpwidth = [50 100];
% Trialselection 
[in2_trials] = pre_processing_pip_trials(attin_dataset,bpfilt,bpwidth,toi)
[out2_trials] = pre_processing_pip_trials(attout_dataset,bpfilt,bpwidth,toi)
[V42_trials] = pre_processing_pip_trials(V4_dataset,bpfilt,bpwidth,toi)

%Performing SSD
[in2ssd_trials] = do_SSD(in2_trials,fs,th)
[out2ssd_trials] = do_SSD(out2_trials,fs,th)
[V42ssd_trials] = do_SSD(V42_trials,fs,th)

parfor ii = 1:16
    in_1fft(ii) = do_fft(in1_trials(ii));
    inssd_1fft(ii) = do_fft(in1ssd_trials(ii));
    in_2fft(ii) = do_fft(in2_trials(ii));
    inssd_2fft(ii) = do_fft(in2ssd_trials(ii));
end 
% Hilbert Angles, instantaneous frequency, filtered instantaneous
% frequency, instantaneous change, summary struct
[hilb_angles.wrapped.in,hilb_angles.in,inst_freq.in,filt_inst_freq1.in,inst_change.in,summary.in]  = pre_processing_pip_hilb(in1ssd_trials,filttype,framelen,filtord);
[hilb_angles.wrapped.out,hilb_angles.out,inst_freq.out,filt_inst_freq1.out,inst_change.out,summary.out]  = pre_processing_pip_hilb(out1ssd_trials,filttype,framelen,filtord);
[hilb_angles.wrapped.V4,hilb_angles.V4,inst_freq.V4,filt_inst_freq1.V4,inst_change.V4,summary.V4]  = pre_processing_pip_hilb(V41ssd_trials,filttype,framelen,filtord);

[hilb_angles.wrapped.in,hilb_angles.in,inst_freq.in,filt_inst_freq2.in,inst_change.in,summary.in]  = pre_processing_pip_hilb(in2ssd_trials,filttype,framelen,filtord);
[hilb_angles.wrapped.out,hilb_angles.out,inst_freq.out,filt_inst_freq2.out,inst_change.out,summary.out]  = pre_processing_pip_hilb(out2ssd_trials,filttype,framelen,filtord);
[hilb_angles.wrapped.V4,hilb_angles.V4,inst_freq.V4,filt_inst_freq2.V4,inst_change.V4,summary.V4]  = pre_processing_pip_hilb(V42ssd_trials,filttype,framelen,filtord);
timebar = -1.3:0.001:5;

% Calculating trial arrays
inst1_in = (struct2matnan(filt_inst_freq1.in,1));
inst1_out = (struct2matnan(filt_inst_freq1.out,1));
inst1_V4 = (struct2matnan(filt_inst_freq1.V4,1));
inst1_inV4 = struct2matnan(filt_inst_freq1.in,2);
inst1_outV4 = (struct2matnan(filt_inst_freq1.out,2));
inst1_inV4dif = inst1_in - inst1_inV4;
inst1_outV4dif = inst1_out - inst1_outV4;
inst1_inV4dif_mean = squeeze(mean(inst1_inV4dif,2,'omitnan'));
inst1_outV4dif_mean = squeeze(mean(inst1_outV4dif,2,'omitnan'));

inst1_in_mean = squeeze(mean(inst1_in,2,'omitnan'));
inst1_out_mean = squeeze(mean(inst1_out,2,'omitnan'));
inst1_V4_mean = squeeze(mean(inst1_V4,2,'omitnan'));
inst1_in_meanmean = mean(inst1_in_mean,1,'omitnan');
inst1_out_meanmean = mean(inst1_out_mean,1,'omitnan');
inst1_V4_meanmean = mean(inst1_V4_mean,1,'omitnan');

inst2_in = (struct2matnan(filt_inst_freq2.in,1));
inst2_out = (struct2matnan(filt_inst_freq2.out,1));
inst2_V4 = (struct2matnan(filt_inst_freq2.V4,1));
inst2_inV4 = struct2matnan(filt_inst_freq2.in,2);
inst2_outV4 = (struct2matnan(filt_inst_freq2.out,2));
inst2_inV4dif = inst2_in - inst2_inV4;
inst2_outV4dif = inst2_out - inst2_outV4;
inst2_inV4dif_mean = squeeze(mean(inst2_inV4dif,2,'omitnan'));
inst2_outV4dif_mean = squeeze(mean(inst2_outV4dif,2,'omitnan'));

inst2_in_mean = squeeze(mean(inst2_in,2,'omitnan'));
inst2_out_mean = squeeze(mean(inst2_out,2,'omitnan'));
inst2_V4_mean = squeeze(mean(inst2_V4,2,'omitnan'));
inst2_in_meanmean = mean(inst2_in_mean,1,'omitnan');
inst2_out_meanmean = mean(inst2_out_mean,1,'omitnan');
inst2_V4_meanmean = mean(inst2_V4_mean,1,'omitnan');
%% Plotting FFT results
figure;
for ii = 1:16
    subplot(2,1,1)
    plot(in_1fft(ii).freq,in_1fft(ii).powspctrm,'r');
    hold on 
    plot(in_2fft(ii).freq,in_2fft(ii).powspctrm,'b');
    hold off 
    legend({'Wide','Small'})
    title('FFT')
    subplot(2,1,2)
    plot(inssd_1fft(ii).freq,inssd_1fft(ii).powspctrm,'r');    
    hold on 
    plot(inssd_2fft(ii).freq,inssd_2fft(ii).powspctrm,'b');
    hold off
    legend({'Wide','Small'})
    title('FFT SSD')
    w = waitforbuttonpress;
    clf;
end 

 %% plotting all together: per session 
figure;
for ii = 1:16
    subplot(2,1,1)
    x = timebar(toi(1)*1000:toi(2)*1000);
    plot(x,inst1_in_mean(ii,:),'r');
    grid on 
    hold on 
    plot(x,inst1_out_mean(ii,:),'b');
    plot(x,inst1_V4_mean(ii,:));
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
    x = timebar(toi(1)*1000:toi(2)*1000);
    plot(x,inst2_in_mean(ii,:),'r');
    grid on 
    hold on 
    plot(x,inst2_out_mean(ii,:),'b');
    plot(x,inst2_V4_mean(ii,:));
    title([sprintf("Session %d: bp %d - %d, filter",ii, bpwidth(1),bpwidth(2)) filttype sprintf("length: %d, order: %d",framelen,filtord)])
    legend('AttIn','AttOut','V4','AutoUpdate','off')
    label = {'MS2','MS3','MS4'};
    xl = xline([ 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
    xlim([toi(1)-1.3 toi(2)-1.3])
    %ylim([62 74])
    hold off
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
    w = waitforbuttonpress;
    clf 
end   





%% Using the fieldtrip averaging
cfg = [];
for ii = 1:16
    in_avg(ii) = ft_timelockanalysis(cfg,filt_inst_freq.out(ii))
end 
%
for ii = 1:16
    in_avg(ii).label{1,1} = 'test';
end 
cfg = [];
cfg.channel = 'test';
in_avg_avg = ft_timelockgrandaverage(cfg,in_avg(1),in_avg(2),in_avg(3),in_avg(4),in_avg(5),in_avg(6),in_avg(7),in_avg(8) ...
    ,in_avg(9),in_avg(10),in_avg(11),in_avg(12),in_avg(13),in_avg(14),in_avg(15),in_avg(16));

