clc
clear 
%% Parameter Settings


saving = true;

% Trial preprocessing parameters
bpfilt = false;
bpwidth = [50 80];
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


%%
% loading inst freq struct 
filename = sprintf("inst_freq_thr%.2f_comp%d_lower%d_upper%d.mat",th,10,lower,upper);
m = matfile(fullfile(matpath,'Inst_freq',filename),'Writable',true);
filt_inst_freq = m.inst;

%%
cur_in = elongate(filt_inst_freq.in);
cur_out = elongate(filt_inst_freq.out);
cur_V4 = elongate_V4(filt_inst_freq.V4)
v4_long = [1 1 2 3 2 3 4 5 6 7 8 9 9 10 11 12 13 14 15 13 14 15 16]
% for ii = 1:length(v4_long)
%     cur_V4(ii) = filt_inst_freq.V4(v4_long(ii))
% end 
elongated = true;
[inst_mean.in,grand_inst.in] = do_grand_avg(cur_in,elongated);
[inst_mean.out,grand_inst.out] = do_grand_avg(cur_out,elongated);
[inst_mean.V4,grand_inst.V4] = do_grand_avg(cur_V4);

%% Plotting all the noncut inst. freqs together 
x = timebar(toi(1)*1000:toi(2)*1000);
f = figure;
f.Units = 'normalized'
f.Position = [0 0 0.5 0.5]
plot(x,grand_inst.in.avg(1,:),'r');
hold on 
plot(x,grand_inst.out.avg(1,:),'b');
plot(x,grand_inst.in.avg(2,:));
plot(x,grand_inst.out.avg(2,:));

title([sprintf("Summary: lower: %d - upper: %d",lower,upper) filttype sprintf("length: %d, order: %d",framelen,filtord)])
legend('AttIn','AttOut','V4a','V4n','AutoUpdate','off')
label = {'MS2','MS3','MS4'};
xl = xline([ 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
xlim([toi(1)-1.3 toi(2)-1.3])
%ylim([62 74])
hold off
xlabel('Time [s]')
ylabel('Frequency [Hz]') 
foldername = fullfile(figpath,'inst_freq_23',sprintf('%d-%d',lower,upper))
if ~exist(foldername,'dir')
    mkdir(foldername)
end 
saveas(gcf,fullfile(foldername,sprintf("inst_freq_summary_low%d_high%d.jpg",lower,upper)));
saveas(gcf,fullfile(foldername,sprintf("inst_freq_summary_low%d_high%d.fig",lower,upper)));

%% plotting all noncut together: per session 
f = figure;
f.Units = ['normalized']
f.Position = [0 0 0.5 0.5]
for ii = 1:length(cur_in)
    x = timebar(toi(1)*1000:toi(2)*1000);
    plot(x,inst_mean.in(ii).avg(1,:),'r');
    grid on 
    hold on 
    plot(x,inst_mean.out(ii).avg(1,:),'b');
    plot(x,inst_mean.V4(ii).avg);

    title([sprintf("Pair %d: lower: %d - upper: %d, filter",ii, lower,upper) filttype sprintf("length: %d, order: %d",framelen,filtord)])
    legend('AttIn','AttOut','V4','AutoUpdate','off')
    label = {'MS2','MS3','MS4'};
    xl = xline([ 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
    xlim([toi(1)-1.3 toi(2)-1.3])
    %ylim([62 74])
    hold off
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
    foldername = fullfile(figpath,'inst_freq_23',sprintf('%d-%d',lower,upper))
    if ~exist(foldername,'dir')
        mkdir(foldername)
    end 
    saveas(gcf,fullfile(foldername,sprintf("inst_freq_pair%d_low%d_high%d.jpg",ii,lower,upper)));
    saveas(gcf,fullfile(foldername,sprintf("inst_freq_pair%d_low%d_high%d.fig",ii,lower,upper)));
    %w = waitforbuttonpress;
    clf 
end 
close all

%%
f = figure;
f.Units = ['normalized']
f.Position = [0 0 0.5 0.5]
for ii = 1:23
    x = timebar(toi(1)*1000:toi(2)*1000);
    plot(x,inst_mean.in(ii).avg(2,:));
    grid on 
    hold on 
    plot(x,inst_mean.out(ii).avg(2,:));
    plot(x,inst_mean.V4(ii).avg)
    title([sprintf("Session %d: lower: %d - upper: %d, filter",ii, lower,upper) filttype sprintf("length: %d, order: %d",framelen,filtord)])
    legend('V41','V42','V4','AutoUpdate','off')
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