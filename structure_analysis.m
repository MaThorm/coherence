% Script showing the general pipeline of my analysis 
clc 
clear all
dir = "/data/projects/V1V4coherence/02_analysis_max/git_repos/structure";
load(fullfile(dir,"attin_trials1150toi04.mat"));
%%
cfg = [];
ft_databrowser(cfg,in_trials(10))


%%
clc 
clear all
dir = "/data/projects/V1V4coherence/02_analysis_max/git_repos/structure";
base_foof = load(fullfile(dir,"osc_mean_bp1150_toi01.mat"));
task_foof = load(fullfile(dir,"osc_mean_bp1150_toi2.34.3.mat"));

%% Random session summary, plot of fooofed signal on log scale
num = 15;
hold on 
plot(log(task_foof.foof_summary.in_all.original(num).freq), log(task_foof.foof_summary.in_all.original(num).powspctrm),'k');
plot(log(task_foof.foof_summary.in_all.original(num).freq), log(task_foof.foof_summary.in_all.fractal(num).powspctrm));
plot(log(task_foof.foof_summary.in_all.original(num).freq), log(task_foof.foof_summary.in_all.osc_alt(num).powspctrm));
xlabel('log-freq'); ylabel('log-power'); grid on;
legend({'original','fractal','oscillatory = spectrum/fractal'},'location','southwest');
title('oscillatory = spectrum / fractal');
hold off
%% Random session summary, plot of foofed signal on non-log scale
num = 15
hold on 
plot(log(task_foof.foof_summary.V4_all.original(num).freq), log(task_foof.foof_summary.V4_all.original(num).powspctrm),'k');
plot(log(task_foof.foof_summary.V4_all.original(num).freq), log(task_foof.foof_summary.V4_all.fractal(num).powspctrm));
plot(log(task_foof.foof_summary.V4_all.original(num).freq), log(task_foof.foof_summary.V4_all.osc_alt(num).powspctrm));
xlabel('log-freq'); ylabel('log-power'); grid on;
legend({'original','fractal','oscillatory = spectrum/fractal'},'location','southwest');
title('oscillatory = spectrum / fractal');
hold off

%% Plotting mean oscillatory component over trials for all conditions 
figure('units','normalized','outerposition',[0 0 1 1]); hold on 
plot(log(task_foof.foof_summary.V4_all.original(1).freq), log(task_foof.foof_summary.in),'r');
plot(log(task_foof.foof_summary.V4_all.original(1).freq), log(task_foof.foof_summary.out),'b');
plot(log(task_foof.foof_summary.V4_all.original(1).freq), log(task_foof.foof_summary.V4));
x = str2double(xticklabels);
xticklabels(exp(x))
legend('Attin','Attout','V4','autoupdate','off')
yline([0])
hold off

%% Plotting mean oscillatory component over trials for all conditions 
figure('units','normalized','outerposition',[0 0 1 1]); hold on 
plot(task_foof.foof_summary.V4_all.original(1).freq, task_foof.foof_summary.in,'r');
plot(task_foof.foof_summary.V4_all.original(1).freq, task_foof.foof_summary.out,'b');
plot(task_foof.foof_summary.V4_all.original(1).freq, task_foof.foof_summary.V4);
legend('Attin','Attout','V4','autoupdate','off')
yline([0])
hold off

%% Loading the hilbert data 
clc 
clear all 
dir = "/data/projects/V1V4coherence/02_analysis_max/git_repos/structure";
load(fullfile(dir,"grand_structbp30100med20.mat"));
%% Session X, trial X Hilbert Angles
num = 1;
cfg = [];
ft_databrowser(cfg,grand_struct.add.in_hilbert(num))

%% Session X, trial X Inst Freq
ft_databrowser(cfg,grand_struct.add.inst_freq(num))

%% Session X, trial X median filter
ft_databrowser(cfg,grand_struct.in_medfiltHilbert(num))

%% Final summary: 
sel = 1:5000;
timebar = -1.3:0.001:5;
x = timebar(sel);
figure
plot(x,grand_struct.insummary.mean(sel),'r');
hold on 
plot(x,grand_struct.outsummary.mean(sel),'b');
plot(x,grand_struct.V4summary.mean(sel));
title("Summary of all conditions")
legend('AttIn','AttOut','V4','AutoUpdate','off')
label = {'Static', 'MS1','MS2','MS3','MS4'};
xl = xline([-0.5 0 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
hold off
xlabel('Time [s]')
ylabel('Frequency [Hz]')