%% Creating struct for trial databrowser
clc 
clear all
dir = "/data/projects/V1V4coherence/02_analysis_max/git_repos/structure";
load("attin_dataset.mat");
load("attout_dataset.mat");
load("V4_dataset.mat");
bpwidth = [1 150];
toi = [-100 100];
[in_trials,out_trials,V4_trials] = pre_processing_pip_trials(attin_dataset,attout_dataset,V4_dataset,bpwidth,toi)
trial_struct.in = in_trials;
trial_struct.out = out_trials;
trial_struct.V4 = V4_trials;
save(fullfile(dir,"trial_struct.mat"),"trial_struct");
%% Creating a fooof struct
clc
clear all
load('attout_dataset.mat')
load('attin_dataset.mat')
load("V4_dataset.mat")
bpwidth = [1 300];
toi = [2.3 4.3];
dir = "/data/projects/V1V4coherence/02_analysis_max/git_repos/structure";
[in_trials,out_trials,V4_trials] = pre_processing_pip_trials(attin_dataset,attout_dataset,V4_dataset,bpwidth,toi)


% Foofing all the things
for ii = 1:length(in_trials)
    [in_foof.fractal(ii), in_foof.original(ii), in_foof.osc(ii) in_foof.osc_alt(ii)] = do_foof(in_trials(ii))
    [out_foof.fractal(ii), out_foof.original(ii), out_foof.osc(ii) out_foof.osc_alt(ii)] = do_foof(out_trials(ii))
    [V4_foof.fractal(ii), V4_foof.original(ii), V4_foof.osc(ii) V4_foof.osc_alt(ii)] = do_foof(V4_trials(ii))
end 
% Taking mean over trials in
cells = struct2cell(in_foof.osc_alt);
cells = squeeze(cells(4,:,:));
mat = cell2matnan(cells');
in_foof.osc_mean = mean(mat,1);
% mean over trials out
cells = struct2cell(out_foof.osc_alt);
cells = squeeze(cells(4,:,:));
mat = cell2matnan(cells');
out_foof.osc_mean = mean(mat,1);
% mean over trials V4
cells = struct2cell(V4_foof.osc_alt);
cells = squeeze(cells(4,:,:));
mat = cell2matnan(cells');
V4_foof.osc_mean = mean(mat,1);
foof_summary.V4 = V4_foof.osc_mean;
foof_summary.in = in_foof.osc_mean;
foof_summary.out = out_foof.osc_mean;
foof_summary.in_all.fractal = in_foof.fractal;
foof_summary.in_all.original = in_foof.original;
foof_summary.in_all.osc_alt = in_foof.osc_alt;
foof_summary.out_all.fractal = out_foof.fractal;
foof_summary.out_all.original = out_foof.original;
foof_summary.out_all.osc_alt = out_foof.osc_alt;
foof_summary.V4_all.fractal = V4_foof.fractal;
foof_summary.V4_all.original = V4_foof.original;
foof_summary.V4_all.osc_alt = V4_foof.osc_alt;
save(fullfile(dir,"foof_struct.mat"),"foof_summary");


%% Creating a Hilbert Struct
clc
clear all
load('attout_dataset.mat')
load('attin_dataset.mat')
load("V4_dataset.mat")
bpwidth = [52 78];
toi = [-100 100];
medianfiltord = 10;
dir = "/data/projects/V1V4coherence/02_analysis_max/git_repos/structure";
save_hilbert = true; 

[in_trials,out_trials,V4_trials] = pre_processing_pip_trials(attin_dataset,attout_dataset,V4_dataset,bpwidth,toi)
[grand_struct,angles,inst_freq] = pre_processing_pip_hilb(in_trials,out_trials,V4_trials,medianfiltord,save_hilbert)

save(fullfile(dir,"Hilbert_struct.mat"),"grand_struct");

%% Showing scripts and things
% First part is just databrowser of cut up sessions
clc
clear all
dir = "/data/projects/V1V4coherence/02_analysis_max/git_repos/structure";
%load(fullfile(dir,"attout_trials1150toi04.mat"));
load(fullfile(dir,"trial_struct.mat"))
cfg = [];
ft_databrowser(cfg,trial_struct.V4(12))


%% Second part: fooof
clc 
clear all
dir = "/data/projects/V1V4coherence/02_analysis_max/git_repos/structure";
%base_foof = load(fullfile(dir,"osc_mean_bp1150_toi01.mat"));
%task_foof = load(fullfile(dir,"osc_mean_bp1150_toi2.34.3.mat"));
task_foof = load(fullfile(dir,"foof_struct.mat"));


%% Random session summary, plot of fooofed signal on log scale
task = task_foof.foof_summary.out_all;
for num = 1:16
plot(log(task_foof.foof_summary.in_all.original(num).freq), log(task.original(num).powspctrm),'k');
hold on 
plot(log(task_foof.foof_summary.in_all.original(num).freq), log(task.fractal(num).powspctrm));
plot(log(task_foof.foof_summary.in_all.original(num).freq), log(task.osc_alt(num).powspctrm));
xlabel('log-freq'); ylabel('log-power'); grid on;
legend({'original','fractal','oscillatory = spectrum/fractal'},'location','southwest');
title(sprintf('oscillatory = spectrum / fractal, sess: %d',num));
hold off
w = waitforbuttonpress;
clf
end 
%% Random session summary, plot of foofed signal on non-log scale
task = task_foof.foof_summary.in_all;
for num = 1:16
plot(log(task_foof.foof_summary.V4_all.original(num).freq), log(task.original(num).powspctrm),'k');
hold on 
plot(log(task_foof.foof_summary.V4_all.original(num).freq), log(task.fractal(num).powspctrm));
plot(log(task_foof.foof_summary.V4_all.original(num).freq), log(task.osc_alt(num).powspctrm));
xlabel('log-freq'); ylabel('log-power'); grid on;
legend({'original','fractal','oscillatory = spectrum/fractal'},'location','southwest');
title(sprintf('oscillatory = spectrum / fractal, sess: %d',num));
hold off
w = waitforbuttonpress;
clf
end 

%% Plotting mean oscillatory component over trials for all conditions 
figure('units','normalized','outerposition',[0 0 1 1]); hold on 
plot(log(task_foof.foof_summary.V4_all.original(1).freq), log(task_foof.foof_summary.in),'r');
plot(log(task_foof.foof_summary.V4_all.original(1).freq), log(task_foof.foof_summary.out),'b');
plot(log(task_foof.foof_summary.V4_all.original(1).freq), log(task_foof.foof_summary.V4));
x = str2double(xticklabels);
xticklabels(exp(x))
legend('Attin','Attout','V4','autoupdate','off')
xlabel('Frequency [Hz]')
ylabel('Power')
title('Mean oscillatory component over all sessions using the fooof method, log scaled')
yline(0)
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
load(fullfile(dir,"Hilbert_struct.mat"));
%% Session X, trial X Hilbert Angles
num = 1;
cfg = [];
ft_databrowser(cfg,grand_struct.in_hilbertData(num))

%% Session X, trial X Inst Freq
ft_databrowser(cfg,grand_struct.in_diffHilbertData(num))

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
