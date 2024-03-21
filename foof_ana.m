clc
clear all
load('attout_dataset.mat')
load('attin_dataset.mat')
load("V4_dataset.mat")
bpwidth = [1 150];
toi = [2.3 4.3];
matpath = '/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files'
[in_trials,out_trials,V4_trials] = pre_processing_pip_trials(attin_dataset,attout_dataset,V4_dataset,bpwidth,toi)


% Foofing all the things
for ii = 1:length(in_trials)
    [in_foof.fractal(ii), in_foof.original(ii), in_foof.osc(ii) in_foof.osc_alt(ii)] = do_foof(in_trials(ii))
    [out_foof.fractal(ii), out_foof.original(ii), out_foof.osc(ii) out_foof.osc_alt(ii)] = do_foof(out_trials(ii))
    [V4_foof.fractal(ii), V4_foof.original(ii), V4_foof.osc(ii) V4_foof.osc_alt(ii)] = do_foof(V4_trials(ii))
end 
% Taking mean over trials in
cells = struct2cell(in_foof.osc_alt)
cells = squeeze(cells(4,:,:))
mat = cell2matnan(cells')
in_foof.osc_mean = mean(mat,1)
% mean over trials out
cells = struct2cell(out_foof.osc_alt)
cells = squeeze(cells(4,:,:))
mat = cell2matnan(cells')
out_foof.osc_mean = mean(mat,1)
% mean over trials V4
cells = struct2cell(V4_foof.osc_alt)
cells = squeeze(cells(4,:,:))
mat = cell2matnan(cells')
V4_foof.osc_mean = mean(mat,1)
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



%% Plottig a random session with fractal, original & oscillatory component
num = 4;
for num = 1:16
hold on 
plot(log(in_foof.original(1).freq), log(in_foof.original(num).powspctrm),'k');
plot(log(in_foof.original(1).freq), log(in_foof.fractal(num).powspctrm));
plot(log(in_foof.original(1).freq), log(in_foof.osc_alt(num).powspctrm));
xlabel('log-freq'); ylabel('log-power'); grid on;
legend({'original','fractal','oscillatory = spectrum/fractal'},'location','southwest');
title('oscillatory = spectrum / fractal');
hold off
w = waitforbuttonpress;
clf
end 
%% Plotting mean oscillatory component over trials for all conditions 
figure('units','normalized','outerposition',[0 0 1 1]); hold on 
plot(log(in_foof.original(1).freq), log(in_foof.osc_mean),'r');
plot(log(in_foof.original(1).freq), log(out_foof.osc_mean),'b');
plot(log(in_foof.original(1).freq), log(V4_foof.osc_mean));
x = str2double(xticklabels);
xticklabels(exp(x))
legend('Attin','Attout','V4','autoupdate','off')
yline([0])
hold off

%% Plotting mean original component of a random session 
figure('units','normalized','outerposition',[0 0 1 1]); hold on 
plot(log(in_foof.original(1).freq), in_foof.original(1).powspctrm,'r');
plot(log(in_foof.original(1).freq), out_foof.original(1).powspctrm,'b');
plot(log(in_foof.original(1).freq), V4_foof.original(1).powspctrm);
x = str2double(xticklabels);
xticklabels(exp(x))
legend('Attin','Attout','V4','autoupdate','off')
yline([0])
hold off

%% Comparing subtracting 1 - 3 from 0 - 1 actually other way around
part2 = load(fullfile(matpath,sprintf("osc_mean_bp%d%d_toi%d%d.mat",bpwidth(1),bpwidth(2),2.3,4.3)));
part1 = load(fullfile(matpath,sprintf("osc_mean_bp%d%d_toi%d%d.mat",bpwidth(1),bpwidth(2),0,1)));
new_V4 = part2.foof_summary.V4 - part1.foof_summary.V4 ;
new_in = part2.foof_summary.in - part1.foof_summary.in ;
new_out = part2.foof_summary.out - part1.foof_summary.out;
figure('units','normalized','outerposition',[0 0 1 1]); hold on 
plot(log(in_foof.original(1).freq),log(new_in),'r');
plot(log(in_foof.original(1).freq),log(new_out),'b');
plot(log(in_foof.original(1).freq),log(new_V4));
x = str2double(xticklabels);
xticklabels(exp(x))
legend('Attin','Attout','V4','autoupdate','off')
yline([0])
hold off

%% Plotting all means in big pic
