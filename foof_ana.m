clc 
clear 
load('attin_dataset.mat')
load('attout_dataset.mat')
load('V4_dataset.mat')
%
matpath = "/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files";
figurepath = "/home/mthormann@brain.uni-bremen.de/V1V4coherence/03_results_max/fooof/"
bpwidth = [1 150];
toilim = [0 1];
%%
parfor ii = 1:length(attin_dataset)
    in_trials(ii) = do_trialselection(attin_dataset(ii).path,attin_dataset(ii).file,attin_dataset(ii).chan,attin_dataset(ii).stimno,bpwidth,toilim);
end 
parfor ii = 1:length(attout_dataset)
    out_trials(ii) = do_trialselection(attout_dataset(ii).path,attout_dataset(ii).file,attout_dataset(ii).chan,attout_dataset(ii).stimno,bpwidth,toilim);
end 
parfor ii = 1:length(V4_dataset)
    V4_trials(ii) = do_trialselection(V4_dataset(ii).path,V4_dataset(ii).file,V4_dataset(ii).chan,V4_dataset(ii).stimno,bpwidth,toilim);
end 

%%
save(fullfile(matpath,sprintf("attin_trials%d%dtoi%d%d.mat",bpwidth(1),bpwidth(2),toilim(1),toilim(2))),"in_trials");
save(fullfile(matpath,sprintf("attout_trials%d%dtoi%d%d.mat",bpwidth(1),bpwidth(2),toilim(1),toilim(2))),"out_trials");
save(fullfile(matpath,sprintf("V4_trials%d%dtoi%d%d.mat",bpwidth(1),bpwidth(2),toilim(1),toilim(2))),"V4_trials");


%%
clearvars -except toilim bpwidth matpath figurepath
attin = load(fullfile(matpath,sprintf("attin_trials%d%dtoi%d%d.mat",bpwidth(1),bpwidth(2),toilim(1),toilim(2))));
attout = load(fullfile(matpath,sprintf("attout_trials%d%dtoi%d%d.mat",bpwidth(1),bpwidth(2),toilim(1),toilim(2))));
V4 = load(fullfile(matpath,sprintf("V4_trials%d%dtoi%d%d.mat",bpwidth(1),bpwidth(2),toilim(1),toilim(2))));
in_trials = attin.in_trials;
out_trials = attout.out_trials;
V4_trials = V4.V4_trials;

%% Foofing all the things
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
save(fullfile(matpath,sprintf("osc_mean_bp%d%d_toi%d%d.mat",bpwidth(1),bpwidth(2),toilim(1),toilim(2))),"foof_summary")



%% Plottig a random session with fractal, original & oscillatory component
num = 2;
hold on 
plot(log(in_foof.original(1).freq), log(in_foof.original(num).powspctrm),'k');
plot(log(in_foof.original(1).freq), log(in_foof.fractal(num).powspctrm));
plot(log(in_foof.original(1).freq), log(in_foof.osc_alt(num).powspctrm));
xlabel('log-freq'); ylabel('log-power'); grid on;
legend({'original','fractal','oscillatory = spectrum/fractal'},'location','southwest');
title('oscillatory = spectrum / fractal');
hold off
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
saveas(gcf,fullfile(figurepath,sprintf('foofmean_bp%d%d_toi%d%d.jpg',bpwidth(1),bpwidth(2),toilim(1),toilim(2))))

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
