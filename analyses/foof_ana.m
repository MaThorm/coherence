%clc
%clear all
load('attout_dataset.mat')
load('attin_dataset.mat')
load("V4_dataset.mat")


% Trial preprocessing parameters
bpfilt = false;
bpwidth = [1 150];
toi = [2.3 4.3];
 
matpath = '/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files/fooof';
% Trialselection 
[in_trials_test] = pre_processing_pip_trials(attin_dataset,bpfilt,bpwidth,toi)
[out_trials_test] = pre_processing_pip_trials(attout_dataset,bpfilt,bpwidth,toi)
[V4_trials_test] = pre_processing_pip_trials(V4_dataset,bpfilt,bpwidth,toi)

% Foofing all the things
for ii = 1:length(in_trials_test)
    [in_foof.fractal(ii), in_foof.original(ii), in_foof.osc(ii) in_foof.osc_alt(ii)] = do_foof(in_trials_test(ii))
    [out_foof.fractal(ii), out_foof.original(ii), out_foof.osc(ii) out_foof.osc_alt(ii)] = do_foof(out_trials_test(ii))
    [V4_foof.fractal(ii), V4_foof.original(ii), V4_foof.osc(ii) V4_foof.osc_alt(ii)] = do_foof(V4_trials_test(ii))
end 
% Taking mean over trials in
cells = struct2cell(in_foof.osc_alt);
cells = squeeze(cells(4,:,:));
mat = cell2matnan(cells',1);
in_foof.osc_mean = mean(mat,1);
% mean over trials out
cells = struct2cell(out_foof.osc_alt);
cells = squeeze(cells(4,:,:));
mat = cell2matnan(cells',1);
out_foof.osc_mean = mean(mat,1);
% mean over trials V4
cells = struct2cell(V4_foof.osc_alt);
cells = squeeze(cells(4,:,:));
mat = cell2matnan(cells',1);
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


% Calculating the mean peak over all conditions
grand_foof = mean([in_foof.osc_mean; out_foof.osc_mean; V4_foof.osc_mean],1);
foof_summary.grand_foof = grand_foof

matpath = "/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files/fooof/"
save(fullfile(matpath,sprintf('fooof_summary_toi%.1f-%.1f.mat',toi(1),toi(2))),'foof_summary')

%% Plottig a random session with fractal, original & oscillatory component
figpath = "/data/projects/V1V4coherence/03_results_max/fooof/";
saving = true;
for num = 1:16
    f = figure;
    f.Units = 'normalized';
    f.Position = [0 0 0.5 0.5];
    plot(log(foof_summary.in_all.original(1).freq), log(foof_summary.in_all.osc_alt(num).powspctrm));
    hold on 
    plot(log(foof_summary.out_all.original(1).freq), log(foof_summary.out_all.osc_alt(num).powspctrm));
    plot(log(foof_summary.V4_all.original(1).freq), log(foof_summary.V4_all.osc_alt(num).powspctrm));
    x = str2double(xticklabels);
    xticklabels(exp(x))
    xlabel('freq [hz]'); ylabel('log-power'); grid on;
    legend({'in','out','V4'},'location','southwest');
    title(sprintf('Fooof results sess %d',num));
    yline(0)
    hold off
    if saving == true
        saveas(gcf,fullfile(figpath,sprintf('Fooof_Sess%i.jpg',num)));
        saveas(gcf,fullfile(figpath,sprintf('Fooof_Sess%i.fig',num)));
    else 
        w = waitforbuttonpress;
        clf
    end 
    close all 
end 
%% Plotting mean oscillatory component over trials for all conditions 
figpath = "/data/projects/V1V4coherence/03_results_max/fooof/";
figure('units','normalized','outerposition',[0 0 1 1]); hold on 
plot(log(foof_summary.in_all.osc_alt(1).freq), log(foof_summary.in),'r');
plot(log(foof_summary.in_all.osc_alt(1).freq), log(foof_summary.out),'b');
plot(log(foof_summary.in_all.osc_alt(1).freq), log(foof_summary.V4));
x = str2double(xticklabels);
xticklabels(exp(x))
legend('V1A','V1N','V4','autoupdate','off')
xlabel('Frequency [hz]')
ylabel('Log Power')
yline([0])
title('Fooof Summary')
hold off
%saveas(gcf,fullfile(figpath,'Fooof_Summary.jpg'));
%saveas(gcf,fullfile(figpath,'Fooof_Summary.fig'));

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

%% Plotting grand foof
figpath = "/data/projects/V1V4coherence/03_results_max/fooof/"
figure('Units','normalized','Position',[0 0 0.5 0.5])
plot(log(foof_summary.in_all.original(1).freq),grand_foof);
hold on
peak = max(grand_foof(:,100:end));
peak_ar = find(grand_foof(:,100:end) >= peak/2);
pmin = min(peak_ar)+100;
pmax = max(peak_ar)+100;
x = str2double(xticklabels);
xline(log(in_foof.original(1).freq(pmin)))
xline(log(in_foof.original(1).freq(pmax)))
text1 = num2str(in_foof.original(1).freq(pmin));
text2 = num2str(in_foof.original(1).freq(pmax));
text(4,5,text1)
text(4.7,5,text2)
xticklabels(exp(x))
title('Combined fooof oscillatory components, lines at FWHM of the peak')
xlabel('Frequency')
ylabel('Power')
yline([0])
hold off
saveas(gcf,fullfile(figpath,'Combined_Fooof_FWHM.jpg'))
saveas(gcf,fullfile(figpath,'Combined_Fooof_FWHM.fig'))
