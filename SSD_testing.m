clear all
matpath = "/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files";

m = matfile(fullfile(matpath,"wavelet/trial_SSDwavelet52-90_red.mat"));
wavSSD.in = m.in;
wavSSD.out = m.out;
wavSSD.V4 = m.V4;

m = matfile(fullfile(matpath,"SSD_testing/testing_th0.01_bounds52-90_ncomp10_toi2.3-4.3.mat"));
testing.in = m.in;
testing.out = m.out;
testing.V4 = m.V4;

m = matfile(fullfile(matpath,"Inst_freq/inst_freq_thr0.01_comp10_lower52_upper90.mat"));
inst.in = m.inst_in;
inst.out = m.inst_out;
inst.V4 = m.inst_V4;


%% Checking average wavelet power of an SSD component in the gamma range
% .freq values 
% 52Hz == 24
% 95Hz == 30  
% 106hz = 31
low = 24;
up = 30;
sess = 1;
trial = 1;
border = 10;
cur_waveSSD = wavSSD.in(sess);  % either in out or V4
cur_inst = inst.in(sess); % either in out or V4
cur_trial = squeeze(cur_waveSSD.powspctrm(trial,:,:,:));

%% PLotting SSD wavelet and average gamma activation 
f = figure;
f.Units = 'normalized';
f.Position = [0 0 0.5 0.5];
for trial = 1:length(cur_waveSSD.trialinfo)
cur_trial = squeeze(cur_waveSSD.powspctrm(trial,:,:,:));
    for ii = 1:length(cur_trial)
        avg_gamma(ii) = mean(cur_trial(low:up,ii),1,'omitnan');
    end 
nans = ~isnan(avg_gamma);
avg_gamma_nonan = avg_gamma(nans); % removing nans for assessing ratio
ratio = length(avg_gamma_nonan(avg_gamma_nonan >= border))/length(avg_gamma_nonan);
subplot(2,1,1)
plot(avg_gamma)
title(sprintf('Session %d, trial %d, ratio %.2f' ,sess,trial, ratio));
%yline(border);
subplot(2,1,2)
ft_plot_matrix(cur_waveSSD.time,cur_waveSSD.freq,cur_trial,'tag', 'cip');
colorbar();
xlim([1 1+length(cur_inst.time{:,trial})/1000])  
axis xy
yline(52)
yline(90)
w = waitforbuttonpress;
clf
end 


%%
f = figure;


%

