clc
clear all
matpath = '/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files';
toi = [2.3 4.3]; % Time region of interest [2.3 4.3] is MC 2&3
fs = 1000; % Sampling frequency of the signal
th = 0.01; % residual variance threshhold   
lower = 52;
upper = 90;
in_ssd_trials = load(fullfile(matpath,'SSD',sprintf('in_ssd_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))));
out_ssd_trials = load(fullfile(matpath,'SSD',sprintf('out_ssd_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))));
V4_SSD_trials = load(fullfile(matpath,'SSD',sprintf('V4_ssd_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))));
load(fullfile(matpath,'SSD',sprintf('inc_in_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))))
load(fullfile(matpath,'SSD',sprintf('inc_out_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))))
load(fullfile(matpath,'SSD',sprintf('inc_V4_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))))
%%

sess = 8

m = matfile(fullfile(matpath,'SSD_testing',sprintf('testing_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))),'Writable',true);
temp = m.in;
testing_in = temp(sess);
temp = m.out;
testing_out = temp(sess);
temp = m.V4
testing_V4 = temp(sess);
clear temp

%%

n_comp = {};
for ii = 1:1:length(testing_out)
    n_comp{ii} = [cellfun(@(x) length(find(x >= 0.49999)),testing_out(ii).q(:,1)); cellfun(@(x) length(find(x >= 0.49999)),testing_in(ii).q(:,1));cellfun(@(x) length(find(x >= 0.49999)),testing_V4(ii).q(:,1))]
end 
histogram(cell2mat(n_comp))
title(sprintf("Number of components with a fraction > 0.5, session %d",sess))
%% Further calculations
in_cur = testing_in;
out_cur = testing_out;
max_in_SSD = nan(length(in_cur.p),length(freq));
max_out_SSD = nan(length(out_cur.p),length(freq));

freq = in_cur.f{1,1};
for ii = 1:length(in_cur.SSD)
    max_one_in(ii) = find(in_cur.bounded{ii,1} == max(in_cur.bounded{ii,1})); % location of max gamme component 
    max_in_SSD(ii,:) = in_cur.p{ii,1}(:,max_one_in(ii)); %In components with the max gamma power per trial    
    peak_in(ii) = find(max_in_SSD(ii,:) == max(max_in_SSD(ii,:))); %Location of peaks per SSD component 
    q_max_in(ii) = in_cur.q{ii,1}(:,max_one_in(ii)); % Fraction of biggest gamma component 
end 

for ii = 1:length(out_cur.SSD)
    max_one_out(ii) = find(out_cur.bounded{ii,1} == max(out_cur.bounded{ii,1}));
    max_out_SSD(ii,:) = out_cur.p{ii,1}(:,max_one_out(ii));
    peak_out(ii) = find(max_out_SSD(ii,:) == max(max_out_SSD(ii,:)));  
    q_max_out(ii) = out_cur.q{ii,1}(:,max_one_out(ii)); 
end 


%% Plottig all sessions fooof vs SSD components

subplot(1,2,1)
plot(log(foof_summary.in_all.osc_alt(1).freq), log(foof_summary.in_all.osc_alt(sess).powspctrm));
hold on 
plot(log(foof_summary.in_all.osc_alt(1).freq), log(foof_summary.out_all.osc_alt(sess).powspctrm));
plot(log(foof_summary.in_all.osc_alt(1).freq), log(foof_summary.V4_all.osc_alt(sess).powspctrm));
x = str2double(xticklabels);
xticklabels(exp(x))
xlabel('log-freq'); ylabel('log-power'); grid on;
legend({'in','out','V4'},'location','southwest');
title(sprintf('Fooof sess %d',num));
hold off
subplot(1,2,2)
histogram(n_comp)

%% Plotting only components above >0.5. Green one has highest fraction 
cur_test = testing_out;
freq = cur_test.f{1,1}
for ii = 1:length(testing_in.SSD)
    cur_pos = find(cur_test.q{ii,1} > 0.4999)
    for i_comp = 1:length(cur_pos)
        plot(freq,cur_test.p{ii,1}(:,cur_pos(i_comp)))
        xline([52])
        xline([110])        
        hold on
        title(sprintf('Trial %d,  % residual variance threshhold   
lower = 5component %d',ii,i_comp))
    end 
    max_one = find(cur_test.q{ii,1} == max(cur_test.q{ii,1}));
    plot(freq,cur_test.p{ii,1}(:,max_one),'g')
    hold off 
    w = waitforbuttonpress
    clf
end 

%% Plotting all, green one has highest fraction
cur_test = testing_in;
freq = cur_test.f{1,1}
for ii = 1:length(testing_in.SSD)
    for i_comp = 1:size(cur_test.p{ii,1},2)
        plot(freq,cur_test.p{ii,1}(:,i_comp))
        xline([52])
        xline([110])
        hold on
    end 
    max_one = find(cur_test.q{ii,1} == max(cur_test.q{ii,1}));
    plot(freq,cur_test.p{ii,1}(:,max_one),'g')
    hold off 
    xlim([1 150])
    w = waitforbuttonpress
    clf
end 


%% Plotting all, green one is biggest
cur_test = testing_in;
freq = cur_test.f{1,1}
saving = true;
for ii = 1:length(cur_test.SSD)
    for i_comp = 1:size(cur_test.p{ii,1},2)
        plot(freq,cur_test.p{ii,1}(:,i_comp))
        xlim([1 150])

        hold on
    end 
    sorted = sort(cur_test.bounded{ii,1});
    max_one = find(cur_test.bounded{ii,1} == sorted(end));
    second_one = find(cur_test.bounded{ii,1} == sorted(end-1));
    %max_one = find(cur_test.bounded{ii,1} == max(cur_test.bounded{ii,1}));
    plot(freq,cur_test.p{ii,1}(:,max_one),'g','LineWidth',2)
    plot(freq,cur_test.p{ii,1}(:,second_one),'r','LineWidth',2)
    xline([52],'--')
    xline([90],'--')
    if inc_in(sess).inc{ii,1} == true
        truth = 'accepted';
    else 
        truth = 'rejected';
    end 
    title(sprintf('Session: %d, trial %d, %s ',sess,ii,truth))
    ylim([0 2*max(cur_test.p{ii,1}(:,max_one))])
    hold off 
    if saving == true
        foldername = sprintf("/data/projects/V1V4coherence/03_results_max/SSD/SSD_powerplots/session%d/%dHzto%dHz",sess,lower,upper)
        if ~exist(foldername,'dir')
            mkdir(foldername)
        end 
        saveas(gcf,fullfile(foldername, sprintf("trial%d.jpg",ii)))
    else 
        w = waitforbuttonpress;
    end 
    clf
    %
end
%% Plotting some example trials 
cur_test = testing_in;
trial5 = [67 101 105 118 106 141 159];
trial8 = [43 47 92 136 152 156 162 175 179 185 191];
cur_trial = trial5;
freq = cur_test.f{1,1}
for ii = 1:length(cur_trial)
    subplot(3,1,1)
    sgtitle(sprintf('Session %d trial %d',sess,cur_trial(ii)))
    sorted = sort(cur_test.bounded{cur_trial(ii),1});
    max_one = find(cur_test.bounded{cur_trial(ii),1} == sorted(end));
    plot(freq,cur_test.p{cur_trial(ii),1}(:,max_one),'g','LineWidth',2);
    xlim([1 150])
    xline([52],'--')
    xline([90],'--')
    subplot(3,1,2)
    plot(cur_test.SSD{cur_trial(ii),1}(max_one,:))
    title('SSD Component')
    subplot(3,1,3)
    plot(in_trials(sess).trial{1,cur_trial(ii)}(1,:))
    title('Original Data')
    w = waitforbuttonpress;
    clf 
end 


%%
sig = in_trials(1).trial{1, 1};
time = in_trials(1).time{1}
[spectrum, ntaper, freqoi] = ft_specest_mtmfft(sig,time,'taper','hanning');

%%
session = 1;
trial = 1;
comp_num = 1;
freq = in(session).f
for trial = 1:length(in(session).SSD)
    for comp_num = 1:size(in(session).SSD{trial,1},1)
        sig = in(session).SSD{trial, 1};
        time = 0.001:0.001:length(sig)/1000;
        [spectrum, ntaper, freqoi] = ft_specest_mtmfft(sig,time,'taper','dpss','tapsmofrq',);
        subplot(2,1,1)
        spectrum = mean(spectrum,1)
        plot(freqoi,squeeze(abs(spectrum(:,comp_num,:))));
        xlim([1 200]) 
        [spectrum, ntaper, freqoi] = ft_specest_mtmfft(sig,time,'taper','hanning','tapsmofrq',1);
        subplot(2,1,2)
        plot(freqoi,squeeze(abs(spectrum(:,comp_num,:))));
        xlim([1 200])         
        w = waitforbuttonpress
        clf 
    end 
end 
