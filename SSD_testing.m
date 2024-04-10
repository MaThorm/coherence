clc
clear all
load('attout_dataset.mat')
load('attin_dataset.mat')
load("V4_dataset.mat")
bpwidth = [30 100];
toi = [2.3 4.3];
medianfiltord = 20;
save_hilbert = false;
fs = 1000;
th = 0.01;
matpath = '/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files';

[in_trials,out_trials,V4_trials] = pre_processing_pip_trials(attin_dataset,attout_dataset,V4_dataset,bpwidth,toi)
[in_trials,testing_struct] = do_SSD(in_trials,fs,th)
%[grand_struct,angles,inst_freq] = pre_processing_pip_hilb(in_trials,out_trials,V4_trials,medianfiltord,save_hilbert)

attin_inst = grand_struct.in_medfiltHilbert;
attout_inst = grand_struct.out_medfiltHilbert;
insummary = grand_struct.insummary;
outsummary = grand_struct.outsummary;
V4_inst = grand_struct.V4_medfiltHilbert;
V4summary = grand_struct.V4summary;
timebar = -1.3:0.001:5;

%% CHecking out how much the max v ariance component differs from others 
min_facts_array = nan(16,252) %252 is the max of the in_trials array, would not work with others 
len_count = 2;
c = 1; % 1 for the V1 channels, 2 for the V4 channels
for i_s = 1:length(testing_struct)
    cur_s = testing_struct(i_s);
    for i_t = 1:length(cur_s.resi_var)
        resi_max = min(cur_s.resi_var(i_t).resi_var{1, c});
        testing_struct(i_s).resi_var(i_t).facts = resi_max ./ cur_s.resi_var(i_t).resi_var{1, c};
        testing_struct(i_s).resi_var(i_t).facts = sort(testing_struct(i_s).resi_var(i_t).facts,'descend');
        if length(testing_struct(i_s).resi_var(i_t).facts) == 1
            testing_struct(i_s).resi_var(i_t).min_facts = nan;
        else
            testing_struct(i_s).resi_var(i_t).min_facts = testing_struct(i_s).resi_var(i_t).facts(2);
            min_facts_array(i_s,i_t) = testing_struct(i_s).resi_var(i_t).min_facts;
        end 
        resi_var_len(len_count) = length(cur_s.resi_var(i_t).resi_var{1, c}); % counts length of residual variance arrays
        len_count = len_count + 1;
    end  
end 
means = median(min_facts_array,2,'omitnan');

%% Plotting boxplots and violins of max variance component ratios
figure
violin(min_facts_array');
title('Highest Variance Component Divided by second highest variance component')
xlabel('Session')
ylabel('Ratio')

%%
min_facts_resh = reshape(min_facts_array.',1,[]);
violin(min_facts_resh');
title('Highest Variance Component Divided by second highest variance component')
ylabel('Ratio')
%% Histogram of number of components per trial
figure
histogram(resi_var_len')
title('Number of SSD components per trial')