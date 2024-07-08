clear 

% Hilbert Angles, instantaneous frequency, filtered instantaneous
% frequency, instantaneous change
load("/data/projects/V1V4coherence/02_analysis_max/git_repos/params.mat")

% loading inst freq struct 
filename = sprintf("inst_freq_toi%.1f-%.1f_lower%i_upper%i_filtord%i.mat",params.toi(1),params.toi(2),params.lower,params.upper,params.filtord);
m = matfile(fullfile(params.matpath,'Inst_freq',filename),'Writable',true);
filt_data = m.inst;

%%
cut_perc = 0.05;
[cut_data.in,cut.in] = do_freqvar_cut(filt_data.in,cut_perc);
[cut_data.out,cut.out] = do_freqvar_cut(filt_data.out,cut_perc);
[cut_data.V4,cut.V4] = do_freqvar_cut(filt_data.V4,cut_perc);
%%
[filt_means.in.s_avg,filt_means.in.g_avg] = do_grand_avg(filt_data.in)
[cut_means.in.s_avg,cut_means.in.g_avg] = do_grand_avg(cut_data.in)
[filt_means.out.s_avg,filt_means.out.g_avg] = do_grand_avg(filt_data.out)
[cut_means.out.s_avg,cut_means.out.g_avg] = do_grand_avg(cut_data.out)

%%
f = figure; 
f.Units = 'normalized';
f.Position = [0 0 1 1]
subplot(2,1,1)
plot(filt_means.in.g_avg.avg)
hold on 
plot(cut_means.in.g_avg.avg)
hold off
legend({'Not Cut','Cut'})
subplot(2,1,2)
plot(filt_means.out.g_avg.avg)
hold on 
plot(cut_means.out.g_avg.avg)
hold off
legend({'Not Cut','Cut'})

%%
big = filt_data.in(11).trial{1,17}(1,:);
cur_trial = big;
for ii = 1:100
    cur_trial = filt_data.in(11).trial{1,ii}(1,:);
    cur_perc = prctile(cur_trial,[25,75]);
    cur_IQR = iqr(cur_perc)
    l_l = cur_perc(1) - 1.5 * cur_IQR;
    l_u = cur_perc(2) + 1.5 * cur_IQR;
    cut_num =  length(cur_trial(cur_trial > l_u | cur_trial < l_l));
    cur_limit = 0.1 * length(cur_trial);
    plot(cur_trial)
    hold on 
    yline(cur_perc(1))
    yline(cur_perc(2))
    yline(l_l)
    yline(l_u)
    w = waitforbuttonpress;
    clf 
end 
%% 
function [cut_data,cut] = do_freqvar_cut(cur_data,outlier_perc)
for i_s = 1:length(cur_data)
    for i_t = 1:length(cur_data(i_s).trial)
        cur_trial = cur_data(i_s).trial{1,i_t}(1,:);
        cur_perc = prctile(cur_trial,[25,75]);
        cur_IQR = iqr(cur_trial);
        l_l = cur_perc(1) - 1.5 * cur_IQR;
        l_u = cur_perc(2) + 1.5 * cur_IQR;
        cut(i_s).num(i_t) = length(cur_trial(cur_trial > l_u | cur_trial < l_l));
        cur_limit = outlier_perc * length(cur_trial);
        if cut(i_s).num(i_t) > cur_limit
            cut(i_s).logi(i_t) = false; 
        else
            cut(i_s).logi(i_t) = true;
        end
    end 
    cfg = [];
    cfg.trials = cut(i_s).logi;
    cut_data(i_s) = ft_selectdata(cfg,cur_data(i_s))
end 
end 