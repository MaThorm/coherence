function [wmean_str, g_wmean_str] = do_weight_mean(ft_struct,ft_env,chan_sel)
% Alternative to do_grand_avg. Instead of taking the regular mean, this one is weighted by the
% squared envelope. 
% Inputs should be: 
% ft_struct: the ft_struct to be averaged (instantaneous frequency)
% ft_env: the corresponding envelope ft_struct
% chan_sel: if the ft_struct is elongated this should be false, otherwise true
% Output is: 
% wmean_str: the weighted average per session 
% g_wmean_str: the weighted grand average
% Right now this cannot compute the corresponding V4 channel in the V1 data

% Turning both ft_structs into 3D arrays
ft_str_array = do_trial_array(ft_struct);
env_str_array = do_trial_array(ft_env);

% Optionally taking envelope to the power of 2
env_str_array = env_str_array.^2;

% Calculating the weighted mean 
ft_str_mean = weighted_mean(ft_str_array,env_str_array,2);
ft_str_gmean = mean(ft_str_mean,1);

[wmean_str, g_wmean_str]  = do_grand_avg(ft_struct,true);
for ii = 1:length(wmean_str)
    wmean_str(ii).avg = squeeze(ft_str_mean(ii,:,:));
end 
g_wmean_str.avg = squeeze(ft_str_gmean);


% Calculations dependingon whether data is elongated or not
cfg = [];
if chan_sel == true
    for ii = 1:length(wmean_str)
        cfg.channel = wmean_str(ii).label{1,1};
        wmean_str(ii) = ft_selectdata(cfg,wmean_str(ii));
        wmean_str(ii).avg = squeeze(wmean_str(ii).avg(1,:));
    end 
cfg = [];
cfg.channel = g_wmean_str.label{1,1};
g_wmean_str = ft_selectdata(cfg,g_wmean_str);
g_wmean_str.avg = squeeze(g_wmean_str.avg(1,:));
elseif chan_sel == false;
    g_wmean_str.avg = squeeze(g_wmean_str.avg(1:2,:));
    for ii = 1:length(wmean_str)
        wmean_str(ii).avg = squeeze(wmean_str(ii).avg(1:2,:));
    end 
end 

function trial_array = do_trial_array(trial_struct)
max_trials = 0;
max_trial_len = 0;
% Finding out max length of trials and trial length
for i_s = 1:length(trial_struct)
    if length(trial_struct(i_s).trial) > max_trials
        max_trials = length(trial_struct(i_s).trial);
    end 
    for i_t = 1:length(trial_struct(i_s).trial)
        if length(trial_struct(i_s).trial{1,i_t}) > max_trial_len
            max_trial_len = length(trial_struct(i_s).trial{1,i_t});
        end 
    end 
end 

% Creating empty nan array with those data
% Importantly this leaves the V4 channel out, so that this cannot be used for V1/V4 calcs
trial_array = nan(length(trial_struct),max_trials,4,max_trial_len);
for i_s = 1:length(trial_struct)
    for i_t = 1:length(trial_struct(i_s).trial)
        trial_array(i_s,i_t,1:size(trial_struct(i_s).trial{1,i_t},1),1:length(trial_struct(i_s).trial{1,i_t})) = trial_struct(i_s).trial{1,i_t};
    end 
end 
end 

function mean_array = weighted_mean(mean_array,weight_array,dim)
% Step 1: Weighted Mean Over Dimension K
% Compute the numerator and denominator for each (i, j)
numerator_k = sum(mean_array .* weight_array, dim,'omitnan');  
denominator_k = sum(weight_array, dim,'omitnan');

% Avoid division by zero by adding a small epsilon where denominator is zero
epsilon = 1e-12;
denominator_k(denominator_k == 0) = epsilon;

% Compute the weighted mean over k
mean_array = squeeze(numerator_k ./ denominator_k);     % Resulting size is (I x J)
end 


end

