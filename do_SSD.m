function [new_trials,testing_struct] = do_SSD(trials,fs,th)
%DO_SSD Function that takes a trialselected fieldtrip file and, using that,
%switches out each trial with the biggest component of the SSD function.
% trials: fieldtrip result of do_trialselection
% Fs: sampling frequency
% th: threshold which defines the variance of the residual (default 0.01),
new_trials = trials;
parfor i_sess = 1:length(trials)
    for ii = 1:length(trials(i_sess).trial)
        for i_c = 1:size(trials(i_sess).trial{1, ii},1)
            cur_trial = trials(i_sess).trial{1, ii}(i_c,:);
            cur_SSD = SSD(cur_trial,fs,th);
            testing_struct(i_sess).SSD(ii).SSD{i_c}(:,:) = cur_SSD; %Storing different channels in different cells since number of components might differ
            nmse = [];
            for comp = 1:size(cur_SSD,1)
                mse = mean((cur_trial - cur_SSD(comp,:)).^2);
                testing_struct(i_sess).resi(ii).resi{i_c}(comp,:)  = mse;
                nmse(comp) = mse / var(cur_trial); 
            end 
            testing_struct(i_sess).resi_var(ii).resi_var{i_c}(:) = nmse;
            [temp,locmax_var] = min(nmse);
            new_trials(i_sess).trial{1, ii}(i_c,:) = cur_SSD(locmax_var,:);
            testing_struct(i_sess).max_SSD(ii).max_SSD{i_c}(:) = cur_SSD(locmax_var,:);
        end 
    end 
end  
end
