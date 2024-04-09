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
            resi_var = [];
            for comp = 1:size(cur_SSD,1)
                resi = cur_trial - cur_SSD(comp,:);
                testing_struct(i_sess).resi(ii).resi{i_c}(comp,:)  = resi;
                resi_var(comp) = var(resi) / var(cur_trial); %Die formel hier für die Residual variance ist nicht unbedingt der korrekte weg, im Lowet Paper die NMSE
            end 
            testing_struct(i_sess).resi_var(ii).resi_var{i_c}(:) = resi_var;
            [temp,locmax_var] = min(resi_var);
            new_trials(i_sess).trial{1, ii}(i_c,:) = cur_SSD(locmax_var,:);
            testing_struct(i_sess).max_SSD(ii).max_SSD{i_c}(:) = cur_SSD(locmax_var,:);
        end 
    end 
end  
end