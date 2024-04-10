function [in_trials,out_trials,V4_trials] = pre_processing_pip(attin_dataset,attout_dataset,V4_dataset,bpwidth,toi)
%PRE_PROCESSING_PIP Gives out a Nx1 Structure that includes the
%trialselection for each sessio (N) using Fieldtrip. Also includes a
%selection of the bandpassfilter width and time of intereset
% Input should be: 
% Attended Dataset
% Unattended Datasete
% V4 Dataset
% bpwdith: Width of the BP filter
% toi: Time period of interest from 0 to 5.3
% Creating trials and bp filtering 
% Attin trials
parfor ii = 1:length(attout_dataset)
    in_trials(ii) = do_trialselection(attin_dataset(ii).path,attin_dataset(ii).file,[attin_dataset(ii).chan attin_dataset(ii).V4chan{1}] ,attin_dataset(ii).stimno,bpwidth,toi);
end 
% Attout trials
parfor ii = 1:length(attout_dataset)
    out_trials(ii) = do_trialselection(attout_dataset(ii).path,attout_dataset(ii).file,[attout_dataset(ii).chan attin_dataset(ii).V4chan{1}] ,attout_dataset(ii).stimno,bpwidth,toi);
end 
% AttV4 trials
parfor ii = 1:length(V4_dataset)
    V4_trials(ii) = do_trialselection(V4_dataset(ii).path,V4_dataset(ii).file,V4_dataset(ii).chan,V4_dataset(ii).stimno,bpwidth,toi);
end 
end 

function [trials] = pre_processing_pip(dataset,bpwidth,toi)
%PRE_PROCESSING_PIP Gives out a Nx1 Structure that includes the
%trialselection for each sessio (N) using Fieldtrip. Also includes a
%selection of the bandpassfilter width and time of intereset
% Input should be: 
% Attended Dataset
% Unattended Datasete
% V4 Dataset
% bpwdith: Width of the BP filter
% toi: Time period of interest from 0 to 5.3
% Creating trials and bp filtering 
% Attin trials
parfor ii = 1:length(attout_dataset)
    in_trials(ii) = do_trialselection(attin_dataset(ii).path,attin_dataset(ii).file,[attin_dataset(ii).chan attin_dataset(ii).V4chan{1}] ,attin_dataset(ii).stimno,bpwidth,toi);
end 
% Attout trials
parfor ii = 1:length(attout_dataset)
    out_trials(ii) = do_trialselection(attout_dataset(ii).path,attout_dataset(ii).file,[attout_dataset(ii).chan attin_dataset(ii).V4chan{1}] ,attout_dataset(ii).stimno,bpwidth,toi);
end 
% AttV4 trials
parfor ii = 1:length(V4_dataset)
    V4_trials(ii) = do_trialselection(V4_dataset(ii).path,V4_dataset(ii).file,V4_dataset(ii).chan,V4_dataset(ii).stimno,bpwidth,toi);
end 
end 