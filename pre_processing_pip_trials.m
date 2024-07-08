function [trials] = pre_processing_pip_trials(dataset,bpfilt,bpwidth)
%PRE_PROCESSING_PIP Gives out a Nx1 Structure that includes the
%trialselection for each sessio (N) using Fieldtrip. Also includes a
%selection of the bandpassfilter width and time of intereset
% Input should be: 
% A dataset
% bpwdith: Width of the BP filter
% toi: Time period of interest from 0 to 5.3
% Creating trials and bp filtering 
parfor ii = 1:length(dataset)
    if ~isfield(dataset,'V4chan') %V1 trials have the corresponding V4channels included for PLV analysis. V4 trials only have one channel
        chan = dataset(ii).chan; 
    else 
        chan = [dataset(ii).chan dataset(ii).V4chan{1}];
    end 
    trials(ii) = do_trialselection(dataset(ii).path,dataset(ii).file, chan, dataset(ii).stimno, bpfilt, bpwidth);
end 
