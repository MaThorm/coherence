function [hilbert_angles_wr,hilb_env] = pre_processing_pip_hilb(trials)
%PRE_PROCESSING_PIP_HILB Perofrms the hilbert transform on fieldtrip
%structs. Gives back the phase angles and envelope
% Input should be
% trials: a trialselected fieldtrip struct

parfor ii = 1:length(trials)
    [hilbert_angles_wr(ii),hilb_env(ii)] = do_hilbert(trials(ii)); 
end 



end