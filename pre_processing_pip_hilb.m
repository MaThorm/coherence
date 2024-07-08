function [hilbert_angles_wr] = pre_processing_pip_hilb(trials)
%PRE_PROCESSING_PIP_HILB Perofrms the hilbert transform on fieldtrip
%structs. Gives back the phase angles
% Input should be
% trials: a trialselected fieldtrip struct

parfor ii = 1:length(trials)
    [hilbert_angles_wr(ii)] = do_hilbert(trials(ii)); 
end 



end