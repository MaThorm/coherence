function [sess_mean,grand_mean] = do_grand_avg(cur,varargin)
%DO_GRAND_AVG Summary of this function goes here
%  Varargin{1} should be true if the struct is elonagated
% Varargin{2} should be the parameter that is averaged
nInputs = numel(varargin)
% Changing labels 
cfg = [];
cfg.channel = 'V1';
if nInputs > 0 %Checking if varargin is given, if yes V4 mean is also calculated. 
    if varargin{1} == true
        cfg.channel = {'V1','V4'};
    end 
end 

for ii = 1:length(cur)
    cur(ii).label{1,1} = 'V1';
    for i_c = 1:length(cur(ii).label)-1
        if nInputs > 0
            if varargin{1} == true
                cur(ii).label{1,2} = 'V4';
            end 
        end 
    end 
    sess_mean(ii) = ft_timelockanalysis(cfg,cur(ii));
end 

if length(cur) == 16 
    grand_mean = ft_timelockgrandaverage(cfg,sess_mean(1),sess_mean(2),sess_mean(3) ...
        ,sess_mean(4),sess_mean(5),sess_mean(6),sess_mean(7),sess_mean(8) ...
        ,sess_mean(9),sess_mean(10),sess_mean(11),sess_mean(12),sess_mean(13) ...
        ,sess_mean(14),sess_mean(15),sess_mean(16));
elseif length(cur) == 23
    grand_mean = ft_timelockgrandaverage(cfg,sess_mean(1),sess_mean(2),sess_mean(3) ...
    ,sess_mean(4),sess_mean(5),sess_mean(6),sess_mean(7),sess_mean(8) ...
    ,sess_mean(9),sess_mean(10),sess_mean(11),sess_mean(12),sess_mean(13) ...
    ,sess_mean(14),sess_mean(15),sess_mean(16),sess_mean(17),sess_mean(18), ...
    sess_mean(19),sess_mean(20),sess_mean(21),sess_mean(22),sess_mean(23));
else
    error('That length is weird dude')
end 
end

