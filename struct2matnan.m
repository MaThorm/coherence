function [nan_array] = struct2matnan(trials)
%Function that takes an n-dimensional fieldtrip struct. And gives back a n
%by x by t matrix. N = number of sessions, x = max number of trials, t =
%max triallength. Matrix is filled with nans where there is no value
for ii = 1:length(trials)
    ncell{ii} = cell2matnan(trials(ii).trial);
end
maxnumrow = cellfun(@(x) size(x,1),ncell,'UniformOutput',false);
maxnumcol = cellfun(@(x) size(x,2),ncell,'UniformOutput',false);
size_row = max(cell2mat(maxnumrow));
size_col = max(cell2mat(maxnumcol));
nan_array = nan(length(trials),size_row,size_col);

for i_s = 1:length(trials)
    nan_array(i_s,1:maxnumrow{i_s},:) = squeeze(cell2matnan(trials(i_s).trial));
end 
