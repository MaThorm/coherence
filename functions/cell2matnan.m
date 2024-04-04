function [aMat] = cell2matnan(cell1)
%CELL_NANMEAN# Summary of this function goes here
%   Detailed explanation goes here
    maxNumCol = max(cellfun(@(c) size(c,2), cell1));
    aMat = cell2mat(cellfun(@(c){[c nan(1,maxNumCol-numel(c))]}, cell1)');
end

