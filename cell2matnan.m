function [aMat] = cell2matnan(cell1,row)
%cell2matnan Transforms a cell with 2-by-X matrices as contents into a matrix filled
%with nans, using only the first row of each matrix.

    % Extract the first row from each 2-by-X matrix in the cell array
    firstRows = cellfun(@(c) c(row,:), cell1, 'UniformOutput', false);
    
    % Find the maximum number of columns in any of the first rows
    maxNumCol = max(cellfun(@(c) size(c,2), firstRows));
    
    % Pad each first row with NaNs to have the same number of columns
    paddedRows = cellfun(@(c) [c, nan(1, maxNumCol - size(c,2))], firstRows, 'UniformOutput', false);
    
    % Concatenate the padded first rows vertically
    aMat = vertcat(paddedRows{:});
end 