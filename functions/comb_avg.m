function [meanmean_array,sdsd_array] = comb_avg(input_str)
%COMB_AVG Function that combines and averages  data types
for ii = 1:length(input_str)
    a = input_str(ii).input
    mean_cell = struct2cell(a);
    mean_cell = squeeze(mean_cell(9,:,:));
    mean_array(ii,:,:) = cell2matnan(mean_cell');
    sd_cell = struct2cell(a);
    sd_cell = squeeze(sd_cell(10,:,:));
    sd_array(ii,:,:) = cell2matnan(sd_cell');
end
sdsd_array = squeeze(mean(sd_array,1));
meanmean_array = squeeze(mean(mean_array,1));
end 
