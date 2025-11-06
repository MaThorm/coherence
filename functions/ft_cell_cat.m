function [concarray] = ft_cell_cat(ft_struct)
%FT_CELL_CAT Concatenates the .trial field of a ft struct into one big cell
%array
my_cell = struct2cell(ft_struct)
new_cell(:,1) = my_cell(5,:,:);
% new_cell(:,2) = my_cell(6,:,:);

concarray = new_cell{1,:};

for ii = 2:length(new_cell)
    concarray = [concarray new_cell{ii,:}];
%     concarray{2,:} = [concarray new_cell{ii,2}]
end 
end

