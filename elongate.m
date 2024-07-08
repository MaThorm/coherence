function [new_struct] = elongate(ft_struct)
%ELONGATE Elongates an FT_struct so that one V1 electrode is always paired
%with a single V4 electrode. Does this by adding both trials and labels
%accordingly. 
counter = 1;

% Probably the dumbest way I could solve the elongate: first reordering
% according to recording days. Followed by which V1 electrodes had multiple
% V4 electrodes. 
new_order = [1 7 2 8 9 3 10 11 12 13 14 15 4 5 16 6];
rep_mat = [1 1 2 2 1 2 1 1 1 1 1 1 1 3 3 1];
% Reordering the structs so that they have the same order as V4
for ii = 1:length(ft_struct)
    re_ft_struct(ii) = ft_struct(new_order(ii))
end 

for i_s = 1:length(ft_struct)
    cur_sess = re_ft_struct(i_s)
    for i_e = 1:rep_mat(i_s)
        new_struct(counter) = cur_sess
        if isfield(cur_sess,'trial') == 1
            for i_t = 1:length(cur_sess.trial)
                tcell{:,i_t}(1,:) = cur_sess.trial{:,i_t}(1,:);
                tcell{:,i_t}(2,:) = cur_sess.trial{:,i_t}(i_e+1,:);
            end 
            new_struct(counter).trial = tcell;
        end 
        if length(cur_sess.label) > 1
            lcell{:,1} = cur_sess.label{:,1};
            lcell{:,2} = cur_sess.label{:,i_e+1};
            new_struct(counter).label = lcell;
        end 
        clear tcell lcell
        counter = counter + 1;
end 
end

