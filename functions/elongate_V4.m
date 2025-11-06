function [cur_V4] = elongate_V4(V4_struct)
%ELONGATES and Orders the V4 Struct
v4_long = [1 1 2 3 2 3 4 5 6 7 8 9 9 10 11 12 13 14 15 13 14 15 16]
for ii = 1:length(v4_long)
    cur_V4(ii) = V4_struct(v4_long(ii))
end 
end

