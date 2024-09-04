function new_ft = do_peak_cut(ft_cut,ft_env,peak_time,peak_prcntl,min_time)

% Function that takes as any ft_struct as an input and cuts it according to
% the peaks found in an ft_envelope structure of the same shape (as calculated by hilbert). 
% Inputs should be: 
% ft_cut: the ft structure to be cut
% ft_env: the corresponding ft_envelope structure
% peak_time: The time used surrounding the peaks
% Peak prctl: The percentile of the envelope peaks used
% min_time: Minimum time of snippets, otherwise they are discarded

for i_s = 1:length(ft_env)
    cut_off = create_peak_cut(ft_env(i_s),peak_prcntl);
    for i_t = 1:length(ft_env(i_s).trial)
        t_reject{i_s,i_t} = find_peak_widths(ft_env(i_s).trial{1,i_t},cut_off,peak_time);
    end 
end 
%
peak_cell = gen_reject(ft_env,t_reject);

%
cfg = [];
cfg.artfctdef.reject = 'partial';
cfg.artfctdef.minaccepttim = min_time;
for ii = 1:length(ft_cut)
    cfg.artfctdef.xxx.artifact = peak_cell{1,ii};
    new_ft(ii) = ft_rejectartifact(cfg,ft_cut(ii));
end 


function r_new_widths = find_peak_widths(cur_data,thresh,width)
% function that finds peaks of an 1byX array, then compares it to some threshold.
% Around all remaining peaks the width is then added to receive a peak area
% If two peaks areas collide, they are added together. 
[pks,locs] = findpeaks(cur_data);
pks = pks(cur_data(locs) > thresh);
locs = locs(cur_data(locs) > thresh);
peak_widths = zeros(2,length(pks));
peak_widths(1,:) = locs - width;
peak_widths(2,:) = locs + width;
peak_widths(peak_widths <= 0) = 0;
peak_widths(peak_widths >= length(cur_data)) = length(cur_data) - 1;

% Removing overlapping peaks 
new_widths = peak_overlap(peak_widths);

% Reversing the top array because I'm stupid
r_new_widths = [];
% Consider time before first peak
if ~isempty(new_widths)
    if new_widths(1,1) > 1
        r_new_widths = [r_new_widths [1; new_widths(1,1)-1]];
    end 
    
    % Between first and last peak
    for ii = 1:size(new_widths,2) - 1
        low = new_widths(2,ii);
        high = new_widths(1,ii + 1);
        r_new_widths = [r_new_widths [(low +1) ;(high - 1)]];
    end 
    
    % Consider time after last peak
    if new_widths(2,end) < length(cur_data)
        r_new_widths = [r_new_widths [(new_widths(2,end) + 1); length(cur_data)]];   
    end 
elseif isempty(new_widths)
    r_new_widths = [1;length(cur_data)];
end 
end 

function new_widths = peak_overlap(peak_widths)
% Checking if peaks overlap
if ~isempty(peak_widths)
    new_widths = peak_widths(:, 1);
    for i = 2:size(peak_widths, 2)
        last_merged_event = new_widths(:, end);
        current_event = peak_widths(:, i);   
        if current_event(1) <= last_merged_event(2)
            new_widths(2, end) = max(last_merged_event(2), current_event(2));
        else
            new_widths = [new_widths, current_event];
        end
    end
else
    new_widths = [];
end
end 


function [peak_cell] = gen_reject(ft_struct,peak_widths)
% Generates a 2 by X rejection array that can be used for ft_rejectartifact for each session
% Input should be: an ft_struct, a session x trial cell array with the starting and end points of each non-rejected times of a trial
for i_s = 1:size(peak_widths,1)
    sess_peaks = [];
    for i_t = 1:length(ft_struct(i_s).trial)
        t_start = ft_struct(i_s).sampleinfo(i_t,1);
        for i_p = 1:size(peak_widths{i_s,i_t},2)
            p_start = t_start + peak_widths{i_s,i_t}(1,i_p);
            p_end = t_start + peak_widths{i_s,i_t}(2,i_p);
            sess_peaks = [sess_peaks' [p_start;p_end]]';
        end 
    end 
    peak_cell{i_s} = sess_peaks;
end 
end 


function cut_off = create_peak_cut(ft_struct,peak_prctile)
% Generates a cut-off amplitude for a single session of an ft_struct using the envelope information of the hilbert transform
% Input should be: single session ft_struct, percentile of accepted peaks
peaks = [];
for ii = 1:length(ft_struct.trial)
    peaks = [peaks findpeaks(ft_struct.trial{1,ii})];
end 
cut_off = prctile(peaks,peak_prctile)
end 

end 