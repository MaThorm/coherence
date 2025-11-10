% Script to test the envelope information of the hilbert transform
clc
clear 
% Trial selection 
% Loading
load("/data/projects/V1V4coherence/02_analysis_max/git_repos/params.mat")
load('attin_dataset.mat')
in_trials = pre_processing_pip_trials(attin_dataset,params.bpfilt,params.bpwidth,params.order);
for ii = 1:length(in_trials)
     cfg = [];
     cfg.channel = in_trials(ii).label{1,1};
     in_trials(ii) = ft_selectdata(cfg,in_trials(ii));
%      in_trials(ii).label = 'V1'
end 
% Hilbert transform and inst. freq
[hilb_angles.wrapped.in,hilb_env.in]  = pre_processing_pip_hilb(in_trials);
[hilb_angles.in,inst_freq.in,filt_data.in]  = pre_processing_pip_filtinst(hilb_angles.wrapped.in,params.filttype,params.framelen,params.filtord,params.toi);

in_trials = do_toi_cut(in_trials,params.toi)
hilb_env.in = do_toi_cut(hilb_env.in,params.toi)

%%
peak_time = params.peak_time;
peak_prcntl = params.peak_prcntl;
min_time = params.min_time;
new_inst_freq = do_peak_cut(filt_data.in,hilb_env.in,peak_time,peak_prcntl,min_time);
[inst_avg.in inst_avg.g_in] = do_grand_avg(new_inst_freq);


%% Adding my own artifact rejection 
peak_time = 30;
for i_s = 1:length(hilb_env.in)
    cut_off = create_peak_cut(hilb_env.in(i_s),80);
    for i_t = 1:length(hilb_env.in(i_s).trial)
        t_reject{i_s,i_t} = find_peak_widths(hilb_env.in(i_s).trial{1,i_t},cut_off,peak_time);
    end 
end 
%
peak_cell = gen_reject(hilb_env.in,t_reject);

%
cfg = [];
cfg.artfctdef.reject = 'partial';
cfg.artfctdef.minaccepttim = 0.1;
for ii = 1:length(in_trials)
    cfg.artfctdef.xxx.artifact = peak_cell{1,ii}
    new_in_trials(ii) = ft_rejectartifact(cfg,hilb_env.in(ii));
end 

    
%%
hold on
plot(in_trials(1).time{1},in_trials(1).trial{1,1});
plot(hilb_env.in(1).time{1},hilb_env.in(1).trial{1,1});
for ii = 1:length(pks)
    xline(hilb_env.in(1).time{1}(locs(ii)));
end 
hold off
%% Testing automatic artifact rejection
cfg = [];
cfg.continuous = 'no';
cfg.artfctdef.threshold.channel   = in_trials(1).label;
cfg.artfctdef.threshold.bpfilter  = 'no'
% cfg.artfctdef.threshold.bpfreq    = [0.3 30]
% cfg.artfctdef.threshold.bpfiltord = 4
cfg.artfctdef.threshold.min       = -5;
[cfg, artifact] = ft_artifact_threshold(cfg, in_trials(1))
% artifact = artifact(:,1:2)
cfg = [];
cfg.artfctdef.reject = 'partial';
cfg.artfctdef.xxx.artifact = artifact;
trials = ft_rejectartifact(cfg, in_trials(1));
%%
%% Functions
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
