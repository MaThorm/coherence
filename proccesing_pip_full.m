clc
clear 

%% Trial selection 
% Loading
load("/data/projects/V1V4coherence/02_analysis_max/git_repos/params.mat")
load('attout_dataset.mat')
load('attin_dataset.mat')
load('V4_dataset.mat')

% Trialselection 
[in_trials] = pre_processing_pip_trials(attin_dataset,params.bpfilt,params.bpwidth)
[out_trials] = pre_processing_pip_trials(attout_dataset,params.bpfilt,params.bpwidth)
[V4_trials] = pre_processing_pip_trials(V4_dataset,params.bpfilt,params.bpwidth)

%Saving
LogicalStr = {'false', 'true'};
m = matfile(fullfile(params.matpath,'trials',sprintf('trials_bpfilt%s.mat',LogicalStr{params.bpfilt+1})),'Writable',true);
m.in = in_trials;
m.out = out_trials;
m.V4 = V4_trials;

%% SSD calculation
if params.bptype == "SSD"
    clear
    load("/data/projects/V1V4coherence/02_analysis_max/git_repos/params.mat")
    
    % Loading
    LogicalStr = {'false', 'true'};
    m = matfile(fullfile(params.matpath,'trials',sprintf('trials_bpfilt%s.mat',LogicalStr{params.bpfilt+1})),'Writable',true);
    in_trials = m.in;
    out_trials = m.out;
    V4_trials = m.V4;
    
    % Taking time of interesdt
    in_trials = do_toi_cut(in_trials,params.toi);
    out_trials = do_toi_cut(out_trials,params.toi);
    V4_trials = do_toi_cut(V4_trials,params.toi);

    % Performing SSD
    [in_SSD_trials,inc_in,testing_in] = do_SSD(in_trials,params.fs,params.th,params.lower,params.upper,params.toi);
    [out_SSD_trials,inc_out,testing_out] = do_SSD(out_trials,params.fs,params.th,params.lower,params.upper,params.toi);
    [V4_SSD_trials,inc_V4,testing_V4] = do_SSD(V4_trials,params.fs,params.th,params.lower,params.upper,params.toi);
    
    % Saving
    m = matfile(fullfile(params.matpath,'SSD',sprintf('ssd_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',params.th,params.lower,params.upper,10,params.toi(1),params.toi(2))),'Writable',true);
    m.in = in_SSD_trials;
    m.out = out_SSD_trials;
    m.V4 = V4_SSD_trials;
    m = matfile(fullfile(params.matpath,'SSD_testing',sprintf('testing_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',params.th,params.lower,params.upper,10,params.toi(1),params.toi(2))),'Writable',true);
    m.in = testing_in;
    m.out = testing_out;
    m.V4 = testing_V4;
end 

%% Hilbert Angles, instantaneous frequency, filtered instantaneous
% frequency, instantaneous change
clear
load("/data/projects/V1V4coherence/02_analysis_max/git_repos/params.mat")

% Loading either trial or SSD struct
if params.bptype == "SSD"
    m = matfile(fullfile(params.matpath,'SSD',sprintf('ssd_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',params.th,params.lower,params.upper,10,params.toi(1),params.toi(2))),'Writable',true);
    in_trials = m.in;
    out_trials = m.out;
    V4_trials = m.V4;
elseif params.bptype == "bpfilt"
    LogicalStr = {'false', 'true'};
    m = matfile(fullfile(params.matpath,'trials',sprintf('trials_bpfilt%s.mat',LogicalStr{params.bpfilt+1})),'Writable',true);
    in_trials = m.in;
    out_trials = m.out;
    V4_trials = m.V4;
end 

% Hilberting
[hilb_angles.wrapped.in]  = pre_processing_pip_hilb(in_trials);
[hilb_angles.wrapped.out]  = pre_processing_pip_hilb(out_trials);
[hilb_angles.wrapped.V4]  = pre_processing_pip_hilb(V4_trials);

% Filtering and derivating : VERY IMPORTANT: Depending on the filter, the filtered instantaneous changes is either inst_freq (in the case of sgolay) or filt_inst_freq (in the case of medfilt)    
[hilb_angles.in,inst_freq.in,filt_data.in]  = pre_processing_pip_filtinst(hilb_angles.wrapped.in,params.filttype,params.framelen,params.filtord,params.toi);
[hilb_angles.out,inst_freq.out,filt_data.out]  = pre_processing_pip_filtinst(hilb_angles.wrapped.out,params.filttype,params.framelen,params.filtord,params.toi);
[hilb_angles.V4,inst_freq.V4,filt_data.V4]  = pre_processing_pip_filtinst(hilb_angles.wrapped.V4,params.filttype,params.framelen,params.filtord,params.toi);


if params.filttype == "sgolay"     
    filt_inst_freq.in = inst_freq.in;
    filt_inst_freq.out = inst_freq.out;
    filt_inst_freq.V4 = inst_freq.V4;
elseif params.filttype == "medfilt"
    filt_inst_freq.in = filt_data.in;
    filt_inst_freq.out = filt_data.out;
    filt_inst_freq.V4 = filt_data.V4;
end 

%Saving the inst freq and angles 
filename = sprintf("inst_freq_%s_toi%.1f-%.1f_lower%i_upper%i_filttype%s_filtord%i_framelen%i.mat",params.bptype,params.toi(1),params.toi(2),params.lower,params.upper,params.filttype,params.filtord,params.framelen);
save(fullfile(params.matpath,'Inst_freq',filename),'params','-v7.3')
m = matfile(fullfile(params.matpath,'Inst_freq',filename),'Writable',true);
m.inst = filt_inst_freq;
m.angle = hilb_angles.wrapped;
m.angle_filt = filt_data;

%% Structures on phase difference between V1 and V4 starting with wrapped phase angles
filename = sprintf("inst_freq_toi%.1f-%.1f_lower%i_upper%i_filtord%i.mat",params.toi(1),params.toi(2),params.lower,params.upper,params.filtord);
m = matfile(fullfile(matpath,'Inst_freq',filename),'Writable',true);

hilb_angles.wrapped = m.angle;


% Elongating to 23 session electrode pairs
angle_long.in = elongate(hilb_angles.wrapped.in);
angle_long.out = elongate(hilb_angles.wrapped.out);


% Swapping out .trial for phase difference
angle_long_dif_wr = angle_long;
for i_s = 1:length(angle_long_dif_wr.in)
   
    angle_long_dif_wr.in(i_s).label{1,1} = [angle_long_dif_wr.in(i_s).label{1,1} '_' angle_long_dif_wr.in(i_s).label{1,2}] 
    for i_t = 1:length(angle_long_dif_wr.in(i_s).trial)
        cur_t = angle_long_dif_wr.in(i_s).trial{:,i_t};
        dif = cur_t(1,:) - cur_t(2,:);
        norm_dif = mod(dif + pi, 2 * pi) - pi;
        angle_long_dif_wr.in(i_s).trial{:,i_t}(1,:) = norm_dif;
    end 

    angle_long_dif_wr.out(i_s).label{1,1} = [angle_long_dif_wr.out(i_s).label{1,1} '_' angle_long_dif_wr.out(i_s).label{1,2}]         
    for i_t = 1:length(angle_long_dif_wr.out(i_s).trial)
        cur_t = angle_long_dif_wr.out(i_s).trial{:,i_t};
        dif = cur_t(1,:) - cur_t(2,:);
        norm_dif = mod(dif + pi, 2 * pi) - pi;        
        angle_long_dif_wr.out(i_s).trial{:,i_t}(1,:) = norm_dif;
    end 

end 
cfg = [];

% Only choosing dif channel
for ii = 1:length(angle_long_dif_wr.in)
    cur_data = angle_long_dif_wr.in(ii);
    cfg.channel = cur_data.label{1,1};
    angle_long_dif_wr.in(ii) = ft_selectdata(cfg,cur_data);
    cur_data = angle_long_dif_wr.out(ii);
    cfg.channel = cur_data.label{1,1};
    angle_long_dif_wr.out(ii) = ft_selectdata(cfg,cur_data);
end 


% Filtering and derivating : VERY IMPORTANT: Depending on the filter, the filtered instantaneous changes is either inst_freq (in the case of sgolay) or filt_inst_freq (in the case of medfilt)    
[angle_long_dif.in,inst_freq_dif.in,filt_data_dif.in]  = pre_processing_pip_filtinst(angle_long_dif_wr.in,filttype,framelen,filtord);
[angle_long_dif.out,inst_freq_dif.out,filt_data_dif.out]  = pre_processing_pip_filtinst(angle_long_dif_wr.out,filttype,framelen,filtord);

if filttype == "sgolay"    
    filt_inst_freq_dif.in = inst_freq_dif.in;
    filt_inst_freq_dif.out = inst_freq_dif.out;
elseif filttype == "medfilt"
    filt_inst_freq_dif.in = filt_data_dif.in;
    filt_inst_freq_dif.out = filt_data_dif.out;
end 

% Saving the inst freq and angles 
mt = [];
filename = sprintf("phase_dif_toi%.1f-%.1fthr%.2f_comp%d_lower%d_upper%d.mat",toi(1),toi(2),th,10,lower,upper);

save(fullfile(matpath,'phase_dif',filename),'mt','-v7.3')
m = matfile(fullfile(matpath,'phase_dif',filename),'Writable',true);
m.angle_wr = angle_long_dif_wr;
m.inst_dif = filt_inst_freq_dif;
m.angle_filt_dif = filt_data_dif;