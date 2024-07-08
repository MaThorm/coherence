clc
clear 
%% loading parameters;
load("/data/projects/V1V4coherence/02_analysis_max/git_repos/params.mat")

saving = true;
% Renaming parameters because I'm too lazy to change it everywhere

% Trial preprocessing parameters
bpfilt = params.bpfilt;
bpwidth = params.bpwidth;
toi = params.toi; % Time region of interest [2.3 4.3] is MC 2&3


% SSD parameters
fs = params.fs; % Sampling frequency of the signal
th = params.th; % residual variance threshhold   
lower = params.lower;
upper = params.upper;

% Hilbert parameters
filttype = params.filttype; %either medfilt or sgolay
framelen = params.framelen;
filtord = params.filtord;

% snipping parameters
sd_mult = params.sd_mult;
outlier_mult = params.outlier_mult;

matpath = params.matpath;
figpath = params.figpath;

% Random variable setting
timebar = -1.3:0.001:5;

%% Trial selection 
% Loading
load('attout_dataset.mat')
load('attin_dataset.mat')
load('V4_dataset.mat')

% Trialselection 
[in_trials] = pre_processing_pip_trials(attin_dataset,params.bpfilt,params.bpwidth,params.toi)
[out_trials] = pre_processing_pip_trials(attout_dataset,params.bpfilt,params.bpwidth,params.toi)
[V4_trials] = pre_processing_pip_trials(V4_dataset,params.bpfilt,params.bpwidth,params.toi)

%Saving
LogicalStr = {'false', 'true'};
save(fullfile(matpath,'trials',sprintf('in_trials_bpfilt%s_toi%.1f-%.1f.mat',LogicalStr{bpfilt+1},toi(1),toi(2))),"in_trials")
save(fullfile(matpath,'trials',sprintf('out_trials_bpfilt%s_toi%.1f-%.1f.mat',LogicalStr{bpfilt+1},toi(1),toi(2))),"out_trials")
save(fullfile(matpath,'trials',sprintf('V4_trials_bpfilt%s_toi%.1f-%.1f.mat',LogicalStr{bpfilt+1},toi(1),toi(2))),"V4_trials") 

%% SSD calculation
LogicalStr = {'false', 'true'};
load(fullfile(params.matpath,'trials',sprintf('in_trials_bpfilt%s_toi%.1f-%.1f.mat',LogicalStr{params.bpfilt+1},params.toi(1),params.toi(2))));
load(fullfile(params.matpath,'trials',sprintf('out_trials_bpfilt%s_toi%.1f-%.1f.mat',LogicalStr{params.bpfilt+1},params.toi(1),params.toi(2))));
load(fullfile(params.matpath,'trials',sprintf('V4_trials_bpfilt%s_toi%.1f-%.1f.mat',LogicalStr{params.bpfilt+1},params.toi(1),params.toi(2))));


% Performing SSD
[in_SSD_trials,inc_in,testing_in] = do_SSD(in_trials,fs,th,lower,upper)
[out_SSD_trials,inc_out,testing_out] = do_SSD(out_trials,fs,th,lower,upper)
[V4_SSD_trials,inc_V4,testing_V4] = do_SSD(V4_trials,fs,th,lower,upper)

% Saving
m = matfile(fullfile(matpath,'SSD',sprintf('ssd_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))),'Writable',true);
m.in = in_SSD_trials;
m.out = out_SSD_trials;
m.V4 = V4_SSD_trials;
m = matfile(fullfile(matpath,'SSD_testing',sprintf('testing_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',th,lower,upper,10,toi(1),toi(2))),'Writable',true);
m.in = testing_in;
m.out = testing_out;
m.V4 = testing_V4;


%% Hilbert Angles, instantaneous frequency, filtered instantaneous
% frequency, instantaneous change
load("/data/projects/V1V4coherence/02_analysis_max/git_repos/params.mat")

m = matfile(fullfile(params.matpath,'SSD',sprintf('ssd_trials_th%.2f_bounds%d-%d_ncomp%d_toi%.1f-%.1f.mat',params.th,params.lower,params.upper,10,params.toi(1),params.toi(2))),'Writable',true);
in_SSD_trials = m.in;
out_SSD_trials = m.out;
V4_SSD_trials = m.V4;


% Hilberting
[hilb_angles.wrapped.in]  = pre_processing_pip_hilb(in_SSD_trials);
[hilb_angles.wrapped.out]  = pre_processing_pip_hilb(out_SSD_trials);
[hilb_angles.wrapped.V4]  = pre_processing_pip_hilb(V4_SSD_trials);

% Filtering and derivating : VERY IMPORTANT: Depending on the filter, the filtered instantaneous changes is either inst_freq (in the case of sgolay) or filt_inst_freq (in the case of medfilt)    
[hilb_angles.in,inst_freq.in,filt_data.in]  = pre_processing_pip_filtinst(hilb_angles.wrapped.in,params.filttype,params.framelen,params.filtord);
[hilb_angles.out,inst_freq.out,filt_data.out]  = pre_processing_pip_filtinst(hilb_angles.wrapped.out,params.filttype,params.framelen,params.filtord);
[hilb_angles.V4,inst_freq.V4,filt_data.V4]  = pre_processing_pip_filtinst(hilb_angles.wrapped.V4,params.filttype,params.framelen,params.filtord);

if params.filttype == "sgolay"     
    filt_inst_freq.in = inst_freq.in;
    filt_inst_freq.out = inst_freq.out;
    filt_inst_freq.V4 = inst_freq.V4
elseif params.filttype == "medfilt"
    filt_inst_freq.in = filt_data.in;
    filt_inst_freq.out = filt_data.out;
    filt_inst_freq.V4 = filt_data.V4;
end 

%Saving the inst freq and angles 
mt = [];
filename = sprintf("inst_freq_toi%.1f-%.1f_lower%i_upper%i_filtord%i.mat",params.toi(1),params.toi(2),params.lower,params.upper,params.filtord);

save(fullfile(params.matpath,'Inst_freq',filename),'mt','-v7.3')
m = matfile(fullfile(params.matpath,'Inst_freq',filename),'Writable',true);
m.inst = filt_inst_freq;
m.angle = hilb_angles.wrapped;
m.angle_filt = filt_data;


%% Snippeting 
% Loading wavelet sctruct 
filename = sprintf('trial_SSDwavelet_toi%.1f-%.1f_%d-%d.mat',toi(1),toi(2),lower,upper);
m = matfile(fullfile(matpath,"wavelet",filename))
wavelet.in = m.in;
wavelet.V4 = m.V4;
wavelet.out = m.out;

% loading inst freq struct 
filename = sprintf("inst_freq_toi%.1f-%.1fthr%.2f_comp%d_lower%d_upper%d.mat",toi(1),toi(2),th,10,lower,upper);
m = matfile(fullfile(matpath,'Inst_freq',filename),'Writable',true);
filt_data = m.inst;
angles = m.angle;
angles_filt = m.angle_filt;

% Cutting up the array
freqax = wavelet.in(1).freq;
min_ax = find(min(abs(freqax-lower))==abs(freqax-lower));
max_ax = find(min(abs(freqax-upper))==abs(freqax-upper));
yrange = [min_ax:max_ax];


[inst_cut,cut_array] = do_snippeting(wavelet,filt_data,yrange,toi*1000,sd_mult,outlier_mult);
[angle_cut,cut_array] = do_snippeting(wavelet,angles,yrange,toi*1000,sd_mult,outlier_mult);
[filt_angle_cut,cut_array] = do_snippeting(wavelet,angles_filt,yrange,toi*1000,sd_mult,outlier_mult);

% Saving
filename = sprintf("inst_freq_cut_toi%.1f-%.1f_sdmult%.1f_lower%d_upper%d.mat",toi(1),toi(2),sd_mult,lower,upper);
save(fullfile(matpath,'Inst_freq_cut',filename),'-v7.3')
m = matfile(fullfile(matpath,'Inst_freq_cut',filename),'Writable',true);
m.inst = inst_cut;
m.angle = angle_cut;
m.angle_filt = filt_angle_cut;


%% Structures on phase difference between V1 and V4 starting with wrapped phase angles
filename = sprintf("inst_freq_toi%.1f-%.1f_lower%i_upper%i_filtord%i.mat",params.toi(1),params.toi(2),params.lower,params.upper,params.filtord);
m = matfile(fullfile(matpath,'Inst_freq',filename),'Writable',true);

hilb_angles.wrapped = m.angle;


% Elongating to 23 session electrode pairs
angle_long.in = elongate(hilb_angles.wrapped.in);
angle_long.out = elongate(hilb_angles.wrapped.out);

% Tried out different way to take phase difference here, didn't work :(
% cfg = [];
% for ii = 1:length(angle_long.in)
%     cfg.channel = angle_long.in(ii).label{1,1}
%     angle_long.inV1(ii) = ft_selectdata(cfg,angle_long.in(ii));
%     cfg.channel = angle_long.in(ii).label{1,2}
%     angle_long.inV4(ii) = ft_selectdata(cfg,angle_long.in(ii));
%     cfg.channel = angle_long.out(ii).label{1,1}
%     angle_long.outV1(ii) = ft_selectdata(cfg,angle_long.out(ii));
%     cfg.channel = angle_long.out(ii).label{1,2}
%     angle_long.outV4 (ii)= ft_selectdata(cfg,angle_long.out(ii));
% end 
%    
% cfg = [];
% cfg.operation = 'subtract'
% cfg.parameter = 'trial'
% for ii = 1:length(angle_long.in)
%     angle_long.inV1V4(ii) = ft_math(cfg,angle_long.inV1(ii),angle_long.inV4(ii))
%     angle_long.outV1V4(ii) = ft_math(cfg,angle_long.outV1(ii),angle_long.outV4(ii))
% end 
    
%
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

if filttype == 'sgolay'     
    filt_inst_freq_dif.in = inst_freq_dif.in;
    filt_inst_freq_dif.out = inst_freq_dif.out;
elseif filttype == 'medfilt'
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


%% Snippeting the phase difference struct
% Loading wavelet sctruct 
filename = sprintf('trial_SSDwavelet_toi%.1f-%.1f_%d-%d.mat',toi(1),toi(2),lower,upper);
m = matfile(fullfile(matpath,"wavelet",filename))
wavelet.in = m.in;
wavelet.out = m.out;
% wavelet = rmfield(wavelet,'V4')
wavelet.in = elongate(wavelet.in)
wavelet.out = elongate(wavelet.out)

% loading inst freq struct 
filename = sprintf("phase_dif_toi%.1f-%.1fthr%.2f_comp%d_lower%d_upper%d.mat",toi(1),toi(2),th,10,lower,upper);
m = matfile(fullfile(matpath,'phase_dif',filename),'Writable',true);
angle_filt_dif  = m.angle_filt_dif;
inst_dif = m.inst_dif;
angle_dif_wr = m.angle_wr;

% Cutting up the array
freqax = wavelet.in(1).freq;
min_ax = find(min(abs(freqax-lower))==abs(freqax-lower));
max_ax = find(min(abs(freqax-upper))==abs(freqax-upper));
yrange = [min_ax:max_ax];


[angle_dif_cut,cut_array] = do_snippeting(wavelet,angle_filt_dif,yrange,toi*1000,sd_mult,outlier_mult);
[inst_dif_cut,cut_array] = do_snippeting(wavelet,inst_dif,yrange,toi*1000,sd_mult,outlier_mult);
[angle_dif_wr_cut,cut_array] = do_snippeting(wavelet,angle_dif_wr,yrange,toi*1000,sd_mult,outlier_mult);

filename = sprintf("phase_dif_cut_toi%.1f-%.1f_sdmult%.1f_lower%d_upper%d.mat",toi(1),toi(2),sd_mult,lower,upper);
mt = [];
save(fullfile(matpath,'phase_dif',filename),'mt','-v7.3')
m = matfile(fullfile(matpath,'phase_dif',filename),'Writable',true);
m.inst_dif_cut = inst_dif_cut;
m.angle_dif_cut = angle_dif_cut;
m.angle_dif_wr_cut = angle_dif_wr_cut;














%% Loading the cut freuquencies
filename = sprintf("inst_freq_cut_toi%.1f-%.1f_sdmult%.1f_lower%d_upper%d.mat",toi(1),toi(2),sd_mult,lower,upper);
m = matfile(fullfile(matpath,'Inst_freq_cut',filename))
inst_cut = m.inst;


[inst_mean_cut.in,grand_inst_cut.in] = do_grand_avg(inst_cut.in);
[inst_mean_cut.out,grand_inst_cut.out] = do_grand_avg(inst_cut.out);
[inst_mean_cut.V4,grand_inst_cut.V4] = do_grand_avg(inst_cut.V4);


%% Loading the non cut inst freq.

filename = sprintf("inst_freq_toi%.1f-%.1fthr%.2f_comp%d_lower%d_upper%d.mat",toi(1),toi(2),th,10,lower,upper);
m = matfile(fullfile(matpath,'Inst_freq',filename),'Writable',true);
filt_data = m.inst;

[inst_mean.in,grand_inst.in] = do_grand_avg(filt_data.in);
[inst_mean.out,grand_inst.out] = do_grand_avg(filt_data.out);
[inst_mean.V4,grand_inst.V4] = do_grand_avg(filt_data.V4);


%% Plotting all the cut inst. freqs together 
x = timebar(toi(1)*1000:toi(2)*1000);
f = figure;
f.Units = 'normalized'
f.Position = [0 0 0.5 0.5]
plot(x,grand_inst_cut.in.avg,'r');
hold on 
plot(x,grand_inst_cut.out.avg,'b');
title([sprintf("Summary: lower: %d - upper: %d",lower,upper) filttype sprintf("length: %d, order: %d",framelen,filtord)])
legend('AttIn','AttOut','AutoUpdate','off')
label = {'MS2','MS3','MS4'};
xl = xline([ 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
xlim([toi(1)-1.3 toi(2)-1.3])
%ylim([62 74])
hold off
xlabel('Time [s]')
ylabel('Frequency [Hz]') 
foldername = fullfile(figpath,'inst_freq_cut_16',sprintf('%d-%d_sdmult%.1f_outmult%.1f/toi%.1f-%.1f/inst_freq_plots',lower,upper,sd_mult,outlier_mult,toi(1),toi(2)))
if ~exist(foldername,'dir')
    mkdir(foldername)
end 
saveas(gcf,fullfile(foldername,sprintf("inst_freq_summary_low%d_high%d.jpg",lower,upper)));
saveas(gcf,fullfile(foldername,sprintf("inst_freq_summary_low%d_high%d.fig",lower,upper)));



%% plotting all cut together: per session 
f = figure;
f.Units = ['normalized']
f.Position = [0 0 0.5 0.5]
for ii = 1:16
    x = timebar(toi(1)*1000:toi(2)*1000);
    plot(x,inst_mean_cut.in(ii).avg,'r');
    grid on 
    hold on 
    plot(x,inst_mean_cut.out(ii).avg,'b');
    title([sprintf("Session %d: lower: %d - upper: %d, filter",ii, lower,upper) filttype sprintf("length: %d, order: %d",framelen,filtord)])
    legend('AttIn','AttOut','AutoUpdate','off')
    label = {'MS2','MS3','MS4'};
    xl = xline([ 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
    xlim([toi(1)-1.3 toi(2)-1.3])
    %ylim([62 74])
    hold off
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
    foldername = fullfile(figpath,'inst_freq_cut_16',sprintf('%d-%d_sdmult%.1f_outmult%.1f/toi%.1f-%.1f/inst_freq_plots',lower,upper,sd_mult,outlier_mult,toi(1),toi(2)))
    if ~exist(foldername,'dir')
        mkdir(foldername)
    end 
     saveas(gcf,fullfile(foldername,sprintf("inst_freq_sess%d_low%d_high%d.jpg",ii,lower,upper)));
     saveas(gcf,fullfile(foldername,sprintf("inst_freq_sess%d_low%d_high%d.fig",ii,lower,upper)));
    %w = waitforbuttonpress;
    clf 
end 

%% Plotting all the noncut inst. freqs together 
x = timebar(toi(1)*1000:toi(2)*1000);
f = figure;
f.Units = 'normalized'
f.Position = [0 0 0.5 0.5]
plot(x,grand_inst.in.avg,'r');
hold on 
plot(x,grand_inst.out.avg,'b');
title([sprintf("Summary: lower: %d - upper: %d",lower,upper) filttype sprintf("length: %d, order: %d",framelen,filtord)])
legend('AttIn','AttOut','AutoUpdate','off')
label = {'MS2','MS3','MS4'};
xl = xline([ 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
xlim([toi(1)-1.3 toi(2)-1.3])
%ylim([62 74])
hold off
xlabel('Time [s]')
ylabel('Frequency [Hz]') 
foldername = fullfile(figpath,'inst_freq',sprintf('%d-%d',lower,upper))
if ~exist(foldername,'dir')
    mkdir(foldername)
end 
% saveas(gcf,fullfile(foldername,sprintf("inst_freq_summary_low%d_high%d.jpg",lower,upper)));
% saveas(gcf,fullfile(foldername,sprintf("inst_freq_summary_low%d_high%d.fig",lower,upper)));



%% plotting all noncut together: per session 
f = figure;
f.Units = ['normalized']
f.Position = [0 0 0.5 0.5]
for ii = 1:16
    x = timebar(toi(1)*1000:toi(2)*1000);
    plot(x,inst_mean.in(ii).avg,'r');
    grid on 
    hold on 
    plot(x,inst_mean.out(ii).avg,'b');
    title([sprintf("Session %d: lower: %d - upper: %d, filter",ii, lower,upper) filttype sprintf("length: %d, order: %d",framelen,filtord)])
    legend('AttIn','AttOut','AutoUpdate','off')
    label = {'MS2','MS3','MS4'};
    xl = xline([ 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
    xlim([toi(1)-1.3 toi(2)-1.3])
    %ylim([62 74])
    hold off
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
    foldername = fullfile(figpath,'inst_freq',sprintf('%d-%d',lower,upper))
    if ~exist(foldername,'dir')
        mkdir(foldername)
    end 
    saveas(gcf,fullfile(foldername,sprintf("inst_freq_sess%d_low%d_high%d.jpg",ii,lower,upper)));
    saveas(gcf,fullfile(foldername,sprintf("inst_freq_sess%d_low%d_high%d.fig",ii,lower,upper)));
    %w = waitforbuttonpress;
    clf 
end 
close all

%% Plotting inst_freq difference 
for ii = 1:size(inst_dif_mat,1)
    plot(inst_dif_mat(ii,:));
    w = waitforbuttonpress;
    clf
end 
%% Plotting mean of inst freq differences
inst_dif_mean = mean(inst_dif_mat,1)
figure
plot(inst_dif_mean)

%%
figure
sess = 1;
for ii = 1:1000
    subplot(2,1,1)
    plot(hilb_angles.wrapped.in(sess).trial{1, ii}(1,:));
    title(sprintf('Session %d, Trial %d',sess,ii))
    subplot(2,1,2)
    plot(inst_freq.in(sess).trial{1, ii}(1,:));
    w = waitforbuttonpress;
    clf;
end 

%% Quick derivative check
a = diff(summary.in.mean)
figure
plot(a)
ylim([-0.12 0.17])

%% Plotting derivative of derivative over all recording sites
sel = 1:5000;
x = timebar(sel);
subplot(3,1,1)
ylim([-2 2])
plot(x,insummary.change_mean(sel),'r')
hold on 
title('AttIn')
subplot(3,1,2)
plot(x,outsummary.change_mean(sel),'b')
title('AttOut')
subplot(3,1,3)
plot(x,V4summary.change_mean(sel))
title('V4')

%% more derivative testing
for ii = 1:length(V4_inst)
    cfg = [];
    cfg.channel = 'all';
    cfg.derivative = 'yes';
    instchange(ii) = ft_preprocessing(cfg,V4_inst(ii));
end 
%% 
sel = 1:5000;
for ii = 1:length(V4_instchange)
    subplot(4,4,ii)
    plot(timebar(sel),V4_instchange(ii).sess_mean(sel));
end 
%% even more derivative testing
for ii = 1:length(V4_instchange)
    mean_dv{ii} = cent_diff_n(V4_inst(ii).sess_mean,1,5);
end 
maxNumCol = max(cellfun(@(c) size(c,2), mean_dv));
aMat = cell2mat(cellfun(@(c){[c nan(1,maxNumCol-numel(c))]}, mean_dv)');
colMeans = mean(aMat,1,'omitnan');
new_instchange.V4 = colMeans;

for ii = 1:length(attin_inst)
    mean_dv{ii} = cent_diff_n(attin_inst(ii).sess_mean,1,5);
end 
maxNumCol = max(cellfun(@(c) size(c,2), mean_dv));
aMat = cell2mat(cellfun(@(c){[c nan(1,maxNumCol-numel(c))]}, mean_dv)');
colMeans = mean(aMat,1,'omitnan');
new_instchange.in = colMeans;

for ii = 1:length(attout_inst)
    mean_dv{ii} = cent_diff_n(attout_inst(ii).sess_mean,1,5);
end 
maxNumCol = max(cellfun(@(c) size(c,2), mean_dv));
aMat = cell2mat(cellfun(@(c){[c nan(1,maxNumCol-numel(c))]}, mean_dv)');
colMeans = mean(aMat,1,'omitnan');
new_instchange.out = colMeans;
%
figure
subplot(3,1,1)
title('AttIn')
sel = 1:5000;
yyaxis left
plot(timebar(sel),new_instchange.in(sel));
ylabel('Change in Frequency')
hold on 
ylim([-0.15 0.15])
yyaxis right
plot(timebar(sel),insummary.mean(sel));
ylabel('Frequency [Hz]')
% 
subplot(3,1,2)
title('AttOut')

yyaxis left
plot(timebar(sel),new_instchange.out(sel));
ylabel('Change in Frequency')
hold on 
ylim([-0.15 0.15])
yyaxis right
plot(timebar(sel),outsummary.mean(sel));
ylabel('Frequency [Hz]')
%
subplot(3,1,3)
title('V4')

yyaxis left
plot(timebar(sel),new_instchange.V4(sel));
ylabel('Change in Frequency')
hold on 
ylim([-0.15 0.15])
yyaxis right
plot(timebar(sel),V4summary.mean(sel));
ylabel('Frequency [Hz]')
 
%%