%% Script to test and plot the V1/V4 correlations 
clc
clear 
% Parameter Settings


load("/data/projects/V1V4coherence/02_analysis_max/git_repos/params.mat")



% Loading for  analysis elongating
filename = sprintf("inst_freq_%s_toi%.1f-%.1f_lower%i_upper%i_filttype%s_filtord%i_framelen%i.mat",params.bptype,params.toi(1),params.toi(2),params.lower,params.upper,params.filttype,params.filtord,params.framelen);
m = matfile(fullfile(params.matpath,'Inst_freq',filename),'Writable',true);
inst_freq = m.inst;


inst_freq.in = elongate(inst_freq.in)
inst_freq.out = elongate(inst_freq.out)
inst_freq.V4 = elongate_V4(inst_freq.V4)


[inst_mean.in, grand_inst_mean.in] = do_grand_avg(inst_freq.in,true)
[inst_mean.out, grand_inst_mean.out] = do_grand_avg(inst_freq.out,true)
[inst_mean.V4, grand_inst_mean.V4] = do_grand_avg(inst_freq.V4)


% Correlation based on sessions with seperated V4
for ii = 1:length(inst_freq.in)
    [in_corr.sess(1,ii) in_corr.sess(2,ii)] = corr(inst_mean.in(ii).avg(1,:)',inst_mean.in(ii).avg(2,:)');
    [out_corr.sess(1,ii) out_corr.sess(2,ii)]  = corr(inst_mean.out(ii).avg(1,:)',inst_mean.out(ii).avg(2,:)');
end 
out_corr.sess_m = mean(out_corr.sess,2)
in_corr.sess_m = mean(in_corr.sess,2)
% % Calculating trial based correlations
% for i_s = 1:length(inst_freq.in)
%     cur_in = inst_freq.in(i_s)
%     cur_out = inst_freq.out(i_s)
%     for i_t = 1:length(cur_in.trial)
%         [in_corr.trial(i_s).r(i_t) in_corr.trial(i_s).p(i_t)] = corr(cur_in.trial{1,i_t}(1,:)',cur_in.trial{1,i_t}(2,:)');
%     end 
%     in_corr.trial(i_s).r_mean = mean(in_corr.trial(i_s).r);
%     in_corr.trial(i_s).rp_mean = mean(in_corr.trial(i_s).r(in_corr.trial(i_s).p > 0.05));
% 
%     for i_t = 1:length(cur_out.trial)
%         [out_corr.trial(i_s).r(i_t) out_corr.trial(i_s).p(i_t)] = corr(cur_out.trial{1,i_t}(1,:)',cur_out.trial{1,i_t}(2,:)');
%     end 
%     out_corr.trial(i_s).r_mean = mean(out_corr.trial(i_s).r);
% end 

%% Plotting V1a V1n, and V4 + correlation 
foldername = fullfile(params.figpath,'inst_freq',params.bptype,sprintf("bpwidth_%i-%i/toi_%.1f-%.1f/V1V4corr/sess_corr_plots",params.lower,params.upper,params.toi(1),params.toi(2)))
f = figure;
f.Units = 'normalized';
f.Position = [0.25 0.25 0.5 0.5]
for ii = 1:length(in_corr.sess)
    subplot(1,2,1)
    plot(inst_mean.in(ii).avg(1,:),'r')
    hold on 
    plot(inst_mean.in(ii).avg(2,:),'Color',"#EDB120")
    title(sprintf("V1a & V4 inst. freq, pair %i. r = %.2f, p = %.2f",ii,in_corr.sess(1,ii),in_corr.sess(2,ii)))
    legend({'AttendIn','V4'})
    xlabel('Time [ms]')
    ylabel('Frequency [Hz]')
    ylim([60 85])    
    hold off

    subplot(1,2,2)
    plot(inst_mean.out(ii).avg(1,:),'b')
    hold on 
    plot(inst_mean.out(ii).avg(2,:),'Color',"#EDB120")
    ylim([60 85])        
    title(sprintf("V1n & V4 inst. freq, pair %i. r = %.2f, p = %.2f",ii,out_corr.sess(1,ii),out_corr.sess(2,ii)))
    legend({'AttendOut','V4'})
    xlabel('Time [ms]')
    ylabel('Frequency [Hz]')
    hold off
    if ~exist(foldername,'dir')
        mkdir(foldername)
    end     
    saveas(f,fullfile(foldername,sprintf('In_out_V4_corr_sess%i.fig',ii)))
    saveas(f,fullfile(foldername,sprintf('In_out_V4_corr_sess%i.png',ii)))
%     w = waitforbuttonpress;
    clf
end 


%% Scatterplot of attend in vs attend out correlations with V4
x = -0.5:0.01:1
y = x;
f = figure; 
f.Units = 'normalized'
f.Position = [0.25 0.25 0.5 0.5]
foldername = fullfile(params.figpath,'inst_freq',params.bptype,sprintf("bpwidth_%i-%i/toi_%.1f-%.1f/V1V4corr",params.lower,params.upper,params.toi(1),params.toi(2)))
if ~exist(foldername,'dir')
    mkdir(foldername)
end 
subplot(9,1,1:8)
scatter(out_corr.sess(1,:),in_corr.sess(1,:))
hold on 
plot(x,y)
xl = xline([0],'--','color',[0.7 0.7 0.7]);
yl = yline([0],'--','color',[0.7 0.7 0.7]);
annotation(f,'textbox',[0.1 0.01 0.8 0.1],'String',sprintf("MC: %s , filttype: %s, \nAttend In cor: %.2f, Attend Out cor: %.2f", ...
    num2str(params.MC),params.bptype,in_corr.sess_m(1,1),out_corr.sess_m(1,1)))
title('In vs Out V1/V4 correlations')
ylabel('V1/V4 AttendIn')
xlabel('V1/V4 AttendOut')
saveas(gcf,fullfile(foldername,sprintf("In_vs_out_Scatter.fig",params.filtord)))
saveas(gcf,fullfile(foldername,sprintf("In_vs_out_Scatter.jpg",params.filtord)))



%% Creating new, detrended averages
inst_mean_det = inst_mean;
in_cell = struct2cell(inst_mean.in);
det_avg = cellfun(@(x) my_detrend(x),in_cell(3,:,:),'UniformOutput',false)
in_cell(3,:,:) = det_avg
inst_mean_det.in = cell2struct(in_cell,{'time','label','avg','var','dof','dimord','cfg'})

out_cell = struct2cell(inst_mean.out);
det_avg = cellfun(@(x) my_detrend(x),out_cell(3,:,:),'UniformOutput',false)
out_cell(3,:,:) = det_avg
inst_mean_det.out = cell2struct(out_cell,{'time','label','avg','var','dof','dimord','cfg'})

% Correlation 
for ii = 1:length(inst_freq.in)
    [in_corr_det.sess(1,ii) in_corr_det.sess(2,ii)] = corr(inst_mean_det.in(ii).avg(1,:)',inst_mean_det.in(ii).avg(2,:)');
    [out_corr_det.sess(1,ii) out_corr_det.sess(2,ii)]  = corr(inst_mean_det.out(ii).avg(1,:)',inst_mean_det.out(ii).avg(2,:)');
end 
in_corr_det.sess_m = mean(in_corr_det.sess,2)
out_corr_det.sess_m = mean(out_corr_det.sess,2)

%% Plotting V1a V1n, and V4 + correlation detrended
foldername = fullfile(params.figpath,'inst_freq',params.bptype,sprintf("bpwidth_%i-%i/toi_%.1f-%.1f/V1V4corr/detrended/sess_corr_plots",params.lower,params.upper,params.toi(1),params.toi(2)))
f = figure;
f.Units = 'normalized';
f.Position = [0 0 1 1]
for ii = 1:length(in_corr.sess)
    subplot(1,2,1)
    plot(inst_mean_det.in(ii).avg(1,:),'r')
    hold on 
    plot(inst_mean_det.in(ii).avg(2,:),'Color',"#EDB120")
    title(sprintf("V1a & V4 inst. freq, pair %i. r = %.2f, p = %.2f",ii,in_corr.sess(1,ii),in_corr.sess(2,ii)))
    legend({'AttendIn','V4'})
    xlabel('Time [ms]')
    ylabel('Frequency [Hz]')
 
    hold off

    subplot(1,2,2)
    plot(inst_mean_det.out(ii).avg(1,:),'b')
    hold on
    plot(inst_mean_det.out(ii).avg(2,:),'Color',"#EDB120")
    title(sprintf("V1n & V4 inst. freq, pair %i. r = %.2f, p = %.2f",ii,out_corr.sess(1,ii),out_corr.sess(2,ii)))
    legend({'AttendOut','V4'})
    xlabel('Time [ms]')
    ylabel('Frequency [Hz]')
    hold off
    if ~exist(foldername,'dir')
        mkdir(foldername)
    end     
    saveas(f,fullfile(foldername,sprintf('In_out_V4_corr_sess%i.fig',ii)))
    saveas(f,fullfile(foldername,sprintf('In_out_V4_corr_sess%i.png',ii)))
%     w = waitforbuttonpress;
%     clf
end 
%% Scatterplot of attend in vs attend out correlations with V4
x = -0.5:0.01:1
y = x;
f = figure; 
f.Units = 'normalized'
f.Position = [0.25 0.25 0.5 0.5]
foldername = fullfile(params.figpath,'inst_freq',params.bptype,sprintf("bpwidth_%i-%i/toi_%.1f-%.1f/V1V4corr/detrended",params.lower,params.upper,params.toi(1),params.toi(2)))
if ~exist(foldername,'dir')
    mkdir(foldername)
end 
subplot(9,1,1:8)
scatter(out_corr_det.sess(1,:),in_corr_det.sess(1,:))
hold on 
plot(x,y)
xl = xline([0],'--','color',[0.7 0.7 0.7]);
yl = yline([0],'--','color',[0.7 0.7 0.7]);
annotation(f,'textbox',[0.1 0.01 0.8 0.1],'String',sprintf("MC: %s , filttype: %s, \nAttend In cor: %.2f, Attend Out cor: %.2f", ...
    num2str(params.MC),params.bptype,in_corr_det.sess_m(1,1),out_corr_det.sess_m(1,1)))
title('In vs Out V1/V4 correlations')
ylabel('V1/V4 AttendIn')
xlabel('V1/V4 AttendOut')
saveas(gcf,fullfile(foldername,sprintf("In_vs_out_Scatter_filtord%i.fig",params.filtord)))
saveas(gcf,fullfile(foldername,sprintf("In_vs_out_Scatter_filtord%i.jpg",params.filtord)))

%% Saving the correlation structures
filename = fullfile(params.matpath,"V1_V4correlations",sprintf("V1_V4correlations_%s_toi%.1f-%.1f_bounds%i-%i.mat",params.bptype,params.toi(1),params.toi(2),params.lower,params.upper))
save(filename,"in_corr","in_corr_det","out_corr","out_corr_det")
%% Functions
function [new_avg] =  my_detrend(cur_avg)
    new_avg(1,:) = detrend(cur_avg(1,:));
    new_avg(2,:) = detrend(cur_avg(2,:));
end 
    
