%% SCript to test the V1 freuency differences 
clear 
load("/data/projects/V1V4coherence/02_analysis_max/git_repos/params.mat")

% Loading for final analysis and taking means
filename = sprintf("inst_freq_%s_toi%.1f-%.1f_lower%i_upper%i_filttype%s_filtord%i_framelen%i.mat",params.bptype,params.toi(1),params.toi(2),params.lower,params.upper,params.filttype,params.filtord,params.framelen);
m = matfile(fullfile(params.matpath,'Inst_freq',filename),'Writable',true);

inst_freq = m.inst;
angles = m.angle;
cur_params = m.params;

[inst_mean.in,grand_inst.in] = do_grand_avg(inst_freq.in);
[inst_mean.out,grand_inst.out] = do_grand_avg(inst_freq.out);
[inst_mean.V4,grand_inst.V4] = do_grand_avg(inst_freq.V4);

for ii = 1:16
    means.in(ii) = mean(inst_mean.in(ii).avg,'omitnan');
%     means.in2(ii) = mean(inst_mean.in(ii).avg(1:1000),'omitnan');
%     means.in3(ii) = mean(inst_mean.in(ii).avg(1001:end),'omitnan');
    means.out(ii) = mean(inst_mean.out(ii).avg,'omitnan');
%     means.out2(ii) = mean(inst_mean.out(ii).avg(1:1000),'omitnan');
%     means.out3(ii) = mean(inst_mean.out(ii).avg(1001:end),'omitnan');
end 
grand_inst.in.MC23 = mean(grand_inst.in.avg,'omitnan');
grand_inst.out.MC23 = mean(grand_inst.out.avg,'omitnan');
grand_inst.V4.MC23 = mean(grand_inst.V4.avg,'omitnan');
% grand_inst.in.MC2 = mean(grand_inst.in.avg(1:1000),'omitnan');
% grand_inst.out.MC2 = mean(grand_inst.out.avg(1:1000),'omitnan');
% grand_inst.V4.MC2 = mean(grand_inst.V4.avg(1:1000),'omitnan');
% grand_inst.in.MC3 = mean(grand_inst.in.avg(1001:2001),'omitnan');
% grand_inst.out.MC3 = mean(grand_inst.out.avg(1001:2001),'omitnan');
% grand_inst.V4.MC3 = mean(grand_inst.V4.avg(1001:2001),'omitnan');

grand_inst.in.SD_MC23 = std(means.in,'omitnan');
grand_inst.out.SD_MC23 = std(means.out,'omitnan');
% grand_inst.in.SD_MC2 = std(means.in2,'omitnan');
% grand_inst.out.SD_MC2 = std(means.out2,'omitnan');
% grand_inst.in.SD_MC3 = std(means.in3,'omitnan');
% grand_inst.out.SD_MC3 = std(means.out3,'omitnan');

% Boxplot of means 
f= figure;
f.Units = 'normalized';
f.Position = [0.25 0.25 0.35 0.6]
subplot(6,1,1:5)
boxplot([means.in', means.out'],'Labels' ,{'Attended','Unattended'})
title(sprintf('Mean inst. freq from %.1f to %.1f',params.toi(1),params.toi(2)))
ylabel('Frequency [Hz]')
[p,h,stats] = signrank(means.in,means.out)
annotation('textbox', [0.13 0.08 0.78 0.1],'String',sprintf('Parameters: MC: %s , filttype: %s, bpwidth: %i - %iHz, Smooth filter: %s, Filter Window: %ims, Filter Order (sgolay): %i  \nH0 rejected: %s, p = %.4f', ...
    num2str(params.MC),params.bptype,params.lower,params.upper,params.filttype,params.framelen,params.filtord,mat2str(h),p))
foldername = fullfile(params.figpath,'inst_freq',params.bptype,sprintf("bpwidth_%i-%i/toi_%.1f-%.1f/mean_plots",params.lower,params.upper,params.toi(1),params.toi(2)))
if ~exist(foldername,'dir')
    mkdir(foldername)
end 
saveas(f,fullfile(foldername,sprintf('boxplot_mean_inst_freq.fig')))
saveas(f,fullfile(foldername,sprintf('boxplot_mean_inst_freq.jpg')))

% Scatterplot 
f = figure;
f.Units = 'normalized';
f.Position = [0.25 0.25 0.5 0.7]
subplot(6,1,1:5)

scatter(means.out,means.in)
c_max = max([means.out means.in]) + 0.2;
c_min = min([means.out means.in]) - 0.2

ylim([c_min c_max])
xlim([c_min c_max])
ylabel('V1a mean frequency [Hz]')
xlabel('V1n mean frequency [Hz]')
hold on 
plot(1:100,1:100)
title(sprintf("V1a vs V1n mean frequencies per session:"))
foldername = fullfile(params.figpath,'inst_freq',params.bptype,sprintf("bpwidth_%i-%i/toi_%.1f-%.1f/mean_plots",params.lower,params.upper,params.toi(1),params.toi(2)))
if ~exist(foldername,'dir')
    mkdir(foldername)
end
annotation('textbox', [0.13 0.08 0.78 0.1],'String',sprintf('Parameters: MC: %s , filttype: %s, bpwidth: %i - %iHz, Smooth filter: %s, Filter Window: %ims, Filter Order (sgolay): %i  \nH0 rejected: %s, p = %.4f', ...
    num2str(params.MC),params.bptype,params.lower,params.upper,params.filttype,params.framelen,params.filtord,mat2str(h),p))
saveas(f,fullfile(foldername,sprintf('scatter_mean_inst_freq.fig')))
saveas(f,fullfile(foldername,sprintf('scatter_mean_inst_freq.jpg')))

close all
%%
f= figure;
f.Units = 'normalized'
f.Position = [0 0 0.5 0.5]

subplot(1,3,1)
boxplot([means.in', means.out'],'Labels' ,{'Attended','Unattended'})
title('Mean inst. freq over MC 2 & 3')

% subplot(1,3,2)
% boxplot([means.in2',means.out2'],'Labels' ,{'Attended','Unattended'})
% ylim([63 76])
% title('Mean inst. freq over MC 2')
% 
% subplot(1,3,3)
% boxplot([means.in3',means.out3'],'Labels' ,{'Attended','Unattended'})
% ylim([63 76])
% title('Mean inst. freq over MC 3')


%% Same with violins
f= figure;
f.Units = 'normalized'
f.Position = [0 0 0.5 0.5]

subplot(1,3,1)
violin([means.in', means.out'],'Labels' ,{'Attended','Unattended'})
title('Mean inst. freq over MC 2 & 3')

% subplot(1,3,2)
% violin([means.in2',means.out2'],'Labels' ,{'Attended','Unattended'})
% title('Mean inst. freq over MC 2')
% 
% subplot(1,3,3)
% violin([means.in3',means.out3'],'Labels' ,{'Attended','Unattended'})
% title('Mean inst. freq over MC 3')

%% Statistics of means 
[p,h,stats] = signrank(means.in,means.out)
[p,h,stats] = signrank(means.in2,means.out2)
[p,h,stats] = signrank(means.in3,means.out3)

%% Statistics of means 
[h,p] = ttest(means.in,means.out)
[h,p] = ttest(means.in2,means.out2)
[h,p] = ttest(means.in3,means.out3)

%% Taking trials as Ns 
incell = ft_cell_cat(inst_freq.in);
outcell = ft_cell_cat(inst_freq.out);
trial_mean.in23 = cell2mat(cellfun(@(x) mean(x(1,:),2,'omitnan'),incell,'UniformOutput',false))
trial_mean.out23 = cell2mat(cellfun(@(x) mean(x(1,:),2,'omitnan'),outcell,'UniformOutput',false))
trial_mean.in2 = cell2mat(cellfun(@(x) mean(x(1,1:992),2,'omitnan'),incell,'UniformOutput',false))
trial_mean.out2 = cell2mat(cellfun(@(x) mean(x(1,1:992),2,'omitnan'),outcell,'UniformOutput',false))
trial_mean.in3 = cell2mat(cellfun(@(x) mean(x(1,1001:end),2,'omitnan'),incell,'UniformOutput',false))
trial_mean.out3 = cell2mat(cellfun(@(x) mean(x(1,1001:end),2,'omitnan'),outcell,'UniformOutput',false))

trial_SD.in23 = cell2mat(cellfun(@(x) std(x(1,:)','omitnan'),incell,'UniformOutput',false));
trial_SD.out23 = cell2mat(cellfun(@(x) std(x(1,:)','omitnan'),outcell,'UniformOutput',false));
trial_SD.in2 = cell2mat(cellfun(@(x) std(x(1,1:992)','omitnan'),incell,'UniformOutput',false));
trial_SD.out2 = cell2mat(cellfun(@(x) std(x(1,1:992)','omitnan'),outcell,'UniformOutput',false));
trial_SD.in3 = cell2mat(cellfun(@(x) std(x(1,1001:end)','omitnan'),incell,'UniformOutput',false));
trial_SD.out3 = cell2mat(cellfun(@(x) std(x(1,1001:end)','omitnan'),outcell,'UniformOutput',false));
%% Boxplots 
g1 = repmat({'In'},length(trial_mean.in23),1);
g2 = repmat({'Out'},length(trial_mean.out23),1);
g = [g1;g2]
f = figure;
f.Units = 'normalized';
f.Position = [0 0 0.5 0.5];

subplot(1,3,1)
boxplot([trial_mean.in23 trial_mean.out23],g)
title('MC 2 & 3')

subplot(1,3,2)
boxplot([trial_mean.in2 trial_mean.out2],g)
title('MC 2')

subplot(1,3,3)
boxplot([trial_mean.in3 trial_mean.out3],g)
title('MC 3')


%% Stats
[h,p] = ttest2(trial_mean.in23,trial_mean.out23)
[h,p] = ttest2(trial_mean.in2,trial_mean.out2)
[h,p] = ttest2(trial_mean.in3,trial_mean.out3)

[p,h,stats] = ranksum(trial_mean.in23,trial_mean.out23)
[p,h,stats] = ranksum(trial_mean.in2,trial_mean.out2)
[p,h,stats] = ranksum(trial_mean.in3,trial_mean.out3)

%% Correlation testing between V1a and V1n FI
% Cannot do between individual session so has to be done on mean
f = figure 
f.Units = 'normalized';
f.Position = [0.25 0.25 0.5 0.5]
for ii = 1:length(inst_mean.in)
    subplot(10,1,1:9)
    plot(inst_mean.in(ii).avg,'r')
    hold on 
    plot(inst_mean.out(ii).avg,'b')
    ylabel('Frequency [Hz]')
    xlabel('Time [ms]')
    title(sprintf('Inst. freq and correlation of V1a and V1n: Site %i',ii));
    ylim([60 80])
    hold off 
    [R,P] = corr(inst_mean.in(ii).avg',inst_mean.out(ii).avg')
    annotation('textbox', [0.13 0.05 0.2 0.1],'String',sprintf('R = %.2f, p = %.2f',R,P))
    foldername = fullfile(params.figpath,'corr',params.bptype,sprintf("bpwidth_%i-%i/toi_%.1f-%.1f/V1a_V1ncorr",params.lower,params.upper,params.toi(1),params.toi(2)))
    if ~exist(foldername,'dir')
        mkdir(foldername)
    end 
    saveas(f,fullfile(foldername,sprintf('correlations_recsite%i.fig',ii)))
    saveas(f,fullfile(foldername,sprintf('correlations_recsite%i.png',ii)))
%     w = waitforbuttonpress;
    clf
end 

%% Correlation of grand mean 
f = figure;
f.Units = 'normalized';
f.Position = [0.25 0.25 0.5 0.5] 
plot(grand_inst.in.avg,'r')
hold on 
plot(grand_inst.out.avg,'b')
hold off
xlabel('Time [ms]')
ylabel('Frequency [Hz]')
ylim([65 75])

%%
for ii = 1:length(inst_mean.in)
    a = inst_mean.in(ii).avg;
    b = inst_mean.out(ii).avg;
    [c, lags] = xcorr(a,b)
    stem(lags,c)
%     ylim([max(c)-1e5,max(c)])
    xline(0,'LineWidth',2)
    w = waitforbuttonpress;
    clf
end 