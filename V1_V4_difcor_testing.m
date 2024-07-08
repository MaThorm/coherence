%% Script to test the frequency difference between V1a and V1n with V4 and correlating it with other parameters such as the correlation with V4

clc
clear 

% Loading paramters
load("/data/projects/V1V4coherence/02_analysis_max/git_repos/params.mat")


% Loading and averaging
% filename = sprintf("phase_dif_cut_toi%.1f-%.1f_sdmult%.1f_lower%d_upper%d.mat",params.toi(1),params.toi(2),params.sd_mult,params.lower,params.upper);
% m = matfile(fullfile(params.matpath,'phase_dif',filename),'Writable',true);
filename = sprintf("phase_dif_toi%.1f-%.1fthr%.2f_comp%d_lower%d_upper%d.mat",params.toi(1),params.toi(2),params.th,10,params.lower,params.upper);
m = matfile(fullfile(params.matpath,'phase_dif',filename),'Writable',true);

% inst_dif = m.inst_dif_cut;
% angle_dif = m.angle_dif_cut;
% angle_dif_wr = m.angle_dif_wr_cut;
inst_dif = m.inst_dif;
angle_dif = m.angle_filt_dif;
angle_dif_wr = m.angle_wr;

% Grand averaging
[inst_dif_mean.in,inst_dif_gmean.in] = do_grand_avg(inst_dif.in);
[inst_dif_mean.out,inst_dif_gmean.out] = do_grand_avg(inst_dif.out);
[angle_dif_mean.in,angle_dif_gmean.in] = do_grand_avg(angle_dif_wr.in);
[angle_dif_mean.out,angle_dif_gmean.out] = do_grand_avg(angle_dif_wr.out);

%changing labels back 
for ii = 1:length(inst_dif_mean.in)
    inst_dif_mean.in(ii).label = inst_dif.in(ii).label
    inst_dif_mean.out(ii).label = inst_dif.out(ii).label
    inst_dif_mean.in(ii).t_avg = mean(inst_dif_mean.in(ii).avg)
    inst_dif_mean.out(ii).t_avg = mean(inst_dif_mean.out(ii).avg)
end

filename = fullfile(params.matpath,"V1_V4correlations",sprintf("V1_V4correlations_toi%.1f-%.1f_filtord%i_bounds%i-%i.mat",params.toi(1),params.toi(2),params.filtord,params.lower,params.upper))
load(filename)
inst_dif_array.in = [inst_dif_mean.in.t_avg].*-1;
inst_dif_array.out = [inst_dif_mean.out.t_avg].*-1;

% Regressing some lines lol 
% Regressing detrended things: 
% V1a/V4 dif - V1a/V4 cor
X_in = inst_dif_array.in
Y_in = in_corr_det.sess(1,:)
mdl_in = fitlm(X_in,Y_in)
% V1n/V4 dif - V1n/V4 corr
X_out = inst_dif_array.out
Y_out = out_corr_det.sess(1,:)
mdl_out = fitlm(X_out,Y_out)

f = figure;
f.Units = 'normalized';
f.Position = [0 0 1 1]
sgtitle('Regressing V1/V4 differences to V1/V4 detrended correlations')

subplot(2,5,1:4)
plot(mdl_in)
title('V1a - V4 freq. difference vs V1a/V4 correlation')
xlabel('V1a - V4 freq difference')
ylabel('V1a/V4 correlation')
subplot(2,5,6:9)
plot(mdl_out)
title('V1n - V4 freq. difference vs V1n/V4 correlation')
xlabel('V1n - V4 freq difference')
ylabel('V1n/V4 correlation')
annotation(f,'textbox',[0.8 0.1 0.1 0.35],'String',sprintf("R² = %.3f\npval = %.3f\nMorph Cycle: %s",mdl_in.Rsquared.Ordinary,mdl_in.Coefficients.pValue(2),num2str(params.MC)),FontSize=20)
annotation(f,'textbox',[0.8 0.575 0.1 0.35],'String',sprintf("R² = %.3f\npval = %.3f\nMorph Cycle: %s",mdl_out.Rsquared.Ordinary,mdl_out.Coefficients.pValue(2),num2str(params.MC)),FontSize=20)

foldername = fullfile(params.figpath,'Regressions','V1V4dif_V1V4corr')
if ~exist(foldername,'dir')
    mkdir(foldername)
end 
saveas(f,fullfile(foldername,sprintf('V1V4dif_V1V4corr_detrended_toi%.1f-%.1f.fig',params.toi(1),params.toi(2))))
saveas(f,fullfile(foldername,sprintf('V1V4dif_V1V4corr_detrended_toi%.1f-%.1f.jpg',params.toi(1),params.toi(2))))
close all 
% Regressing non-detrended things:
% V1a/V4 dif - V1a/V4 corr
X_in = inst_dif_array.in
Y_in = in_corr.sess(1,:)
mdl_in = fitlm(X_in,Y_in)
% V1n/V4 dif - V1n/V4 corr
X_out = inst_dif_array.out
Y_out = out_corr.sess(1,:)
mdl_out = fitlm(X_out,Y_out)

f = figure;
f.Units = 'normalized';
f.Position = [0 0 1 1]
sgtitle('Regressing V1/V4 differences to V1/V4 correlations')

subplot(2,5,1:4)
plot(mdl_in)
title('V1a - V4 freq. difference vs V1a/V4 correlation')
xlabel('V1a - V4 freq difference')
ylabel('V1a/V4 correlation')
subplot(2,5,6:9)
plot(mdl_out)
title('V1n - V4 freq. difference vs V1n/V4 correlation')
xlabel('V1n - V4 freq difference')
ylabel('V1n/V4 correlation')
annotation(f,'textbox',[0.8 0.1 0.1 0.35],'String',sprintf("R² = %.3f\npval = %.3f\nMorph Cycle: %s",mdl_in.Rsquared.Ordinary,mdl_in.Coefficients.pValue(2),num2str(params.MC)),FontSize=20)
annotation(f,'textbox',[0.8 0.575 0.1 0.35],'String',sprintf("R² = %.3f\npval = %.3f\nMorph Cycle: %s",mdl_out.Rsquared.Ordinary,mdl_out.Coefficients.pValue(2),num2str(params.MC)),FontSize=20)

foldername = fullfile(params.figpath,'Regressions','V1V4dif_V1V4corr')
if ~exist(foldername,'dir')
    mkdir(foldername)
end 
saveas(f,fullfile(foldername,sprintf('V1V4dif_V1V4corr_toi%.1f-%.1f.fig',params.toi(1),params.toi(2))))
saveas(f,fullfile(foldername,sprintf('V1V4dif_V1V4corr_toi%.1f-%.1f.jpg',params.toi(1),params.toi(2))))


%% Scatterplotting frequency difference between V1/V4 versus the respective correlation 
scatter(inst_dif_array.in,in_corr.sess(1,:),'filled','r')
hold on 
scatter(inst_dif_array.out,out_corr.sess(1,:),'filled','b')
xl = xline([0],'--','color',[0.7 0.7 0.7]);
yl = yline([0],'--','color',[0.7 0.7 0.7]);
hold off
ylabel('Correllation V1/V4')
xlabel('Frequency difference V1 - V4');
legend({'V1a - V4','V1n - V4'});
