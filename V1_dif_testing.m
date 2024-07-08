%% Script to test the difference between the V1 conditions and correlating this to other parameters such as the correlation between V1 and V4

clc
clear 
% Parameter Settings


load("/data/projects/V1V4coherence/02_analysis_max/git_repos/params.mat")



% Loading for final analysis and taking means
filename = sprintf("inst_freq_toi%.1f-%.1f_lower%i_upper%i_filtord%i.mat",params.toi(1),params.toi(2),params.lower,params.upper,params.filtord);
m = matfile(fullfile(params.matpath,'Inst_freq',filename),'Writable',true);

inst_freq = m.inst;

inst_freq.in = elongate(inst_freq.in)
inst_freq.out = elongate(inst_freq.out)

cfg = [];
for ii = 1:length(inst_freq.in)
    cfg.channel = inst_freq.in(ii).label{1}
    inst_freq.in(ii) = ft_selectdata(cfg,inst_freq.in(ii))
    cfg.channel = inst_freq.out(ii).label{1}
    inst_freq.out(ii) = ft_selectdata(cfg,inst_freq.out(ii)) 
end 



[inst_mean.in,grand_inst.in] = do_grand_avg(inst_freq.in);
[inst_mean.out,grand_inst.out] = do_grand_avg(inst_freq.out);
cfg = [];
cfg.operation = 'x1 - x2'
cfg.parameter = 'avg'
for ii = 1:length(inst_mean.in)
    inst_mean.dif(ii) = ft_math(cfg,inst_mean.in(ii),inst_mean.out(ii))
end 
for ii = 1:length(inst_mean.in)
    inst_mean.dif(ii).t_avg = mean(inst_mean.dif(ii).avg)
end 

filename = fullfile(params.matpath,"V1_V4correlations",sprintf("V1_V4correlations_toi%.1f-%.1f_filtord%i_bounds%i-%i.mat",params.toi(1),params.toi(2),params.filtord,params.lower,params.upper))
load(filename)
% Regression detrended: V1a/V1n dif - V1/V4 correlation


X_in = [inst_mean.dif.t_avg]
Y_in = in_corr_det.sess(1,:)
X_out = [inst_mean.dif.t_avg]
Y_out = out_corr_det.sess(1,:)

mdl_in = fitlm(X_in,Y_in)
mdl_out = fitlm(X_out,Y_out)

%%
f = figure;
f.Units = 'normalized';
f.Position = [0 0 1 1]
sgtitle('Regressing difference of V1a - V1n on correlation between detrended V1 & V4')
% Regressions detrended: V1a/V1n dif - V1a/V4 corr
subplot(2,5,1:4)
plot(mdl_in)
title('V1a-V1n vs V1a/V4 corr')
xlabel('V1a - V1n difference')
ylabel('V1a - V4 correlation')

% Regressions detrended: V1a/V1n dif - V1n/V4 corr
subplot(2,5,6:9)
plot(mdl_out)
title('V1a-V1n vs V1n/V4 corr')
xlabel('V1a - V1n difference')
ylabel('V1a - V4 correlation')

annotation(f,'textbox',[0.8 0.1 0.1 0.35],'String',sprintf("R² = %.3f\npval = %.3f\nMorph Cycle: %s",mdl_in.Rsquared.Ordinary,mdl_in.Coefficients.pValue(2),num2str(params.MC)),FontSize=15)
annotation(f,'textbox',[0.8 0.575 0.1 0.35],'String',sprintf("R² = %.3f\npval = %.3f\nMorph Cycle: %s",mdl_out.Rsquared.Ordinary,mdl_out.Coefficients.pValue(2),num2str(params.MC)),FontSize=15)

% %saving
% foldername = fullfile(params.figpath,'Regressions','V1dif_V1V4cor','detrended')
% if ~exist(foldername,'dir')
%     mkdir(foldername)
% end 
% saveas(f,fullfile(foldername,sprintf('V1dif_V4cor_toi%.1f-%.1f_detrended.jpg',params.toi(1),params.toi(2))))
% saveas(f,fullfile(foldername,sprintf('V1dif_V4cor_toi%.1f-%.1f_detrended.fig',params.toi(1),params.toi(2))))
% close
%% Regression: V1a/V1n dif - V1/V4 correlation


X_in = [inst_mean.dif.t_avg]
Y_in = in_corr.sess(1,:)
X_out = [inst_mean.dif.t_avg]
Y_out = out_corr.sess(1,:)

mdl_in = fitlm(X_in,Y_in)
mdl_out = fitlm(X_out,Y_out)


f = figure;
f.Units = 'normalized';
f.Position = [0 0 1 1]
sgtitle('Regressing difference of V1a - V1n on correlation between V1 & V4')
% Regressions detrended: V1a/V1n dif - V1a/V4 corr
subplot(2,5,1:4)
plot(mdl_in)
title('V1a-V1n vs V1a/V4 corr')
xlabel('V1a - V1n difference')
ylabel('V1a - V4 correlation')

% Regressions detrended: V1a/V1n dif - V1n/V4 corr
subplot(2,5,6:9)
plot(mdl_out)
title('V1a-V1n vs V1n/V4 corr')
xlabel('V1a - V1n difference')
ylabel('V1a - V4 correlation')

annotation(f,'textbox',[0.8 0.1 0.1 0.35],'String',sprintf("R² = %.3f\npval = %.3f\nMorph Cycle: %s",mdl_in.Rsquared.Ordinary,mdl_in.Coefficients.pValue(2),num2str(params.MC)),FontSize=15)
annotation(f,'textbox',[0.8 0.575 0.1 0.35],'String',sprintf("R² = %.3f\npval = %.3f\nMorph Cycle: %s",mdl_out.Rsquared.Ordinary,mdl_out.Coefficients.pValue(2),num2str(params.MC)),FontSize=15)

%saving
foldername = fullfile(params.figpath,'Regressions','V1dif_V1V4cor')
if ~exist(foldername,'dir')
    mkdir(foldername)
end 
saveas(f,fullfile(foldername,sprintf('V1dif_V4cor_toi%.1f-%.1f.jpg',params.toi(1),params.toi(2))))
saveas(f,fullfile(foldername,sprintf('V1dif_V4cor_toi%.1f-%.1f.fig',params.toi(1),params.toi(2))))
close