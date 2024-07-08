%% Script to investigate V1a/V1n frequency differences and correlating this to PLV
clear
load("/data/projects/V1V4coherence/02_analysis_max/git_repos/params.mat")

% loading structures
filename = fullfile(params.matpath,"PLV",sprintf("PLV_toi%.1f-%.1f_bounds%i-%i.mat",params.toi(1),params.toi(2),params.lower,params.upper));
load(filename)


filename = sprintf("inst_freq_toi%.1f-%.1f_lower%i_upper%i_filtord%i.mat",params.toi(1),params.toi(2),params.lower,params.upper,params.filtord);
m = matfile(fullfile(params.matpath,'Inst_freq',filename),'Writable',true);
inst_freq = m.inst;

filename = sprintf("phase_dif_toi%.1f-%.1fthr%.2f_comp%d_lower%d_upper%d.mat",params.toi(1),params.toi(2),params.th,10,params.lower,params.upper);
m = matfile(fullfile(params.matpath,'phase_dif',filename),'Writable',true);
inst_V14_dif = m.inst_dif;

inst_freq.in = elongate(inst_freq.in)
inst_freq.out = elongate(inst_freq.out)

%
cfg = [];
for ii = 1:length(inst_freq.in)
    cfg.channel = inst_freq.in(ii).label{1}
    inst_freq.in(ii) = ft_selectdata(cfg,inst_freq.in(ii))
    cfg.channel = inst_freq.out(ii).label{1}
    inst_freq.out(ii) = ft_selectdata(cfg,inst_freq.out(ii)) 
end 



[inst_mean.in,grand_inst.in] = do_grand_avg(inst_freq.in);
[inst_mean.out,grand_inst.out] = do_grand_avg(inst_freq.out);
[inst_mean.V1a4_dif,grand_inst.V1a4_dif] = do_grand_avg(inst_V14_dif.in)
[inst_mean.V1n4_dif,grand_inst.V1n4_dif] = do_grand_avg(inst_V14_dif.out)

cfg = [];
cfg.operation = 'x1 - x2'
cfg.parameter = 'avg'
for ii = 1:length(inst_mean.in)
    inst_mean.V1dif(ii) = ft_math(cfg,inst_mean.in(ii),inst_mean.out(ii))
end 
for ii = 1:length(inst_mean.in)
    inst_mean.V1dif(ii).t_avg = mean(inst_mean.V1dif(ii).avg)
    inst_mean.V1a4_dif(ii).t_avg = mean(inst_mean.V1a4_dif(ii).avg)
    inst_mean.V1n4_dif(ii).t_avg = mean(inst_mean.V1n4_dif(ii).avg)
end 


%% Testing whether V1a/V1n frequency difference explains PLVs
mdl_in = fitlm([inst_mean.V1dif.t_avg],[PLV.in.s_mean])
mdl_out = fitlm([inst_mean.V1dif.t_avg],[PLV.out.s_mean])

f = figure;
f.Units = 'normalized';
f.Position = [0 0 1 1];
sgtitle("Regression fit of V1 attended and unattended frequency difference to respective PLV values")
subplot(2,5,1:4)
plot(mdl_in)
title('V1a - V1n vs V1a - V4 PLV')
xlabel('V1a/V1n frequency difference [Hz]')
ylabel("V1a/V4 PLV")
% curxlim = ylim;
subplot(2,5,6:9)
plot(mdl_out)
title(sprintf('V1a - V1n vs V1n - V4 PLV'))
xlabel('V1a/V1n frequency difference [Hz]')
ylabel("V1n/V4 PLV")
annotation(f,'textbox',[0.8 0.1 0.1 0.35],'String',sprintf("R² = %.3f\npval = %.3f\nMorph Cycle: %s",mdl_in.Rsquared.Ordinary,mdl_in.Coefficients.pValue(2),num2str(params.MC)),FontSize=20)
annotation(f,'textbox',[0.8 0.575 0.1 0.35],'String',sprintf("R² = %.3f\npval = %.3f\nMorph Cycle: %s",mdl_out.Rsquared.Ordinary,mdl_out.Coefficients.pValue(2),num2str(params.MC)),FontSize=20)
%Saving 
foldername = fullfile(params.figpath,'Regressions/V1freqdif_PLV')
if ~exist(foldername,"dir")
    mkdir(foldername)
end 
saveas(f,fullfile(foldername,sprintf("V1freqdif_PLV_toi%.1f-%.1f.fig",params.toi(1),params.toi(2))))
saveas(f,fullfile(foldername,sprintf("V1freqdif_PLV_toi%.1f-%.1f.png",params.toi(1),params.toi(2))))


% ylim(curxlim)

%% Testing whether V1 - V4 differences explain PLVs
mdl_in = fitlm([inst_mean.V1n4_dif.t_avg],[PLV.in.s_mean])
mdl_out = fitlm([inst_mean.V1n4_dif.t_avg],[PLV.out.s_mean])

f = figure;
f.Units = 'normalized';
f.Position = [0 0 1 1];
sgtitle("Regression fit of V1/V4 frequency differences to respective PLV values")

subplot(2,5,1:4)
plot(mdl_in)
title('V1a - V4 vs V1a - V4 PLV')
xlabel('V1a/V4 frequency difference [Hz]')
ylabel("V1a/V4 PLV")
subplot(2,5,6:9)
plot(mdl_out)
title('V1n - V4 vs V1n - V4 PLV')
xlabel('V1n/V4 frequency difference [Hz]')
ylabel("V1n/V4 PLV")
annotation(f,'textbox',[0.8 0.1 0.1 0.35],'String',sprintf("R² = %.3f\npval = %.3f\nMorph Cycle: %s",mdl_in.Rsquared.Ordinary,mdl_in.Coefficients.pValue(2),num2str(params.MC)),FontSize=20)
annotation(f,'textbox',[0.8 0.575 0.1 0.35],'String',sprintf("R² = %.3f\npval = %.3f\nMorph Cycle: %s",mdl_out.Rsquared.Ordinary,mdl_out.Coefficients.pValue(2),num2str(params.MC)),FontSize=20)
%Saving 
foldername = fullfile(params.figpath,'Regressions/V1V4freqdif_PLV')
if ~exist(foldername,"dir")
    mkdir(foldername)
end 
saveas(f,fullfile(foldername,sprintf("V1V4freqdif_PLV_toi%.1f-%.1f.fig",params.toi(1),params.toi(2))))
saveas(f,fullfile(foldername,sprintf("V1V4freqdif_PLV_toi%.1f-%.1f.png",params.toi(1),params.toi(2))))