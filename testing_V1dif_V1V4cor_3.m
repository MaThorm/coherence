%% Script to test the frequency difference between V1a and V1n and correlating it with other parameters such as the correlation with V4

clc
clear 

% Loading 
load("/data/projects/V1V4coherence/02_analysis_max/git_repos/params.mat")
filename = sprintf("inst_freq_%s_toi%.1f-%.1f_lower%i_upper%i_filttype%s_filtord%i_framelen%i.mat",params.bptype,params.toi(1),params.toi(2),params.lower,params.upper,params.filttype,params.filtord,params.framelen);
m = matfile(fullfile(params.matpath,'Inst_freq',filename),'Writable',true);
inst_freq = m.inst;
angles = m.angle;
cur_params = m.params;

filename = fullfile(params.matpath,"V1_V4correlations",sprintf("V1_V4correlations_%s_toi%.1f-%.1f_bounds%i-%i.mat",params.bptype,params.toi(1),params.toi(2),params.lower,params.upper))
load(filename)

filename = fullfile(params.matpath,"PLV",sprintf("PLV_%s_toi%.1f-%.1f_bounds%i-%i.mat",params.bptype,params.toi(1),params.toi(2),params.lower,params.upper));
load(filename)


inst_freq.in = elongate(inst_freq.in);
inst_freq.out = elongate(inst_freq.out);
inst_freq.V4 = elongate_V4(inst_freq.V4)

[inst_mean.in,grand_inst.in] = do_grand_avg(inst_freq.in,true);
[inst_mean.out,grand_inst.out] = do_grand_avg(inst_freq.out,true);
[inst_mean.V4,grand_inst.V4] = do_grand_avg(inst_freq.V4,true);

%selecting only first channel
for ii = 1:length(inst_mean.in)
    cfg = [];
    cfg.channel = inst_mean.in(ii).label{1,1};
    inst_mean.in(ii) = ft_selectdata(cfg,inst_mean.in(ii))
    cfg = [];
    cfg.channel = inst_mean.out(ii).label{1,1};
    inst_mean.out(ii) = ft_selectdata(cfg,inst_mean.out(ii))    
end 

% Calculating inst_freq difference 
for ii = 1:length(inst_mean.in)
    cfg = [];
    cfg.operation = 'x1-x2';
    cfg.parameter = 'avg';
    inst_mean.V1dif(ii) = ft_math(cfg,inst_mean.in(ii),inst_mean.out(ii));
end 
cfg = [];
cfg.parameter = 'avg';
grand_inst.V1dif = ft_timelockgrandaverage(cfg,inst_mean.V1dif(1),inst_mean.V1dif(2),inst_mean.V1dif(3),inst_mean.V1dif(4),inst_mean.V1dif(5),inst_mean.V1dif(6),inst_mean.V1dif(7),inst_mean.V1dif(8),inst_mean.V1dif(9),inst_mean.V1dif(10),inst_mean.V1dif(11),inst_mean.V1dif(12),inst_mean.V1dif(13),inst_mean.V1dif(14),inst_mean.V1dif(15),inst_mean.V1dif(16),inst_mean.V1dif(17),inst_mean.V1dif(18),inst_mean.V1dif(19),inst_mean.V1dif(20),inst_mean.V1dif(21),inst_mean.V1dif(22),inst_mean.V1dif(23))

for ii = 1:length(inst_mean.V1dif)
    inst_mean.V1dif(ii).t_avg = mean(inst_mean.V1dif(ii).avg);
end 

% Calculating correlations between V1 differences and V1/V4 correlations
[in_dif_corr(1), in_dif_corr(2)] = corr([inst_mean.V1dif.t_avg]',in_corr.sess(1,:)');
[out_dif_corr(1), out_dif_corr(2)] = corr([inst_mean.V1dif.t_avg]',out_corr.sess(1,:)');

% same for detrended
[in_dif_corr_det(1), in_dif_corr_det(2)] = corr([inst_mean.V1dif.t_avg]',in_corr_det.sess(1,:)');
[out_dif_corr_det(1), out_dif_corr_det(2)] = corr([inst_mean.V1dif.t_avg]',out_corr_det.sess(1,:)');

[in_dif_PLV(1), in_dif_PLV(2)] = corr([inst_mean.V1dif.t_avg]',PLV.in.s_mean);
[out_dif_PLV(1), out_dif_PLV(2)] = corr([inst_mean.V1dif.t_avg]',PLV.out.s_mean);

% Scatterplotting v1dif vs in_corr
f = figure;
f.Units = 'normalized';
f.Position = [0.25 0.25 0.75 0.5]

% Scatterplotting in correaltion 
subplot(1,9,1:4)
scatter([inst_mean.V1dif.t_avg],in_corr.sess(1,:),70,'filled','MarkerFaceColor','red')
xlabel('V1a - V1n difference [Hz]')
ylabel('V1a/V4 correlation')
title(sprintf('V1a/V1n difference vs V1a/V4 correlation: r = %.2f, p = %.2f',in_dif_corr(1),in_dif_corr(2)))

% Scatterplotting out correlation 
subplot(1,9,5:8)
scatter([inst_mean.V1dif.t_avg],out_corr.sess(1,:),70,'filled','MarkerFaceColor','blue')
xlabel('V1a - V1n difference [Hz]')
ylabel('V1n/V4 correlation')
title(sprintf('V1a/V1n difference vs V1n/V4 correlation: r = %.2f, p = %.2f',out_dif_corr(1),out_dif_corr(2)))

t = annotation("textbox",[0.85 0.11 0.12 0.83],'String',sprintf('Parameters:\nBPtype: %s BPbounds: %i-%iHz\nMCs: %s \nSmoothing: %s \nSmooth Window: %ims\nFilt order (sgolay): %i',params.bptype,params.lower,params.upper,num2str(params.MC),params.filttype,params.framelen,params.filtord))


foldername = fullfile(params.figpath,'inst_freq',params.bptype,sprintf("bpwidth_%i-%i/toi_%.1f-%.1f/V1dif_V1V4corr",params.lower,params.upper,params.toi(1),params.toi(2)))
if ~exist(foldername,'dir')
    mkdir(foldername)
end 
saveas(f,fullfile(foldername,'V1dif_V1V4corr_scatter.fig'))
saveas(f,fullfile(foldername,'V1dif_V1V4corr_scatter.jpg'))

% Scatterplotting v1dif vs in_corr: dettrended
f = figure;
f.Units = 'normalized';
f.Position = [0.25 0.25 0.75 0.5]

% Scatterplotting in correaltion 
subplot(1,9,1:4)
scatter([inst_mean.V1dif.t_avg],in_corr_det.sess(1,:),70,'filled','MarkerFaceColor','red')
xlabel('V1a - V1n difference')
ylabel('V1a/V4 correlation (detrended)')
title(sprintf('V1a/V1n diff vs detrended V1a/V4 cor: r = %.2f, p = %.2f',in_dif_corr_det(1),in_dif_corr_det(2)))

% Scatterplotting out correlation 
subplot(1,9,5:8)
scatter([inst_mean.V1dif.t_avg],out_corr_det.sess(1,:),70,'filled','MarkerFaceColor','blue')
xlabel('V1a - V1n difference')
ylabel('V1n/V4 correlation (detrended)')
title(sprintf('V1a/V1n diff vs detrended V1n/V4 cor: r = %.2f, p = %.2f',out_dif_corr_det(1),out_dif_corr_det(2)))

t = annotation("textbox",[0.85 0.11 0.12 0.83],'String',sprintf('Parameters:\nBPtype: %s BPbounds: %i-%iHz\nMCs: %s \nSmoothing: %s \nSmooth Window: %ims\nFilt order (sgolay): %i',params.bptype,params.lower,params.upper,num2str(params.MC),params.filttype,params.framelen,params.filtord))


foldername = fullfile(params.figpath,'inst_freq',params.bptype,sprintf("bpwidth_%i-%i/toi_%.1f-%.1f/V1dif_V1V4corr/detrended",params.lower,params.upper,params.toi(1),params.toi(2)))
if ~exist(foldername,'dir')
    mkdir(foldername)
end 
saveas(f,fullfile(foldername,'V1dif_V1V4corrdet_scatter.fig'))
saveas(f,fullfile(foldername,'V1dif_V1V4corrdet_scatter.jpg'))


% V1dif vs PLV
%%%Scatterplotting 
f = figure;
f.Units = 'normalized';
f.Position = [0.25 0.25 0.75 0.5]

% Scatterplotting in PLV 
subplot(1,9,1:4)
scatter([inst_mean.V1dif.t_avg],PLV.in.s_mean,70,'filled','MarkerFaceColor','red')
xlabel('V1a - V1n difference')
ylabel('V1a/V4 PLV')
title(sprintf('V1a/V1n difference vs V1a/V4 PLV: r = %.2f, p = %.2f',in_dif_PLV(1),in_dif_PLV(2)))

% Scatterplotting out PLV 
subplot(1,9,5:8)
scatter([inst_mean.V1dif.t_avg],PLV.out.s_mean,70,'filled','MarkerFaceColor','blue')
xlabel('V1a - V1n difference')
ylabel('V1a/V4 PLV')
title(sprintf('V1a/V1n difference vs V1a/V4 PLV: r = %.2f, p = %.2f',out_dif_PLV(1),out_dif_PLV(2)))

t = annotation("textbox",[0.85 0.11 0.12 0.83],'String',sprintf('Parameters:\nBPtype: %s BPbounds: %i-%iHz\nMCs: %s \nSmoothing: %s \nSmooth Window: %ims\nFilt order (sgolay): %i',params.bptype,params.lower,params.upper,num2str(params.MC),params.filttype,params.framelen,params.filtord))

foldername = fullfile(params.figpath,'inst_freq',params.bptype,sprintf("bpwidth_%i-%i/toi_%.1f-%.1f/V1dif_V1V4PLV/",params.lower,params.upper,params.toi(1),params.toi(2)))
if ~exist(foldername,'dir')
    mkdir(foldername)
end 
saveas(f,fullfile(foldername,'V1dif_V1V4PLV_scatter.fig'))
saveas(f,fullfile(foldername,'V1dif_V1V4PLV_scatter.jpg'))

close all