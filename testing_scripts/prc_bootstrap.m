% Script to plot how the instantaneous V1/V4 frequency difference relates to the instantaneous V1/V4 phase difference 
clear 
load("/data/projects/V1V4coherence/02_analysis_max/git_repos/params.mat")
filename = sprintf("phase_dif_%s_toi%.1f-%.1f_lower%i_upper%i_filttype%s_filtord%i_framelen%i.mat",params.bptype,params.toi(1),params.toi(2),params.lower,params.upper,params.filttype,params.filtord,params.framelen);
m = matfile(fullfile(params.matpath,'phase_dif',filename),'Writable',true);


inst_dif = m.inst_dif;
angle_dif = m.angle_filt_dif;
angle_dif_wr = m.angle_wr;

cfg = [];
cfg.toilim = params.toi;
for ii = 1:length(angle_dif_wr.in)
    angle_dif_wr.in(ii) = ft_redefinetrial(cfg,angle_dif_wr.in(ii));
    angle_dif_wr.out(ii) = ft_redefinetrial(cfg,angle_dif_wr.out(ii));
end 


% Calculating PRC 
phase_bins = linspace(-pi,pi,13);
phase_bins = (phase_bins(1:end-1) + phase_bins(2:end)) / 2;
mean_in = calc_prc(inst_dif.in,angle_dif_wr.in,phase_bins);
mean_out = calc_prc(inst_dif.out,angle_dif_wr.out,phase_bins);




% Doing some bootstrapping 
clear mean_boot
for ii = 1:10
    new_phase = ftstr_boots(angle_dif_wr.in);
    mean_boot.in(ii,:) = calc_prc(inst_dif.in,new_phase,phase_bins);
    new_phase = ftstr_boots(angle_dif_wr.out);
    mean_boot.out(ii,:) = calc_prc(inst_dif.out,new_phase,phase_bins);
    sprintf('Num %i',ii)
end 
%

for ii = 1:size(mean_boot.in,1)
    cur_it = mean_boot.in(ii,:);
    mean_boot.in_dif(ii) = max(cur_it) - min(cur_it);
    mean_boot.in_dif = sort(mean_boot.in_dif);
    cur_it = mean_boot.out(ii,:);
    mean_boot.out_dif(ii) = max(cur_it) - min(cur_it);
    mean_boot.out_dif = sort(mean_boot.out_dif);
end 

in_dif = max(mean_in) - min(mean_in);
out_dif = max(mean_out) - min(mean_out);
mean_boot.in_len = length(mean_boot.in_dif(mean_boot.in_dif < in_dif))
mean_boot.out_len = length(mean_boot.in_dif(mean_boot.out_dif < out_dif))

% plotting prc
f = myplot()
f.Units = 'normalized'
f.Position = params.pl_l
plot(phase_bins,mean_in,'r','LineWidth',1.5)
hold on 
plot(phase_bins,mean_out,'b','LineWidth',1.5)
for ii = 1:length(mean_boot.in_dif)
    plot(phase_bins,mean_boot.in(ii,:),"Color",[1 0 0 0.005],LineWidth=0.01)
    plot(phase_bins,mean_boot.out(ii,:),"Color",[0 0 1 0.005],LineWidth=0.01)
end 
hold off
xlabel('Phase Difference [Radians]')
ylabel('Frequency Difference [Hz]')
xlim([(-pi -0.1) (pi + 0.1)])
legend({'V1a - V4','V1n - V4'},'Location','best','Interpreter','latex')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
foldername = fullfile(params.inst_figpath,"PRC")
if ~exist(foldername,'dir')
    mkdir(foldername)
end 
saveas(gcf,fullfile(foldername,"PRC.fig"))
exportgraphics(gcf,fullfile(foldername,"PRC.png"),'Resolution',600)

% Text file
filename = fullfile(params.inst_figpath,'PRC','PRC_modulation_stats.txt')
fileID = fopen(filename, 'w');
fprintf(fileID,sprintf('~~~~~~ File on the bootstrapping statistics of the PRC modulation ~~~~~~~ \n'));
fprintf(fileID,sprintf('In: Observed difference is larger than %i out of %i values. Leading to a p value of ?? \n',mean_boot.in_len,length(mean_boot.in_dif)));
fprintf(fileID,sprintf('Out: Observed difference is larger than %i out of %i values. Leading to a p value of ??',mean_boot.out_len,length(mean_boot.out_dif)));
fclose(fileID);
%% Functions 
function boots_array_new = do_bootstrap(boots_array)
% Creates a randomized version of an array, with zurückleging
    random_indices = randi(length(boots_array), size(boots_array));
    boots_array_new = boots_array(random_indices);
end 

function phase_struct = ftstr_boots(phase_struct)
    for i_s = 1:length(phase_struct)
        for i_t = 1:length(phase_struct(i_s).trial)
            phase_struct(i_s).trial{:,i_t} = do_bootstrap(phase_struct(i_s).trial{:,i_t});
        end
    end 
end 

function mean_prc = calc_prc(inst_struct,phase_struct,phase_bins)
% Function that calculates the PRC for two ft_structs, round_num should be a linspace from -pi to pi with nbins
avg_sess_array_mean = nan(length(phase_struct),length(phase_bins));

for i_s = 1:length(phase_struct)
    sess_array_mean = nan(length(phase_struct(i_s).trial),length(phase_bins));
    for i_t = 1:length(phase_struct(i_s).trial)
        cur_trial = phase_struct(i_s).trial{:,i_t};
        rounded_phases = roundToClosest(cur_trial,phase_bins);
        cur_inst = inst_struct(i_s).trial{:,i_t};
        [sess_array_mean(i_t,:) sess_array_std(i_t,:)] = averageBForA(cur_inst,rounded_phases,phase_bins);
    end 
    avg_sess_array_mean(i_s,:) = mean(sess_array_mean,1,'omitnan');
end 
mean_prc = mean(avg_sess_array_mean,1);
end 


function [avg_inst std_inst] = averageBForA(inst, phase,phase_bins)
% Gives back the average value for inst, for each phase value. Phase values should be rounded 
%     phase_bins = unique(phase);
    avg_inst = zeros(length(phase_bins),1);
    for ii = 1:length(phase_bins)
        phase_val = phase_bins(ii);
        freq_vals = inst(phase == phase_val);
        avg_inst(ii,:) = mean(freq_vals);
        std_inst(ii,:) = std(freq_vals);

    end 
end


function roundedX = roundToClosest(X, Y)
% Function that rounds values in X to closest value in Y
    roundedX = zeros(size(X));
    for i = 1:length(X)
        [~, idx] = min(abs(Y - X(i)));       
        roundedX(i) = Y(idx);
    end
end
