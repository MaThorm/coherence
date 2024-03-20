clc
clear all
dir = "/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files";
toi = [2.3 4.3];
bp = [55 95];
medfiltord = 20;
load(fullfile(dir,sprintf("grand_structbp%d%dmed%dtoi%f-%f.mat",bp(1),bp(2),medfiltord,toi(1),toi(2))));
%%
data = grand_struct.add.out_hilbert(1);
cfg = [];
cfg.method = 'plv';
cfg.channel = 'all';
ft_connectivityanalysis(cfg, data)

%% Creating trials and bp filtering
in_plv = nan(16,252) %creating nan array of all sessions and max trials
out_plv = nan(16,231) % differs from in because different max 
for i_s = 1:length(grand_struct.add.out_hilbert)
    for ii = 1:length(grand_struct.add.out_hilbert(i_s).trial)
        signal_1 = grand_struct.add.out_hilbert(i_s).trial{1,ii}(1,:);
        signal_2 = grand_struct.add.out_hilbert(i_s).trial{1,ii}(2,:);
        out_phase_diff{:,ii} = signal_1-signal_2;
        out_plv(i_s,ii) = abs(mean(exp(1j * out_phase_diff{:,ii})));
    end 
end 
for i_s = 1:length(grand_struct.add.in_hilbert)
    for ii = 1:length(grand_struct.add.in_hilbert(i_s).trial)
        signal_1 = grand_struct.add.in_hilbert(i_s).trial{1,ii}(1,:);
        signal_2 = grand_struct.add.in_hilbert(i_s).trial{1,ii}(2,:);
        in_phase_diff{:,ii} = signal_1-signal_2;
        in_plv(i_s,ii) = abs(mean(exp(1j * in_phase_diff{:,ii})));
    end 
end 
%%
in_plv_m = mean(in_plv,2,'omitnan');
out_plv_m = mean(out_plv,2,'omitnan');
in_plv_fin = mean(in_plv_m);
out_plv_fin = mean(out_plv_m);
