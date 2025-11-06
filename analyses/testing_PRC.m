

%% Script to plot how the instantaneous V1/V4 frequency difference relates to the instantaneous V1/V4 phase difference 
clear 
load("/data/projects/V1V4coherence/02_analysis_max/git_repos/params.mat")
filename = sprintf("phase_dif_%s_toi%.1f-%.1f_lower%i_upper%i_filttype%s_filtord%i_framelen%i.mat",params.bptype,params.toi(1),params.toi(2),params.lower,params.upper,params.filttype,params.filtord,params.framelen);
m = matfile(fullfile(params.matpath,'phase_dif',filename),'Writable',true);


inst_dif = m.inst_dif;
angle_dif = m.angle_filt_dif;
angle_dif_wr = m.angle_wr;

% cfg = 
% inst_dif.in = 

% Grand averaging
[inst_dif_mean.in,inst_dif_gmean.in] = do_grand_avg(inst_dif.in);
[inst_dif_mean.out,inst_dif_gmean.out] = do_grand_avg(inst_dif.out);
[angle_dif_mean.in,angle_dif_gmean.in] = do_grand_avg(angle_dif_wr.in);
[angle_dif_mean.out,angle_dif_gmean.out] = do_grand_avg(angle_dif_wr.out);

%changing labels back 
for ii = 1:length(inst_dif_mean.in)
    inst_dif_mean.in(ii).label = inst_dif.in(ii).label
    inst_dif_mean.out(ii).label = inst_dif.out(ii).label
end

% 
cfg = [];
cfg.toilim = params.toi;
for ii = 1:length(angle_dif_wr.in)
    angle_dif_wr.in(ii) = ft_redefinetrial(cfg,angle_dif_wr.in(ii));
    angle_dif_wr.out(ii) = ft_redefinetrial(cfg,angle_dif_wr.out(ii));
end 
    
%% Getting the inst. freq dif as a function of inst. phase
[inst_dif_perphase.in,phase_dif_key.in] = process_condition(angle_dif_wr.in, inst_dif.in);
[inst_dif_perphase.out,phase_dif_key.out] = process_condition(angle_dif_wr.out, inst_dif.out);
sz = [length(inst_dif_perphase.in{1}),length(inst_dif_perphase.in)];
inst_dif_perphase.in = cell2mat(inst_dif_perphase.in);
inst_dif_perphase.out = cell2mat(inst_dif_perphase.out);
inst_dif_perphase.in = reshape(inst_dif_perphase.in,sz);
inst_dif_perphase.out = reshape(inst_dif_perphase.out,sz);

%% Plotting inst. freq difference as a function of inst. phase per session
for ii = 1:size(inst_dif_perphase.in,1)
    plot(phase_dif_key.in{1},inst_dif_perphase.in(:,ii),'r')
    hold on 
    plot(phase_dif_key.out{1},inst_dif_perphase.out(:,ii),'b')
    title(sprintf('Inst. freq. difference as a function of inst. phase difference per condition: Pair %i',ii))
    legend({'AttendIn','AttendOut'})
    w = waitforbuttonpress
    clf
end 

%% PLotting inst. freq difference as a function of inst. phase grand mean 
f = figure
f.Units = 'normalized'
f.Position = [0.25 0.25 0.5 0.5]
plot(phase_dif_key.in{1},mean(inst_dif_perphase.in,2),'r')
hold on 
plot(phase_dif_key.out{1},mean(inst_dif_perphase.out,2),'b')
xlabel('Phase')
ylabel('Frequency [Hz]')
title(sprintf('Inst. freq. difference as a function of inst. phase difference per condition: toi %.1f-%.1f',params.toi(1),params.toi(2)));
legend({'V1a/V4dif','V1n/V4dif'})
foldername = fullfile(params.figpath,'inst_freq',params.bptype,sprintf("bpwidth_%i-%i/toi_%.1f-%.1f/PRC",params.lower,params.upper,params.toi(1),params.toi(2)))
if ~exist(foldername,'dir')
    mkdir(foldername)
end 
saveas(gcf,fullfile(foldername,"PRC.fig"))
saveas(gcf,fullfile(foldername,"PRC.png"))

%% Inst dif per sess
f = figure
f.Units = 'normalized'
f.Position = [0.25 0.25 0.5 0.5]
for ii = 1:length(inst_dif_mean.in)
    plot(inst_dif_mean.in(ii).avg)
    hold on 
    plot(inst_dif_mean.out(ii).avg)
    title(sprintf('Inst frequency difference, pair %i, electrodes %s',ii,inst_dif_mean.in(ii).label{1}))
    hold off
    w = waitforbuttonpress
    clf 
end 
%% Inst dif grand average
figure
plot(inst_dif_gmean.in.avg,'r')
hold on
plot(inst_dif_gmean.out.avg,'b')
hold off





%%

%% Scatterplotting them against each other
f = figure
f.Units = 'normalized'
f.Position = [0.25 0.25 0.5 0.5]
for ii = 1:length(inst_dif.in_cell)
    scatter(angle_dif_wr.in_cell{1,ii},inst_dif.in_cell{1,ii},0.6)
    w = waitforbuttonpress;
    clf 
end 





















%%
for i_t = 1:length(angles_in)
    dif_cell = struct2cell(angles_in(i_t));
    dif_cell = dif_cell(11,:);
    dif_arr = cell2mat(dif_cell{1,1});
    dif_arr = reshape(dif_arr,[length(dif_cell{1,1}) params.tlength]);
    angles_in(i_t).dif_avg = mean(dif_arr,1,'omitnan');

    dif_cell = struct2cell(angles_out(i_t));
    dif_cell = dif_cell(11,:);
    dif_arr = cell2mat(dif_cell{1});
    dif_arr = reshape(dif_arr,[length(dif_cell{1,1}) params.tlength]);
    angles_out(i_t).dif_avg = mean(dif_arr,1,'omitnan');
end 

%% Plotting avg difference per session 
f = figure;
f.Units = 'normalized';
f.Position = [0.25 0.25 0.5 0.5]
for i_t = 1:length(angles_in)
    plot(angles_in(i_t).dif_avg,'r');
    hold on 
    plot(angles_out(i_t).dif_avg,'b');
    xlabel('Time [ms]');
    ylabel('Phase Angle difference')
    w = waitforbuttonpress
    clf
end 

%% Plotting the total difference 
grand_dif = struct2cell(angles_in);
grand_dif = squeeze(grand_dif(12,:,:));
grand_dif = cell2mat(grand_dif);
grand_dif_in = mean(grand_dif,1,'omitnan');

grand_dif = struct2cell(angles_out);
grand_dif = squeeze(grand_dif(12,:,:));
grand_dif = cell2mat(grand_dif);
grand_dif_out = mean(grand_dif,1,'omitnan');

plot(grand_dif_in,'r');
hold on
plot(grand_dif_out,'b');
hold off

sess_num = 1;
cur_sess = angles_in(sess_num);
f = figure;
for ii = 1:length(cur_sess.dif)
    plot(cur_sess.dif{1,ii})
     title(sprintf('Session %i, trial %i',sess_num,ii))
    w = waitforbuttonpress
    clf 
end 


%% Just plotting V1 and V4 inst. freqs
filename = sprintf("inst_freq_cut_toi%.1f-%.1f_sdmult%.1f_lower%d_upper%d.mat",params.toi(1),params.toi(2),params.sd_mult,params.lower,params.upper);
m = matfile(fullfile(params.matpath,'Inst_freq_cut',filename),'Writable',true);
inst_freq = m.inst;
angles = m.angle;

[inst_mean.in,grand_inst.in] = do_grand_avg(inst_freq.in);
[inst_mean.out,grand_inst.out] = do_grand_avg(inst_freq.out);
[inst_mean.V4,grand_inst.V4] = do_grand_avg(inst_freq.V4);
%%
figure;
plot(grand_inst.in.avg,'r');
hold on 
plot(grand_inst.out.avg,'b');
plot(grand_inst.V4.avg,'y');

%% Function to pad arrays with nans so I can take the mean with cell2mat
function trial = padwithnan(trial,tlength)
    trial(end+1:tlength) = nan;
end



function [new_mean_cell,all_keys] = process_condition(angle_dif_wr, inst_dif)
    % Function that takes all phase angle difference values and all inst.freq
    % difference values from a single session as to allow scatterplotting these
    % next to each other. Next it averages over phase values. 
    for i_s = 1:length(angle_dif_wr)
        cur = [];
        for i_t = 1:length(angle_dif_wr(i_s).trial)
            cur = [cur angle_dif_wr(i_s).trial{1,i_t}];
        end 
        angle_dif_cell{i_s} = cur;
        cur = [];
        for i_t = 1:length(inst_dif(i_s).trial)
            cur = [cur inst_dif(i_s).trial{1,i_t}];
        end 
        inst_dif_cell{i_s} = cur;
    end 
    
    inst_dif_redcell = cellfun(@(x) x(~isnan(x)),inst_dif_cell,'UniformOutput',false);
    angle_dif_redcell = cellfun(@(x) x(~isnan(x)),angle_dif_cell,'UniformOutput',false);
%     angle_dif_redcell = cellfun(@(x) round(x,1),angle_dif_redcell,'UniformOutput',false)
%   Rounding to next quarter integer
    angle_dif_redcell = cellfun(@(x) round(x * 4)/4,angle_dif_redcell,'UniformOutput',false)
    for i_s = 1:length(angle_dif_redcell)

        groupedValues = containers.Map('KeyType', 'double', 'ValueType', 'any');
        for i = 1:length(angle_dif_redcell{1,i_s})
%             fprintf(sprintf('\n i_s = %i, i = %i',i_s,i))
            key = angle_dif_redcell{1,i_s}(1,i);
            value = inst_dif_redcell{1,i_s}(1,i);
            
            % If the key already exists in the map, append the value to the list
            if isKey(groupedValues, key)
                groupedValues(key) = [groupedValues(key), value];
            else
                % If the key does not exist, create a new entry with the value
                groupedValues(key) = value;
            end
        end
        
        % Display the grouped values
        keys = groupedValues.keys;
        means = [];
        for i = 1:length(keys)
            key = keys{i};
            means(i) = mean(groupedValues(key));
        end
        all_keys{i_s,:} = cell2mat(keys);
        all_means{i_s,:} = means;
    
    end 
    % Checking if all keys are the same and if not pad with nans
    key_lengths = cellfun(@length, all_keys);
    
    % Find longest set of keys, set as reference
    [max_length, max_key_pos] = max(key_lengths);
    reference_keys = all_keys{max_key_pos};
    
    new_mean_cell = cell(1, length(all_keys));
    fprintf('Oh yeah\n')
    for ii = 1:length(all_keys)
        new_mean = nan(1, max_length);
        [is_present, loc] = ismember(reference_keys, all_keys{ii});   
        % Fill new_mean with corresponding values from all_means
        new_mean(is_present) = all_means{ii};
        new_mean_cell{ii} = new_mean;
    end
end 
