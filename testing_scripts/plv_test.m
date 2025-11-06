clc
clear 

% Loading paramters
load("/data/projects/V1V4coherence/02_analysis_max/git_repos/params.mat")


% loading inst freq struct 
filename = sprintf("inst_freq_toi%.1f-%.1f_lower%i_upper%i_filtord%i.mat",params.toi(1),params.toi(2),params.lower,params.upper,params.filtord);
m = matfile(fullfile(params.matpath,'Inst_freq',filename),'Writable',true);
angles = m.angle;
t_length = (params.toi(2) - params.toi(1)) * 1000 + 1


angles.in = elongate(angles.in);
angles.out = elongate(angles.out);

% Calculating PC the Iris way: elongated
comb_data = angles.in;
counter = 1;
for i_sess = 1:length(comb_data)
    var = comb_data(i_sess);
    contrib = zeros(1,t_length);
    dCos = zeros(1,t_length);
    dSin = zeros(1,t_length);
    for i_trial = 1:length(var.trial)
        temp_array = NaN(1,t_length);
        temp_array(:,1:length(var.trial{i_trial})) = 1;
        phase_diff = var.trial{i_trial}(1,:) - var.trial{i_trial}(2,:);
        contrib = contrib + ~isnan(temp_array);
        dCos(:,1:length(var.trial{i_trial})) = dCos(:,1:length(var.trial{i_trial})) + cos(phase_diff);
        dSin(:,1:length(var.trial{i_trial})) = dSin(:,1:length(var.trial{i_trial})) + sin(phase_diff);
    end 
    dCos = dCos./contrib;
    dSin = dSin./contrib;
    VL_in(i_sess,:) = sqrt((dCos).^2+(dSin).^2);
    VLEXP = sqrt(pi)./(2.*sqrt(contrib));
    VL_in(i_sess,:) = VL_in(i_sess,:) - VLEXP;
    meanAngle = acos((dCos)./VL_in(i_sess,:));
    label{i_sess,1} = var.label{1,1};
    label{i_sess,2} = var.label{1,2}
end 
t_meanin = mean(VL_in,1)
sing_mean_in = mean(VL_in,2)
grand_mean_in = mean(sing_mean_in)

comb_data = angles.out;
counter = 1;
for i_sess = 1:length(comb_data)
    var = comb_data(i_sess);
    contrib = zeros(1,t_length);
    dCos = zeros(1,t_length);
    dSin = zeros(1,t_length);
    for i_trial = 1:length(var.trial)
        temp_array = NaN(1,t_length);
        temp_array(:,1:length(var.trial{i_trial})) = 1;
        phase_diff = var.trial{i_trial}(1,:) - var.trial{i_trial}(2,:);
        contrib = contrib + ~isnan(temp_array);
        dCos(:,1:length(var.trial{i_trial})) = dCos(:,1:length(var.trial{i_trial})) + cos(phase_diff);
        dSin(:,1:length(var.trial{i_trial})) = dSin(:,1:length(var.trial{i_trial})) + sin(phase_diff);
    end 
    dCos = dCos./contrib;
    dSin = dSin./contrib;
    VL_out(i_sess,:) = sqrt((dCos).^2+(dSin).^2);
    VLEXP = sqrt(pi)./(2.*sqrt(contrib));
    VL_out(i_sess,:) = VL_out(i_sess,:) - VLEXP;
    meanAngle = acos((dCos)./VL_out(i_sess,:));
    label{i_sess,1} = var.label{1,1};
    label{i_sess,2} = var.label{1,2}
end 
t_meanout = mean(VL_out,1)
sing_mean_out = mean(VL_out,2)
grand_mean_out = mean(sing_mean_out)

% Saving the structure 
PLV.label = label
PLV.in.VL = VL_in;
PLV.out.VL = VL_out;
PLV.in.t_mean = t_meanin;
PLV.out.t_mean = t_meanout;
PLV.in.s_mean = sing_mean_in;
PLV.out.s_mean = sing_mean_out;
PLV.in.g_mean = grand_mean_in;
PLV.out.g_mean = grand_mean_out;
foldername = fullfile(params.matpath,"PLV")
if ~exist(foldername,'dir')
    mkdir(foldername)
end     
filename = fullfile(params.matpath,"PLV",sprintf("PLV_toi%.1f-%.1f_bounds%i-%i.mat",params.toi(1),params.toi(2),params.lower,params.upper));
save(filename,"PLV");
%%
new = [sing_mean_in, sing_mean_out] 
[p, value] = signrank(sing_mean_in', sing_mean_out')
figure;
boxplot(new,'Labels',{'PLV attended','PLV unattended'})
ylabel('PLV')
ylim([0 0.21])
title(sprintf('PLVs per trial, SSD Boundaries: %d, %d, Time: %.1f - %.1f s',params.lower,params.upper,params.toi(1),params.toi(2)))
foldername = fullfile(params.figpath,sprintf('PLV/bounds%d-%d/toi%.1f-%.1f',params.lower,params.upper,params.toi(1),params.toi(2)));
if ~exist(foldername,'dir')
    mkdir(foldername)
end 
% saveas(gcf,fullfile(foldername,'PLV_boxplots.fig'))
% saveas(gcf,fullfile(foldername,'PLV_boxplots.jpg'))

% Plotting coherence over time
f = figure;
f.Units = "normalized";
f.Position = [0 0 0.5 0.5]
plot(t_meanin,'r');
xlabel('Time [ms]');
ylabel('PLV');
hold on
plot(t_meanout,'b');
legend('AttIn','AttOut')
title(sprintf('Average PLVs over time, SSD Boundaries: %d, %d, Time: %.1f - %.1f s',params.lower,params.upper,params.toi(1),params.toi(2)))
foldername = fullfile(params.figpath,sprintf('PLV/bounds%d-%d/toi%.1f-%.1f',params.lower,params.upper,params.toi(1),params.toi(2)));
if ~exist(foldername,'dir')
    mkdir(foldername)
end 
% saveas(gcf,fullfile(foldername,'PLV_time.fig'))
% saveas(gcf,fullfile(foldername,'PLV_time.jpg'))

%% Calculating angles using wavelet 
filename = sprintf("angle_SSDwavelet%d-%d.mat",lower,upper);
m = matfile(fullfile(matpath,'wavelet',filename),'Writable',true);
wangles.in = m.in;
wangles.out = m.out;
%%
% Calculating PC the Iris way
comb_data = wangles.in;
counter = 1;
for i_sess = 1:length(comb_data)
    var = comb_data(i_sess);
    for i_chan = 1:size(var.label,2)-1
        contrib = zeros(35,5301);
        dCos = zeros(35,5301);
        dSin = zeros(35,5301);
        for i_trial = 1:length(var.trialinfo)
            temp_array = NaN(35,5301);
            temp_array(:,1:length(var.time)) = 1;[in_SSD_trials,inc_in,testing_in] = do_SSD(in_trials,fs,th,lower,upper)

            phase_diff = squeeze(var.fourierspctrm(i_trial,1,:,:))- squeeze(var.fourierspctrm(i_trial,i_chan+1,:,:));
            contrib = contrib + ~isnan(temp_array);
            dCos = dCos + cos(phase_diff);
            dSin = dSin + sin(phase_diff);
        end 
        dCos = dCos./contrib;
        dSin = dSin./contrib;
        VL_in = sqrt((dCos).^2+(dSin).^2);
        VLEXP = sqrt(pi)./(2.*sqrt(contrib));
        VL_in = VL_in - VLEXP;
        meanAngle = acos((dCos)./VL_in);
        in_VL(counter,:,:) = VL_in;
        counter = counter + 1;
    end 
end 

comb_data = wangles.out;
counter = 1;
for i_sess = 1:length(comb_data)
    var = comb_data(i_sess);
    for i_chan = 1:size(var.label,2)-1
        contrib = zeros(35,5301);
        dCos = zeros(35,5301);
        dSin = zeros(35,5301);
        for i_trial = 1:length(var.trialinfo)
            temp_array = NaN(35,5301);
            temp_array(:,1:length(var.time)) = 1;
            phase_diff = squeeze(var.fourierspctrm(i_trial,1,:,:))- squeeze(var.fourierspctrm(i_trial,i_chan+1,:,:));
            contrib = contrib + ~isnan(temp_array);
            dCos = dCos + cos(phase_diff);
            dSin = dSin + sin(phase_diff);
        end 
        dCos = dCos./contrib;
        dSin = dSin./contrib;
        VL_in = sqrt((dCos).^2+(dSin).^2);
        VLEXP = sqrt(pi)./(2.*sqrt(contrib));
        VL_in = VL_in - VLEXP;
        meanAngle = acos((dCos)./VL_in);
        out_VL(counter,:,:) = VL_in;
        counter = counter + 1;
    end 
end 
mean_out = squeeze(mean(out_VL,1));
mean_in = squeeze(mean(in_VL,1));
%% Plotting the wavelet coherences
out = wangles.in(2);
lfax = log(out.freq');
nYtick = 17;
dyt = (lfax(end)-lfax(1))/(nYtick-1);
lyt = (0:nYtick-1)*dyt+lfax(1);
yt = exp(lyt);
lfax2 =  1:35;
dyt2   = (lfax2(end)-lfax2(1))/(nYtick-1);
lyt2   = (0:nYtick-1)*dyt2+lfax2(1);
yt2    = (lyt2);

subplot(1,2,1)
imagesc(2300:4300,1:35,mean_out(:,2300:4300))
ax = gca;
set(ax,... 
'color','none',...
'TickDir','out',...
'YDir','normal',...
'YTick',yt2,...
'YTickLabel',num2str(yt','%3.2f'));
colormap(jet(256));
colorbar();

subplot(1,2,2)
imagesc(2300:4300,1:35,mean_in(:,2300:4300))
ax = gca;
set(ax,... 
'color','none',...
'TickDir','out',...
'YDir','normal',...
'YTick',yt2,...
'YTickLabel',num2str(yt','%3.2f'));
colormap(jet(256));
colorbar();
%%

in_plv = nan(23,252); %creating nan array of all sessions and max trials
out_plv = nan(23,231); % differs from in because different max 
count = 1;
for i_sess = 1:length(angles.out)
    out_phase_diff = cell(1,length(angles.out(i_sess).trial));
    for i_chan = 1:length(angles.out(i_sess).label)-1
        for i_trial = 1:length(angles.out(i_sess).trial)
            signal_1 = angles.out(i_sess).trial{1,i_trial}(1,:);
            signal_2 = angles.out(i_sess).trial{1,i_trial}(i_chan+1,:);
            out_phase_diff{:,i_trial} = signal_1-signal_2;
            out_plv(count,i_trial) = abs(mean(exp(1j * out_phase_diff{:,i_trial})));
        end 
        count = count + 1;
    end 
end 
count = 1;
for i_sess = 1:length(angles.in)
    in_phase_diff = cell(1,length(angles.in(i_sess).trial));
    for i_chan = 1:length(angles.in(i_sess).label)-1
        for i_trial = 1:length(angles.in(i_sess).trial)
            signal_1 = angles.in(i_sess).trial{1,i_trial}(1,:);
            signal_2 = angles.in(i_sess).trial{1,i_trial}(i_chan+1,:);
            in_phase_diff{:,i_trial} = signal_1-signal_2;
            in_plv(count,i_trial) = abs(mean(exp(1j * in_phase_diff{:,i_trial})));
        end 
        count = count + 1;
    end 
end 
in_plv_m = mean(in_plv,2,'omitnan');
out_plv_m = mean(out_plv,2,'omitnan');
in_plv_fin = mean(in_plv_m);
out_plv_fin = mean(out_plv_m);


%% 
new = cat(2,in_plv_m,out_plv_m);
figure;
boxplot(new,'Labels',{'PLV attended','PLV unattended'})
ylabel('PLV')
title('PLV in the frequency range from 52 to 78 Hz using conjugate method')
ylim([0.1 0.25])
[p, value] = signrank(in_plv_m, out_plv_m)
