%% Script to do analyis on the PLVs 
filename = fullfile(params.matpath,"PLV",sprintf("PLV_%s_toi%.1f-%.1f_bounds%i-%i.mat",params.bptype,params.toi(1),params.toi(2),params.lower,params.upper));
load(filename)

% Boxplots of PLV
new = [PLV.in.s_mean, PLV.out.s_mean] 
[p, value] = signrank(PLV.in.s_mean', PLV.out.s_mean')
f = figure;
f.Units = 'normalized';
f.Position = [0.25 0.25 0.6 0.7];

subplot(1,9,1:8)
boxplot(new,'Labels',{'PLV attended','PLV unattended'})
ylabel('PLV')
ylim([0 0.21])
title(sprintf('PLVs per trial',params.bptype,params.lower,params.upper,num2str(params.MC)))

bool_vec = {'false','true'}
t = annotation("textbox",[0.85 0.11 0.12 0.83],'String',sprintf('Parameters:\nBPtype: %s BPbounds: %i-%iHz\nMCs: %s \nSmoothing: %s \nSmooth Window: %ims\nFilt order (sgolay): %i \nSign. Diff: %s',params.bptype,params.lower,params.upper,num2str(params.MC),params.filttype,params.framelen,params.filtord,bool_vec{value+1}))

foldername = fullfile(params.figpath,'inst_freq',params.bptype,sprintf("bpwidth_%i-%i/toi_%.1f-%.1f/PLV",params.lower,params.upper,params.toi(1),params.toi(2)))
if ~exist(foldername,'dir')
    mkdir(foldername)
end 
saveas(gcf,fullfile(foldername,'PLV_boxplots.fig'))
saveas(gcf,fullfile(foldername,'PLV_boxplots.jpg'))

% Scatterplots of PLV
f = figure;
f.Units = 'normalized';
f.Position = [0.25 0.25 0.6 0.7];
subplot(1,9,1:8)

scatter(PLV.in.s_mean,PLV.out.s_mean)
ylabel('V1a/V4 PLV per session')
xlabel('V1n/V4 PLV per session')
hold on 
xl = xlim;
ylim(xl)
title(sprintf('PLVs per trial',params.bptype,params.lower,params.upper,num2str(params.MC)))
bool_vec = {'false','true'}
t = annotation("textbox",[0.85 0.11 0.12 0.83],'String',sprintf('Parameters:\nBPtype: %s BPbounds: %i-%iHz\nMCs: %s \nSmoothing: %s \nSmooth Window: %ims\nFilt order (sgolay): %i \nSign. Diff: %s',params.bptype,params.lower,params.upper,num2str(params.MC),params.filttype,params.framelen,params.filtord,bool_vec{value+1}))

foldername = fullfile(params.figpath,'inst_freq',params.bptype,sprintf("bpwidth_%i-%i/toi_%.1f-%.1f/PLV",params.lower,params.upper,params.toi(1),params.toi(2)))
if ~exist(foldername,'dir')
    mkdir(foldername)
end 
saveas(gcf,fullfile(foldername,'PLV_scatter.fig'))
saveas(gcf,fullfile(foldername,'PLV_scatter.jpg'))

% Plotting coherence over time
f = figure;
f.Units = "normalized";
f.Position = [0.25 0.25 0.7 0.6]

subplot(1,9,1:8)
plot(PLV.in.t_mean,'r');
xlabel('Time [ms]');
ylabel('PLV');
hold on
plot(PLV.out.t_mean,'b');
legend('V1a/V4','V1n/V4')
title(sprintf('Average PLVs over time, %s boundaries: %d, %d Hz, MC: %s',params.bptype,params.lower,params.upper,num2str(params.MC)))
t = annotation("textbox",[0.85 0.11 0.12 0.83],'String',sprintf('Parameters:\nBPtype: %s BPbounds: %i-%iHz\nMCs: %s \nSmoothing: %s \nSmooth Window: %ims\nFilt order (sgolay): %i',params.bptype,params.lower,params.upper,num2str(params.MC),params.filttype,params.framelen,params.filtord))


foldername = fullfile(params.figpath,'inst_freq',params.bptype,sprintf("bpwidth_%i-%i/toi_%.1f-%.1f/PLV",params.lower,params.upper,params.toi(1),params.toi(2)))
if ~exist(foldername,'dir')
    mkdir(foldername)
end 
saveas(gcf,fullfile(foldername,'PLV_time.fig'))
saveas(gcf,fullfile(foldername,'PLV_time.jpg'))

close all
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
