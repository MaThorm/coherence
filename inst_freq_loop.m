% Script that allows for looping through multiple bp-filter width 
clc
clear all
load('attout_dataset.mat')
load('attin_dataset.mat')
load("V4_dataset.mat")
bp_low = [30 40 50 60 70 80 90];
bp_add = [10 20 30 40 50 60 70];
toi = [-100 100];
medianfiltord = 20;
matpath = '/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files';
timebar = -1.3:0.001:5;
for i_low = 1:length(bp_low)
    for i_add = 1:length(bp_add)
        bpwidth = [bp_low(i_low) bp_low(i_low) + bp_add(i_add)];
        if bpwidth(2) > 100
            continue
        end 
        [in_trials,out_trials,V4_trials] = pre_processing_pip_trials(attin_dataset,attout_dataset,V4_dataset,bpwidth,toi)
        [grand_struct,angles,inst_freq] = pre_processing_pip_hilb(in_trials,out_trials,V4_trials,medianfiltord)
        attin_inst = grand_struct.in_medfiltHilbert;
        attout_inst = grand_struct.out_medfiltHilbert;
        insummary = grand_struct.insummary;
        outsummary = grand_struct.outsummary;
        V4_inst = grand_struct.V4_medfiltHilbert;
        V4summary = grand_struct.V4summary;
        plotting_inst(insummary,outsummary,V4summary,timebar,bpwidth)
    end 
end 

%%
function [] = plotting_inst(insummary,outsummary,V4summary,timebar,bpwidth)
% Plotting AttIn Summary
dir = sprintf('/home/mthormann@brain.uni-bremen.de/V1V4coherence/03_results_max/inst_freq_loop/%d/%d',bpwidth(1),bpwidth(2));
status = mkdir(dir);
figure('units','normalized','outerposition',[0 0 1 1],'visible','off')
sel = 1:5000;
x = timebar(sel)
y = insummary.mean(sel);
sd = insummary.std(sel);
patch([x fliplr(x)], [y-sd  fliplr(y+sd)], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
title(sprintf("AttIn Summary - bpwidth of %d & %d",bpwidth(1),bpwidth(2)))
xlabel('Time [s]')
ylabel('Freqeuncy [Hz]')
hold on 
label = {'Static', 'MS1','MS2','MS3','MS4'}
xl = xline([-0.5 0 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
plot(x,y,'r');
hold off
saveas(gcf,fullfile(dir,sprintf('AttIn_Summary_bp%d-%d.png',bpwidth(1),bpwidth(2))))

% Plotting attout Summary
sel = 1:5000;
x = timebar(sel)
figure('units','normalized','outerposition',[0 0 1 1],'visible','off')
y = outsummary.mean(sel);
sd = outsummary.std(sel);
patch([x fliplr(x)], [y-sd  fliplr(y+sd)], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
label = {'Static', 'MS1','MS2','MS3','MS4'}
xl = xline([-0.5 0 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
xlabel('Time [s]')
ylabel('Freqeuncy [Hz]')
title(sprintf("AttOut Summary - bpwidth of %d & %d",bpwidth(1),bpwidth(2)))
hold on 
plot(x,y,'r');
hold off
saveas(gcf,fullfile(dir,sprintf('AttOut_Summary_bp%d-%d.png',bpwidth(1),bpwidth(2))))


% plotting V4 summary 
figure('units','normalized','outerposition',[0 0 1 1],'visible','off')
sel = 1:5000;
x = timebar(sel)
y = V4summary.mean(sel);
sd = V4summary.std(sel);
patch([x fliplr(x)], [y-sd  fliplr(y+sd)], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
label = {'Static', 'MS1','MS2','MS3','MS4'}
xl = xline([-0.5 0 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
xlabel('Time [s]')
ylabel('Freqeuncy [Hz]')
title(sprintf("V4 Summary - bpwidth of %d & %d",bpwidth(1),bpwidth(2)))
hold on 
plot(x,y,'r');
hold off
saveas(gcf,fullfile(dir,sprintf('V4_Summary_bp%d-%d.png',bpwidth(1),bpwidth(2))))


% plotting all together 
sel = 1:5000;
x = timebar(sel);
figure('units','normalized','outerposition',[0 0 1 1],'visible','off')
plot(x,insummary.mean(sel),'r');
hold on 
plot(x,outsummary.mean(sel),'b');
plot(x,V4summary.mean(sel));
title(sprintf("Summary of all conditions - bpwidth of %d & %d",bpwidth(1),bpwidth(2)))
legend('AttIn','AttOut','V4','AutoUpdate','off')
label = {'Static', 'MS1','MS2','MS3','MS4'};
xl = xline([-0.5 0 1 2 3],'--',label,'color',[0.7 0.7 0.7]);
hold off
xlabel('Time [s]')
ylabel('Frequency [Hz]')
saveas(gcf,fullfile(dir,sprintf('Grand_Summary_bp%d-%d.png',bpwidth(1),bpwidth(2))))

end 