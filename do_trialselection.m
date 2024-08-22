function [trials] = do_trialselection(path,filename,CoI,StimNo,bpfilt,bpwidth,order,BlCor,IncOut,offset)
%DO_TRIALSELECTION 
%   path = path of the file
%   filename = name of the file
%   CoI = Channels of Interest 
%   StimNo = Included Stimulus Numbers 
%   bpwidth = widht of bp filter in Hz
%   toilim = times of interest in seconds 
%   BlCor = Baseline correction (default = false)
%   IncOut = Included Outcomes (default = [0])
%   offset = value in s by which x-axis gets reoriented (default = 1.3)

if ~exist('BlCor'), BlCor = false;end
if ~exist('IncOut'), IncOut = [0];end
if ~exist('offset'), offset = 1.3;end

filename = fullfile(path,filename);
vblTrigger = 15;
trialdef.alignTrigger = 5;
trialdef.startTrigger = 5;
trialdef.endTrigger = -11;
trialdef.tPreStartTrigger = 0;
trialdef.tPostEndTrigger = 0;

data = dh.ft_read_cont(filename, CoI);
header = dh.ft_read_header(filename, CoI);
events = dh.ft_read_event(filename, header);
eventsWithoutVBL = events(~arrayfun(@(x) strcmp('trigger', x.type) && abs(x.value) == vblTrigger, events));


% trial definition
cfg = [];
cfg.event = eventsWithoutVBL;
cfg.hdr = header;
cfg.trialdef = trialdef;
[trl, event] = dh.trialfun_general(cfg);
cfg = [];
cfg.trl = trl;
trials = ft_redefinetrial(cfg, data);

% trial selection
selectCfg = [];
selectCfg.trials = ismember(trials.trialinfo(:,3), IncOut) & ...
    ismember(trials.trialinfo(:,2), StimNo);
trials = ft_selectdata(selectCfg, trials);

% % BP filtering
if bpfilt == true
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = bpwidth;
    cfg.padding = 4;
%     cfg.padtype = 'mirror'
    cfg.bpfilttype = 'fir'; 
    cfg.bpfiltord = order;
    trials = ft_preprocessing(cfg,trials);
end 

%new, manual bp-filtering
% fs = 1000;
% Fbp = bpwidth;
% type = 'fir';
% dir = 'twopass';
% if bpfilt == true
%     for i_t = 1:length(trials.trial)
%         for i_c = 1:size(trials.trial{1,i_t},1)
%             cur_trial = trials.trial{1,i_t}(i_c,:);
%             filt_trial = ft_preproc_bandpassfilter(cur_trial,fs,Fbp,order,type,dir);
%             trials.trial{1,i_t}(i_c,:) = filt_trial;
%         end 
%     end 
% end 

% % selecting certain time window
% if exist('toilim') 
%     cfg = [];
%     cfg.toilim = toilim;
%     cfg.minlength = 0.5; % minimum length of trials in seconds
%     trials = ft_redefinetrial(cfg, trials);
% end 

% shift 0 by offset
cfg = [];
cfg.offset = -round(offset*data.hdr.Fs);
trials = ft_redefinetrial(cfg, trials);


if BlCor == true
    preprocCfg = [];
    preprocCfg.demean = 'yes';
    preprocCfg.baselinewindow = [-1.3 -0.5]
    trials = ft_preprocessing(preprocCfg, trials);
end 


%
end 