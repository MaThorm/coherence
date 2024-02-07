function [trials] = do_trialselection(path,filename,CoI,StimNo,bpwidth,toilim,BlCor,IncOut,offset)
%DO_TRIALSELECTION 
%   path = path of the file
%   filename = name of the file
%   CoI = Channels of Interest 
%   StimNo = Included Stimulus Numbers 
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

% selecting certain time window
if exist('toilim') 
    cfg = [];
    cfg.toilim = toilim;
end 

% shift 0 by offset
cfg = [];
cfg.offset = -round(offset*data.hdr.Fs);
trials = ft_redefinetrial(cfg, trials);


%{  
remove average
if BlCor == true
    preprocCfg = [];
    preprocCfg.demean = 'yes';
    preprocCfg.baselinewindow = [-1.3 -0.5]
    trials = ft_preprocessing(preprocCfg, trials);
end 
%}
%
% BP filtering
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = bpwidth;
trials = ft_preprocessing(cfg,trials);
%
end 