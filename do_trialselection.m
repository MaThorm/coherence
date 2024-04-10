function [trials] = do_trialselection(path,filename,CoI,StimNo,bpfilt,varargin)
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

% Input parsing

p = inputParser;

dbpwidth = [1 150];
dtoilim = [-100 100];
dblcor = false;
dincout = [0];
doffset = 1.3;

addRequired(p,'path', @(x) ischar(x));
addRequired(p,'filename', @(x) ischar(x));
addRequired(p,'CoI',@(x) isnumeric(x));
addRequired(p,'StimNo',@(x) isnumeric(x));
addRequired(p,'bpfilt',@(x) islogical(x));
addParameter(p,'bpwidth',dbpwidth,@(x) isnumeric(x));
addParameter(p,'toilim',dtoilim,@(x) isnumeric(x));
addParameter(p,'blcor',dblcor,@(x) islogical(x)); 
addParameter(p,'IncOut',dincout,@(x) isnumeric(x));
addParameter(p,'offset',doffset, @(x) isnumeric(x));
parse(p,path,filename,CoI,StimNo,bpfilt,varargin{:});

path = p.Results.path;
filename = p.Results.filename;
CoI = p.Results.CoI;
StimNo = p.Results.StimNo;
bpfilt = p.Results.bpfilt;
toilim = p.Results.toilim;
blcor = p.Results.blcor;
IncOut = p.Results.IncOut;
offset = p.Results.offset;
 

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
cfg = [];
cfg.toilim = toilim;
trials = ft_redefinetrial(cfg, trials);

% shift 0 by offset
cfg = [];
cfg.offset = -round(offset*data.hdr.Fs);
trials = ft_redefinetrial(cfg, trials);

if blcor == true
    preprocCfg = [];
    preprocCfg.demean = 'yes';
    preprocCfg.baselinewindow = [-1.3 -0.5]
    trials = ft_preprocessing(preprocCfg, trials);
end 

% BP filtering
if bpfilt == true
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = p.Results.bpwidth;
    cfg.bpfilttype = 'fir'; 
    trials = ft_preprocessing(cfg,trials);
end 
%
end 