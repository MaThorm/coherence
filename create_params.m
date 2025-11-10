%% Parameter Settings

% Trial preprocessing parameters
params.bptype = "bpfilt"; % Should be either 'SSD' or 'bpfilt' 
params.bpwidth = [48 82];
params.toi = [1 3]; % Time region of interest [1 3] is MC 2&3
params.order = 1000;

% SSD parameters
params.fs = 1000; % Sampling frequency of the signal
params.th = 0.01; % residual variance threshhold   
params.lower = params.bpwidth(1); % Changed it now so that they are the same as the bp width filters
params.upper = params.bpwidth(2);


% Hilbert parameters
params.filttype = 'medfilt'; %either 'medfilt' or 'sgolay'
params.framelen = 31;
params.filtord = 1;

% snipping parameters
params.snipping = false;
params.peak_time = 3;
params.peak_prcntl = 20;
params.min_time = 0.00001;

% averaging parameters
params.avg_type = 'avg'; %either 'weight' or 'avg' 

% color plotting paramers
params.blue = [0 0.4470 0.7410];
params.red = [0.8500 0.3250 0.0980];
params.orange = [0.9290 0.6940 0.1250];
params.pl_l = [0.25 0.25 0.25 0.25];
params.pl_sq = [0.25 0.25 0.21 0.3];
params.pl_trip = [0.25 0.25 0.5 0.5];
params.pl_half = [0.25 0.25 0.10 0.3];

% Configuration that follows
if params.bptype == "bpfilt"
    params.bpfilt = true;
elseif params.bptype == "SSD"
    params.bpfilt = false;
end 

if isequal(params.toi, [1 3])
    params.MC = [2 3];
elseif isequal(params.toi, [1 2])
    params.MC = [2];
elseif isequal(params.toi, [2 3])
    params.MC = [3];
end

params.matpath = '/data/projects/V1V4coherence/02_analysis_max/git_repos/mat_files';
if params.snipping == false;
    params.figpath = '/data/projects/V1V4coherence/03_results_max/results_summary';
else 
    params.figpath = '/data/projects/V1V4coherence/03_results_max/results_summary/snippeted';
end 

params.inst_figpath = fullfile(params.figpath,params.bptype,sprintf("bpwidth_%i-%i/toi_%.1f-%.1f/%s_%i/",params.lower,params.upper,params.toi(1),params.toi(2),params.filttype,params.framelen),params.avg_type)


if strcmp(params.filttype, 'medfilt') && params.filtord > 1
    error('In case of medilt, filtord should be 1 for consistency')
end 

save("/data/projects/V1V4coherence/02_analysis_max/git_repos/params.mat","params")
