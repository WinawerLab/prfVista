function results = prfVistasoft(stimfiles, datafiles, stimradius, varargin)
% Create a temporary, hidden vistasoft session and solve pRF model
%
% Inputs 
%   stimfiles:  path to one or more files with stimulus description, char or
%                   cell of chars
%   datafiles:  path to one or more files with BOLD time series, char or
%                   cell of chars
%   stimradius: radius in degrees, scalar 
%   varargin:   input pairs for <model>, <wsearch>, <detrend>, <keepallPoints>, <numberStimulusGridPoints>
%                   These are optional inputs to rmMain. 
%   TODO: Update this function to simply pass all paired arguments to rmMain. No need to parse them here.
%
% Outputs
%   results:    mrVista data structure with model results and model parameters
%
% Examples:
%   results = prfVistasoft(fullfile(prfRootPath, 'local', 'vista2test_Stim.nii.gz'), fullfile(prfRootPath, 'local', 'vista2test_Data.nii.gz'), 10);
%   
%   See README.md for several more detailed examples


% First clear the workspace just in case
mrvCleanWorkspace

% Read the inputs
varargin = mrvParamFormat(varargin);
p = inputParser;
p.addRequired('stimfiles'                       , @(x) or(ischar(x), iscell(x)));
p.addRequired('datafiles'                       , @(x) or(ischar(x), iscell(x)));
p.addRequired('stimradius'                      , @isnumeric);
p.addParameter('model'        , 'one gaussian'  , @ischar);
p.addParameter('wsearch'      , 'coarse to fine', @ischar);
p.addParameter('detrend'      , 1               , @isnumeric);
p.addParameter('keepAllPoints', true            , @islogical);
p.addParameter('numberStimulusGridPoints', 50   , @isnumeric);
p.addParameter('tr', 1                          , @isnumeric);
p.addParameter('hrfparams', 'two gammas (SPM style)', @ischar);
p.addParameter('decimate'      , 2               , @isnumeric);
p.addParameter('calcPC'        , true            , @islogical);

p.parse(stimfiles, datafiles, stimradius, varargin{:});

% Assign it
model          = p.Results.model;
wSearch        = p.Results.wsearch;
detrend        = p.Results.detrend;
keepAllPoints  = p.Results.keepAllPoints;
numGridPoints  = p.Results.numberStimulusGridPoints;
tr             = p.Results.tr;
hrfparams      = p.Results.hrfparams;
decimatefactor = p.Results.decimate;
calcPC         = p.Results.calcPC;

% How many scans?
if iscell(stimfiles)
    numscans = length(stimfiles);
else
    numscans = 1;
    stimfiles = {stimfiles};
    datafiles = {datafiles};
end



%% Set up files and directories
homedir = fullfile(tempdir, 'vistaPRF');
mkdir(homedir)
cd(homedir);

if exist(fullfile(homedir,'Raw'),'dir');warning('RAW DIR EXISTS');end
mkdir(fullfile(homedir, 'Raw'));
mkdir(fullfile(homedir, 'Stimuli'));

%% move the data into the new directories
for ii = 1:numscans
    fprintf('[pmVistasoft] Stim path %d: %s\n',ii, stimfiles{ii})
    fprintf('[pmVistasoft] Data path %d: %s\n',ii, datafiles{ii})
    copyfile(datafiles{ii}, fullfile(homedir, 'Raw'));
    copyfile(stimfiles{ii}, fullfile(homedir, 'Stimuli'));
end
fprintf('\n[pmVistasoft] This is stimradius: %i\n',stimradius)

%% convert stim files to .mat format that vistasoft can read
stimfileMat = cell(1, numscans);
for ii = 1:numscans
    [~, ~, e] = fileparts(stimfiles{ii});
    switch e
        case {'.nii' '.gz'}
            ni = niftiRead(stimfiles{ii});
            
            images = squeeze(ni.data);
            pixdim = niftiGet(ni, 'pixdim');
            tr     = pixdim(end);
        case '.mat'
            s = load(stimfiles{ii});
            if ~isfield(s, 'stim'), error('Cannot find stim field in .mat stim file'); end
            images = s.stim;
    end            
    
    sprintf('\n\n USING TR:%2.2f \n\n',tr)
    stimulus.seq = 1:size(images,3);
    stimulus.seqtiming = (stimulus.seq-1) * tr;
    
    stimfileMat{ii} = fullfile('.', 'Stimuli', sprintf('images_and_params_scan%d', ii));
    save(stimfileMat{ii}, 'images', 'stimulus');
end

%% create a pseudo inplane underlay, required by vistasoft, by averaging the
%   time series for each voxel
fmri        = niftiRead(datafiles{1});
ippath      = fullfile('.', 'Raw', 'inplane.nii.gz');
ip          = fmri; 
ip.data     = mean(fmri.data, length(size(fmri.data)));
ip.dim      = size(ip.data);
ip.ndim     = numel(ip.dim);
niftiWrite(ip, ippath);

% A = niftiRead(ippath)


%% Set up the vistasoft session
params = mrInitDefaultParams;

params.sessionDir   = homedir;
params.vAnatomy     = [];
params.inplane      = ippath;

for ii = 1:numscans
    [~, f, e] = fileparts(datafiles{ii});
    params.functionals{ii}  = fullfile('.','Raw', sprintf('%s%s', f,e));
end
% Run it:
ok = mrInit(params); if ~ok, error('mrInit failed'); end
dir(fullfile('.', filesep,'Raw'))


dir(fullfile('.', filesep,'Raw'))



%% Check it
%{
vw = initHiddenInplane();
sz =  viewGet(vw,'anatsize');
[a, b, c] = ind2sub(sz, 1:prod(sz));
coords = [a; b; c];
vw = newROI(vw,'all',true,'w',coords, 'all voxels');

% plot raw signal
newGraphWin();
plotMeanTSeries(vw, viewGet(vw, 'current scan'), [], true);

% plot percent signal modulation
newGraphWin();
plotMeanTSeries(vw, viewGet(vw, 'current scan'), [], false);
%}

%% Set up prf model

vw = initHiddenInplane();
% edit GLU: dataTYPES is not found, but it was stablished as global in mrInit()
% Load mrSESSION in here to see if this solves it
load(fullfile(homedir,'mrSESSION.mat'), 'dataTYPES')
disp(dataTYPES)
dataTYPES.scanParams
% Set default retinotopy stimulus model parameters
sParams = rmCreateStim(vw);

for ii = 1:length(sParams)
    sParams(ii).stimType   = 'StimFromScan'; % This means the stimulus images will
    % be read from a file.
    sParams(ii).stimSize   = stimradius;     % stimulus radius (deg visual angle)
    sParams(ii).nDCT       = detrend;        % detrending frequeny maximum (cycles
    % per scan): 1 means 3 detrending
    % terms, DC (0 cps), 0.5, and 1 cps
    sParams(ii).imFile     = stimfileMat{ii};    % file containing stimulus images
    sParams(ii).paramsFile = stimfileMat{ii};    % file containing stimulus parameters
    % 'thresholdedBinary',  whenreading in images, treat any pixel value
    %                       different from background as a 1, else 0
    sParams(ii).imFilter   = 'none';
    % this is a string. see hrfGet for possible values.
    sParams(ii).hrfType    = hrfparams;
    % pre-scan duration will be stored in frames for the rm, but was stored in
    % seconds in the stimulus file
    sParams(ii).prescanDuration = 0;
end

dataTYPES = dtSet(dataTYPES, 'rm stim params', sParams);

saveSession();
% Check it
%   vw = rmLoadParameters(vw);
%   [~, M] = rmStimulusMatrix(viewGet(vw, 'rmparams'), [], [], 1, false);


%% Solve prf Models
vw = initHiddenInplane();

vw = rmMain(vw, [], wSearch, ...
            'model', {model}, ...
            'matFileName', 'tmpResults', ...
            'keepAllPoints', keepAllPoints, ...
            'numberStimulusGridPoints', numGridPoints, ...
            'decimate', decimatefactor, ...
            'calcPC', calcPC);

% Load the results        
d = dir(fullfile(dataDir(vw), sprintf('%s*', 'tmpResults')));
[~,newestIndex] = max([d.datenum]);
results = load(fullfile(dataDir(vw), d(newestIndex).name));

%% Delete all global variables created by mrVista
mrvCleanWorkspace

% Delete the temp folder with all the tmp files, we only want the results
cd(homedir)
cd('../')
rmdir(homedir, 's')

%% TODO: Convert results to mgz or some other standardized format?



end
