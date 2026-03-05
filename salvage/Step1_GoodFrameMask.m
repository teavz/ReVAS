function salvage = Step1_GoodFrameMask(inputHint, varargin)
% Step1_GoodFrameMask  Compute a conservative good-frame mask from existing outputs.
%
% Usage
%   salvage = Step1_GoodFrameMask(inputHint, 'minPeak',0.55,'maxMotion',0.12,
%                                 'minSegFrames',10,'maxInterpGap',8)
%
% Inputs
%   inputHint : path to any related .mat output (e.g., *_strip.mat,
%               *_position.mat, *_reference.mat) or a folder containing them.
%   Name-Value options:
%     - minPeak      : minimum acceptable correlation peak (default 0.55)
%     - maxMotion    : maximum inter-strip motion as a fraction of frame
%                      height per frame ("%/fr" in plots, default 0.12)
%     - verbose      : true to produce a QC figure (default true)
%     - outdir       : directory to write outputs (default: alongside mats)
%
% Output struct fields (also saved to <base>_salvage_step1.mat):
%   .goodMask          logical Nx1 mask on strip samples (timeSec-aligned)
%   .timeSec           Nx1 timestamps (s)
%   .peakValueArray    Nx1 peaks if available, else []
%   .deltaPos          Nx1 motion fraction per frame
%   .thresholds        struct with minPeak, maxMotion
%   .percentKept       scalar percent of good samples
%   .notes             cellstr notes about any fallbacks
%   .sourceFiles       struct with paths used
%
% This function does NOT modify any pipeline file. It reads existing mats,
% computes the mask, and writes a compact QC report + PNG.

ip = inputParser;
ip.addParameter('minPeak', 0.55, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
ip.addParameter('maxMotion', 0.12, @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=1);
ip.addParameter('verbose', true, @(x)islogical(x)||ismember(x,[0 1]));
ip.addParameter('outdir','', @(x)ischar(x)||isstring(x));
ip.parse(varargin{:});
opts = ip.Results;

notes = {};

% Discover related files
if isfolder(inputHint)
    cand = dir(fullfile(inputHint,'*.mat'));
    files = {cand.name};
    root = inputHint;
else
    root = fileparts(inputHint);
    files = {dir(fullfile(root,'*.mat')).name};
end

% Helper to find first file matching pattern
findfirst = @(pat) string(fullfile(root, files{find(~cellfun(@isempty,regexpi(files,pat)),1,'first')})) ;

stripMat = ""; posMat = ""; refMat = "";
if any(~cellfun(@isempty, regexpi(files,'_strip(\.\w+)?\.mat$'))) 
    stripMat = findfirst('_strip(\.\w+)?\.mat$');
end
if posMat == "" && any(~cellfun(@isempty, regexpi(files,'_position(\.\w+)?\.mat$'))) 
    posMat = findfirst('_position(\.\w+)?\.mat$');
end
if any(~cellfun(@isempty, regexpi(files,'_reference\.mat$'))) 
    refMat = findfirst('_reference\.mat$');
end

if stripMat == "" && isfile(inputHint)
    % Use the provided mat directly if it looks relevant
    [~,nm,ext] = fileparts(inputHint);
    if endsWith(ext,'.mat') && contains(nm, {'strip','position','reference'})
        switch true
            case contains(nm,'strip')
                stripMat = string(inputHint);
            case contains(nm,'position')
                posMat = string(inputHint);
            case contains(nm,'reference')
                refMat = string(inputHint);
        end
    end
end

% Load what we can
peakValueArray = [];
rawPosition = [];
timeSec = [];
params = struct; 
source = struct('strip', stripMat, 'position', posMat, 'reference', refMat);

if stripMat ~= "" && isfile(stripMat)
    s = load(stripMat);
    source.strip = stripMat;
    if isfield(s,'peakValueArray'), peakValueArray = s.peakValueArray; end
    if isfield(s,'rawPosition'), rawPosition = s.rawPosition; end
    if isfield(s,'timeSec'), timeSec = s.timeSec; end
    if isfield(s,'params'), params = s.params; end
end

% If missing, try position mat
if isempty(rawPosition) && posMat ~= "" && isfile(posMat)
    s = load(posMat);
    source.position = posMat;
    if isfield(s,'rawPosition')
        rawPosition = s.rawPosition;
    elseif isfield(s,'position')
        % Fall back to masked position if raw not saved
        rawPosition = s.position; %#ok<NASGU>
        notes{end+1} = 'rawPosition missing; used position as proxy'; %#ok<AGROW>
    end
    if isempty(timeSec) && isfield(s,'timeSec'), timeSec = s.timeSec; end
    if isempty(fieldnames(params)) && isfield(s,'params'), params = s.params; end
end

% If peaks still missing, try reference mat (ReReference/MakeReference)
if isempty(peakValueArray) && refMat ~= "" && isfile(refMat)
    s = load(refMat);
    source.reference = refMat;
    if isfield(s,'peakValues')
        peakValueArray = s.peakValues; % may differ in length; handle below
        notes{end+1} = 'Using peakValues from reference mat (length may differ).'; %#ok<AGROW>
    end
    if isempty(fieldnames(params)) && isfield(s,'params'), params = s.params; end
end

% Validate essentials
assert(~isempty(timeSec) && ~isempty(rawPosition), ...
    'Step1_GoodFrameMask: could not locate timeSec/rawPosition in provided files.');

% Recover height and stripsPerFrame
if isfield(params,'referenceFrame')
    [height, ~] = size(params.referenceFrame);
else
    % Conservative fallback: infer from rowNumbers spacing vs. time step
    height = 512; % common default; only affects labeling, not gating limits when using fraction
    notes{end+1} = 'referenceFrame missing; assumed height=512 for motion normalization.'; %#ok<AGROW>
end

% Compute deltaPos (motion fraction per frame)
N = size(rawPosition,1);
stripsPerFrame = round(params.samplingRate / params.frameRate);
if ~isfield(params,'samplingRate') || ~isfield(params,'frameRate') || ~isnan(stripsPerFrame)==0
    % Fallback: estimate stripsPerFrame from rowNumbers if available
    if isfield(params,'rowNumbers') && numel(params.rowNumbers)>1
        stripsPerFrame = numel(params.rowNumbers);
    else
        stripsPerFrame = 20; % typical value
        notes{end+1} = 'samplingRate/frameRate missing; assumed stripsPerFrame=20.'; %#ok<AGROW>
    end
end

% deltaPos as in StripAnalysis: normalize by frame height and per-strip dt
deltaPos = nan(N,1); deltaPos(1) = 0;
if N>1
    dtPerScanLine = (timeSec(2)-timeSec(1)) / (params.rowNumbers(2)-params.rowNumbers(1));
    for i=2:N
        dr = (rawPosition(i,:)-rawPosition(i-1,:))/height;
        dt = (timeSec(i)-timeSec(i-1)) / (dtPerScanLine * height / stripsPerFrame);
        deltaPos(i) = sqrt(sum(dr.^2,2)) * stripsPerFrame / max(dt,eps);
    end
end
% The original code multiplies by stripsPerFrame; above already scales. Keep behavior consistent:
deltaPos = sqrt([0; sum(diff(rawPosition,1,1).^2,2)] )/height; % per-strip normalized
% Convert to per-frame equivalent like MakeReference: multiply by stripsPerFrame
if stripsPerFrame>0
    deltaPos = deltaPos * stripsPerFrame;
end

% Align peak array length, if necessary
if ~isempty(peakValueArray) && numel(peakValueArray)~=N
    % Try to interpolate or truncate
    if numel(peakValueArray) > 2 && N > 2
        x1 = linspace(0,1,numel(peakValueArray));
        x2 = linspace(0,1,N);
        peakValueArray = interp1(x1, peakValueArray, x2, 'linear','extrap');
        notes{end+1} = 'Resampled peak array to match timeSec length.'; %#ok<AGROW>
    else
        peakValueArray = peakValueArray(:);
        peakValueArray = peakValueArray(1:min(end,N));
        if numel(peakValueArray)<N
            peakValueArray(end+1:N,1) = NaN; %#ok<AGROW>
        end
        notes{end+1} = 'Truncated/NaN-padded peak array to match length.'; %#ok<AGROW>
    end
end

% Compute good mask
if isempty(peakValueArray)
    goodMask = deltaPos <= opts.maxMotion;
    notes{end+1} = 'Peaks unavailable; using motion-only gating.'; %#ok<AGROW>
else
    goodMask = (peakValueArray >= opts.minPeak) & (deltaPos <= opts.maxMotion);
end

% Summaries
percentKept = 100 * nnz(goodMask) / numel(goodMask);
medPeakGood = []; 
p90MotionGood = prctile(deltaPos(goodMask),90);
if ~isempty(peakValueArray)
    medPeakGood = median(peakValueArray(goodMask),'omitnan');
end

salvage = struct();
salvage.goodMask = logical(goodMask(:));
salvage.timeSec = timeSec(:);
salvage.peakValueArray = peakValueArray(:);
salvage.deltaPos = deltaPos(:);
salvage.thresholds = struct('minPeak',opts.minPeak,'maxMotion',opts.maxMotion);
salvage.percentKept = percentKept;
salvage.medianPeakGood = medPeakGood;
salvage.p90MotionGood = p90MotionGood;
salvage.notes = notes(:);
salvage.sourceFiles = source;

% Write outputs
if isempty(opts.outdir)
    outdir = root;
else
    outdir = char(opts.outdir);
end
base = guess_base_name(source);
outMat = fullfile(outdir, base + "_salvage_step1.mat");
save(outMat, 'salvage');

% QC Figure
if opts.verbose
    fh = figure('visible','off');
    tiledlayout(fh,2,1,'Padding','compact','TileSpacing','compact');

    % Peak scatter
    nexttile;
    if ~isempty(peakValueArray)
        scatter(timeSec, peakValueArray, 8, 'filled'); hold on;
        yline(opts.minPeak,'-','linewidth',1.5);
        ylabel('peak value'); xlim([min(timeSec) max(timeSec)]);
        ylim([0 1]); grid on; title('Peaks and threshold');
    else
        text(0.5,0.5,'Peaks unavailable','HorizontalAlignment','center'); axis off;
    end

    % Motion scatter
    nexttile;
    scatter(timeSec, 100*deltaPos, 8, 'filled'); hold on;
    yline(100*opts.maxMotion,'-','linewidth',1.5);
    ylabel('motion (%/fr)'); xlabel('time (s)');
    xlim([min(timeSec) max(timeSec)]);
    ylim([0 100*max(opts.maxMotion, prctile(deltaPos,99)+0.02)]);
    grid on; title('Motion and threshold');

    outPng = fullfile(outdir, base + "_salvage_step1_qc.png");
    exportgraphics(fh, outPng, 'Resolution', 150);
    close(fh);
    salvage.qcFigure = outPng;
end

end

function base = guess_base_name(source)
% Create a reasonable base from whichever source exists
f = '';
if isfield(source,'position') && strlength(source.position)>0
    f = char(source.position);
elseif isfield(source,'strip') && strlength(source.strip)>0
    f = char(source.strip);
elseif isfield(source,'reference') && strlength(source.reference)>0
    f = char(source.reference);
else
    base = "salvage"; return;
end
[~,nm,~] = fileparts(f);
base = string(regexprep(nm, '(_strip|_reference|_position(_deg)?|_filtered(_sacsdrifts)?)$', ''));
end
