% Run Step 1 salvage on the included bad sample without touching pipeline
inp = fullfile(fileparts(mfilename('fullpath')),'..','ReVAS','samples','bad');
% Prefer the position mat if present
cand = dir(fullfile(inp,'*_position.mat'));
if isempty(cand)
    error('Could not locate a *_position.mat under %s', inp);
end
matPath = fullfile(cand(1).folder, cand(1).name);

salvage = Step1_GoodFrameMask(matPath, 'minPeak',0.55,'maxMotion',0.12,'verbose',true);

fprintf('Kept %.1f%% of samples (p90 motion=%.2f%%/fr).\n', salvage.percentKept, 100*salvage.p90MotionGood);
fprintf('Saved: %s\n', fullfile(fileparts(matPath), [""+extractBefore(cand(1).name, strlength(cand(1).name)-3) "_salvage_step1.mat"]));
