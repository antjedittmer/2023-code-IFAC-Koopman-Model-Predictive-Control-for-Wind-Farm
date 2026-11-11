function files = plotAdaptiveAICFlowFields(tSel, dirFig, rerun)
%PLOTADAPTIVEAICFLOWFIELDS Wind fields of the adaptive AIC closed loop (thesis Section 5.2)
%   plotAdaptiveAICFlowFields() plots the longitudinal wind fields of the
%   adaptive AIC closed-loop run of mainGeneratePlotThesisAdaptiveAIC
%   (sol_update = 1, default update settings, UrFromPower = 1) with
%   WFSim_flow_tiled: two figures with two snapshots each, flow from left
%   to right, deviation u - V_inf from the inflow speed of each snapshot
%   with one colour scale for all tiles (white: free stream), titles with
%   time, P_ref and C_T' of both turbines.
%
%   plotAdaptiveAICFlowFields(tSel, dirFig, rerun)
%   - tSel: times in s, one row per figure (default [240 2100; 2520 4990]:
%     start, after the V_inf ramp from 7 to 9 m/s at 2000 s, after the
%     P_ref step at 2500 s, steady state); multiples of 5 s
%   - dirFig: output folder (default figPaper in the main folder)
%   - rerun: true repeats the closed-loop run (default: only if the
%     snapshot file does not exist)
%   files: names of the saved figures (png, eps)
%
%   WFSim_demo saves only the time series of the run. The wind fields are
%   therefore taken from a repetition of the run and saved every 5 s in
%   AdaptiveAIC_flowSnapshots.mat in the folder of the run. The saved
%   results of the run (figures, _plotData.mat) are restored afterwards,
%   so the thesis data are not changed.

if nargin < 1 || isempty(tSel)
    tSel = [240 2100; 2520 4990];
end
mainDir = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
if nargin < 2 || isempty(dirFig)
    dirFig = fullfile(mainDir, 'figPaper');
end
if nargin < 3 || isempty(rerun)
    rerun = false;
end
if ~isfolder(dirFig)
    mkdir(dirFig);
end

dirCL = fullfile(mainDir, 'DataT2OLWFSim', 'Vinf7dot0To9dot0_CL_K06_P0');
runBase = 'Stairs20R1dot0e_neg06_AICnoV1_UrP_Update1';                    % adaptive run, UrFromPower = 1
snapFile = fullfile(dirCL, 'AdaptiveAIC_flowSnapshots.mat');
if rerun || ~isfile(snapFile)
    runSnapshots(mainDir, dirCL, runBase, snapFile);
end
S = load(snapFile);
tmp = load(fullfile(dirCL, [runBase, '_plotData.mat']), 'pd');
pd = tmp.pd;
tSnap = round([S.snap.time]);

% One colour scale for all figures: deviation u - V_inf from the inflow
% speed of each snapshot (7 or 9 m/s), colormap hot (WFSim_flow_tiled)
iAll = arrayfun(@(tt) find(tSnap == tt, 1), tSel(:));
if numel(iAll) < numel(tSel)
    error('plotAdaptiveAICFlowFields:time', 'Snapshots only every 5 s (times %s).', mat2str(tSel));
end
dU = arrayfun(@(s) s.u(5:end-4,5:end-4) - s.u(1, ceil(end/2)), S.snap(iAll), 'UniformOutput', false);
opts = struct('relative', true, ...
    'cLim', [min(cellfun(@(d) min(d, [], 'all'), dU)), max(cellfun(@(d) max(d, [], 'all'), dU))]);

files = cell(size(tSel,1), 1);
for f = 1:size(tSel,1)
    solP = cell(1, size(tSel,2));
    opts.titles = solP;
    for k = 1:size(tSel,2)
        s = S.snap(tSnap == tSel(f,k));
        solP{k} = struct('u', s.u, 'turbInput', s.turbInput);
        opts.titles{k} = sprintf('t = %d s: P_{ref} = %.2f MW, C_T'' = %.2f / %.2f', ...
            tSel(f,k), pd.Pref(tSel(f,k))/1e6, s.CT_prime);
    end
    nF = figure;
    pos = get(nF, 'Position');
    set(nF, 'Position', [pos(1:3), pos(4)*0.5*size(tSel,2)]);            % default MATLAB width, height of the tiles
    WFSim_flow_tiled(S.Wp, solP, [], [], [], nF, opts);
    files{f} = fullfile(dirFig, ['adaptiveAIC', sprintf('_t%d', tSel(f,:)), '_FlowTiled']);
    set(nF, 'Name', files{f});
    print(nF, files{f}, '-dpng');
    print(nF, files{f}, '-depsc');
end
end

function runSnapshots(mainDir, dirCL, runBase, snapFile)
% Repeat the adaptive closed loop (settings as mainGeneratePlotThesisAdaptiveAIC)
% and save the wind field every 5 s; the saved results are restored afterwards
oldDir = cd(mainDir);
restoreDir = onCleanup(@() cd(oldDir));
addpath('WFSimCode'); addpath('KoopmanIODMD'); addpath(fullfile('WFSimCode','bin','Kmpc'));

% Back up the results of the run (WFSim_demo does not repeat a run with
% existing _plotData.mat and overwrites the figures)
bakDir = tempname;
mkdir(bakDir);
bak = dir(fullfile(dirCL, [runBase, '*']));
for i = 1:numel(bak)
    copyfile(fullfile(dirCL, bak(i).name), bakDir);
end
restoreRun = onCleanup(@() restoreFiles(bakDir, dirCL, bak)); % also after an error
delete(fullfile(dirCL, [runBase, '_plotData.mat']));

warnState = warning('off', 'MATLAB:nearlySingularMatrix');
restoreWarn = onCleanup(@() warning(warnState));
R = 1e-6; refstairs = 20; measured = 0; KoopmanStates = 6; PolyLiftingFunction = 0;
controller = 2; ControlSetStr = 'sowfa_2turb_aic_5000'; vinfStr = ''; kTimestep = 2;
OptYaw = 0; KoopAIC = 1; Rvec = [1,1,0.001]; RU = 100; KoopFeat = 0; AICnoV1 = 1;
Vinf = min(7 + 0.1*max((1:5001) - 2000, 0), 9);
clear qLPVMPCcontroller_dU_quadprog % reset the persistent buffers of the controller
rng(0); % excitation signal of the adaptive scheme
[sol_array, JR, JQ] = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller, ...
    ControlSetStr,Vinf,vinfStr,kTimestep,1,OptYaw,KoopAIC,Rvec,RU,KoopFeat,AICnoV1,[],1);
close all;

kSave = [1, 5:5:numel(sol_array)];
snap = struct('k', num2cell(kSave));
for i = 1:numel(kSave)
    s = sol_array(kSave(i));
    snap(i).time = s.time;
    snap(i).u = s.u;
    snap(i).turbInput = s.turbInput;
    snap(i).power = s.turbine.power;
    snap(i).CT_prime = s.turbine.CT_prime;
end
% Mesh and turbines as in WFSim_demo
Wp = layoutSet_sowfa_2turb_alm_turbl();
Wp.site.u_Inf = Vinf(1);
Wp = InitWFSim(Wp, solverSet_default(Wp), 0);
save(snapFile, 'snap', 'Wp', 'JR', 'JQ');
fprintf('Adaptive AIC run repeated: TE = %.1f kW, 1000*AA = %.1f\n', JQ*1e3, JR*1e3);
end

function restoreFiles(bakDir, dirCL, bak)
% Restore the backed up results of the run
for i = 1:numel(bak)
    copyfile(fullfile(bakDir, bak(i).name), dirCL);
end
rmdir(bakDir, 's');
end
