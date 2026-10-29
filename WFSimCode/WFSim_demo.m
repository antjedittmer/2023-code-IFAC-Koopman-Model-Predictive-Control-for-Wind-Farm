function [sol_array,JR,JQ,fileName] = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,...
    controller,ControlSetStr,Vinf, vinfStr,kTimestep,sol_update,OptYaw,KoopAIC,Rvec,RU,KoopFeat,AICnoV1,ad)
% WFSim_demo is the WIND FARM SIMULATOR (WFSIM) by S. Boersma and B. Doekemeijer
% from https://github.com/TUDelft-DataDrivenControl/WFSim
%
% A. Dittmer, B. Sharan: Reformulated as a function to enable running
% as a test suite with different input combinations.
%
% Inputs
% - R: Weights on control input changed (J = sum(e'Qe + dU'Rdu)
% - refstairs: Set reference to stairs
% - measured (for Koopman model): Use measured values as feedback
% - KoopmanStates (for Koopman model):  Number of Koopman states
% - PolyLiftingFunction(for Koopman model): Higher order lifting function
% - controller: integer switch: Open-loop: 0,  Wfsim NMPC: 1, KMPC: 2
% - KoopFeat (KMPC with KoopAIC = 0): 1 uses the input-only Koopman power
%   model with input features (KoopFeatModel.mat, exportPowerFeatureModel)
%   instead of the model with the measured wind speed V1. Default: 0
% - AICnoV1 (KMPC with KoopAIC = 1): 1 uses the effective wind speed model
%   with the inputs CT1, CT2 only (identifyAICnoV1, no V1), updated online
%   with sol_update = 1 (adaptive AIC). Default: 0
% - ad (AICnoV1): settings of the online update (koopmanAdaptStep), e.g.
%   struct('nu',0.1,'h',10); fields not given as in MPCinitv (h = 5,
%   nu = 0.15 m/s, resetState, trigger on Vr1). Default []: MPCinitv.
%   Results are saved with the suffix _nu<100*nu>_h<h>
%
% Closed loop with KMPC (controller = 2, no library update): the plot data of
% each run are saved as <figure name>_plotData.mat in the output folder. If
% this file exists, the run is not repeated and the figures are recreated
% from it (plotWFSimTimeseriesTex); delete the file to run it again.

rng('default');
if ~nargin
    R = 1e-6;
    refstairs = -1;
    measured = 0;
    KoopmanStates = 6;
    PolyLiftingFunction = 0;% 2,4,6,10,12,14,18,24
end

if nargin < 6
    % Define the control to be used
    % 0. Run the simulation open loop
    % 1. NMPC Controller by TU Delft
    % 2. KLMPC/ qLMPC controller solved by quadprog solver
    controller = 2;
end

if nargin < 7
    ControlSetStr = 'sowfa_2turb_alm_turbl'; %'sowfa_2turb_yaw_alm_turbl';
    %sowfa_2turb_yaw_alm_combined1.mat sowfa_2turb_yaw_alm_turbl_AllComb
    %ControlSetStr = 'sowfa_2turb_yaw_alm_uniform';
    %ControlSetStr = 'sowfa_2turb_yaw_alm_combined';
    %ControlSetStr = 'sowfa_2turb_yaw_alm_combined1';
    %ControlSetStr = 'sowfa_2turb_yaw_alm_combined2'; %! Warning: Matrix is singular to working precision.
    %ControlSetStr = 'sowfa_2turb_yaw_alm_turbl_AllComb1'; %! Warning: Matrix is singular to working precision.
    %ControlSetStr = 'sowfa_2turb_yaw_alm_turbl_AllComb';
    %ControlSetStr = 'sowfa_2turb_yaw_steps';
    %ControlSetStr = 'sowfa_2turb_yaw_steps_Ct_comb';
    %ControlSetStr = 'sowfa_2turb_yaw_alm_turbl_AllComb';
    %ControlSetStr = 'sowfa_2turb_yaw_noise_step';'sowfa_2turb_yaw_noise2';
    %ControlSetStr = 'sowfa_2turb_yaw_noise_step3';
end

if nargin < 8
    Vinf = 8;
end

if nargin < 9 || isempty(vinfStr)
    vinfStr = '';%'_Vinf'; %_Vinf'; %alternativ '';
end

kVinf = 0;
if nargin <10
    kTimestep =  2;
    if refstairs
        kTimestep =  200;
    end
end

if nargin < 11
    sol_update = 0;
end

if nargin < 12
    OptYaw = 1;
end

if nargin < 13
    KoopAIC = 1;
end

sol_VinfInput = Vinf;

if nargin < 14
    %Rvec = [1,1,0.1]; %[1, 1,0.1];
    Rvec = [1,1,0.001];
end

if nargin < 15
    RU = 100; %0.75; 0
end

if nargin < 16
    KoopFeat = 0;
end

if nargin < 17
    AICnoV1 = 0;
end

if nargin < 18
    ad = [];
end


ControlSetStr2 = 'sowfa_2turb_yaw_alm_turbl_AllComb1';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%      WIND FARM SIMULATOR (WFSIM) by S. Boersma and B. Doekemeijer
%                 Delft University of Technology, 2018
%          Repo: https://github.com/TUDelft-DataDrivenControl/WFSim
%
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%
%%   Quick use:
%     1. Specify the wind farm you would like to simulate on line 73.
%        A list of all the wind farm scenarios can be found in
%        the 'layoutDefinitions' folder. You can also create your own
%        wind farm scenario here.
%     2. Either load a predefined timeseries of control inputs using line
%        74, or alternatively add your turbine/farm controllers in lines
%        115-119.
%     3. Setup the model solver settings in line 75. The default selection
%        is 'solverSet_default(Wp)', as defined in 'solverDefintions'.
%     4. Setup the simulation settings in lines 78-81.
%     5. Press start.
%
%%   Relevant input/output variables:
%     - modelOptions: this struct contains simulation settings
%     related to the wind farm itself (solution methodology, etc.)
%     - scriptOptions: this struct contains simulation settings, not
%     related to the wind farm itself (outputs, etc.)
%
%     - Wp: this struct contains all the simulation settings related to the
%           wind farm, the turbine inputs, the atmospheric properties, etc.
%         Wp.Nu:      Number of model states concerning longitudinal flow.
%         Wp.Nv:      Number of model states concerning lateral flow.
%         Wp.Np:      Number of model states concerning pressure terms.
%         Wp.sim:     Substruct containing timestep and simulation length.
%         Wp.turbine: Substruct containing turbine properties and settings.
%         Wp.site:    Substruct containing freestream atmospheric properties.
%         Wp.mesh:    Substruct containing topology and meshing settings.
%
%     - sol: this struct contains the system states at a certain timestep.
%         sol.k:     Discrete timestep  to which these system states belong
%         sol.time:  Actual time (in s) to which these system states belong
%         sol.x:     True system state (basically flow field excluding bcs)
%         sol.u:     Instantaneous longitudinal flow field over the mesh (in m/s)
%         sol.v:     Instantaneous longitudinal flow field over the mesh (in m/s)
%         sol.p:     Instantaneous pressure field over the mesh (in Pa)
%         sol.uu:    Same as sol.u, used for convergence
%         sol.vv:    Same as sol.v, used for convergence
%         sol.pp:    Same as sol.p, used for convergence
%         sol.turbine: a struct containing relevant turbine outputs such as
%         the ax. ind. factor, the generated power, and the ct coefficient
%         sol.turbInput: turbine inputs at current time
%
%     - sys: this struct contains the system matrices at a certain timestep.
%         sys.A:     System matrix A in the grand picture: A*sol.x = b
%         sys.b:     System vector b in the grand picture: A*sol.x = b
%         sys.pRCM:  Reverse Cuthill-McKee algorithm for solving A*x=b faster.
%         sys.B1:    Submatrix in the matrix sys.A.
%         sys.B2:    Submatrix in the matrix sys.A.
%         sys.bc:    Vector with boundary conditions of the continuity equation (part of sys.b).
%
%%   Debugging and contributing:
%     - First, try to locate any errors by turning all possible outputs
%       on (printProgress, printConvergence, Animate, plotMesh).
%     - If you cannot solve your problems, reach out on the Github.
%
%%   Nonlinear or stochastic MPC provinding active power control:
%     - In line 104 set if sol.k>=0 for open-loop (no control)
%       and sol.k<=0 for closed-loop (wind farm will track a power reference).
%     - Call in line 117 NMPCcontroller.m or SMPCcontroller.m for the nonlinear or stochastic controller, respectively.
%     - When the simulation is finished, run MPC_animation.m to see the results
%     - Note that to run a controller, installation of YALMIP and a solver
%       is required. Current controller settings belong to the sowfa_9turb_apc_alm_turbl
%       case and CPLEX 12.8 solver.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
profile ON
tic

%% Define simulation settings: layout, control inputs and simulation duration
[DataKoopman,DataOut,S] = WFSim_addpaths(controller,PolyLiftingFunction); % Add essential paths to MATLABs environment
eval(['Wp = layoutSet_sowfa_2turb_alm_turbl',vinfStr,';']); %alternativ; %layoutSet_sowfa_9turb_apc_alm_turbl(); % Choose which scenario to simulate. See 'layoutDefinitions' folder for the full list.

if ~isempty(Vinf) && length(Wp.site.u_Inf) == 1 %Vinf is set by the user
    Wp.site.u_Inf = Vinf(1); % Wp.turbine.Crx(2) = Wp.turbine.Crx(1) + 3 * Wp.turbine.Drotor;
end
%Wp.turbine.Crx(2) = Wp.turbine.Crx(1)+  diff(Wp.turbine.Crx)/5*5;

turbInputSet = eval(['controlSet_',ControlSetStr,'(Wp)']); % turbInputSet = controlSet_greedy(Wp);% for greedy control input
turbInputSetID =  eval(['controlSet_',ControlSetStr2,'(Wp)']);

WturbInputSetID.CT_prime = turbInputSetID.CT_prime - (mean(turbInputSetID.CT_prime'))';

nNoise = 10;
lNoise = ceil(length(turbInputSetID.CT_prime)/nNoise);
rng(0);
shN = randn(lNoise,1);
sHMat = repmat(shN',nNoise,1);
noiseVec1 = sHMat(:);

shN = randn(lNoise,1);
sHMat = repmat(shN',nNoise,1);
noiseVec2 = sHMat(:);

turbInputSetID.CT_prime = [noiseVec1'; noiseVec2'];

% figure; plot(noiseVec1(1:3600));´hold on; plot(noiseVec2(1:3600));
% figure; plot(turbInputSetID.CT_prime(1,1:3600));hold on; plot(turbInputSetID.CT_prime(2,1:3600));

modelOptions = solverSet_default(Wp); % Choose model solver options

% Simulation length, display and visualization settings
NN = floor(turbInputSet.t(end))/Wp.sim.h; %min(1200,floor(turbInputSet.t(end)/Wp.sim.h)); % Number of timesteps in simulation
verboseOptions.printProgress = 0;    % Print progress in cmd window every timestep. Default: true.
verboseOptions.Animate       = 0;    % Plot flow fields every [X] iterations (0: no plots). Default: 10.
verboseOptions.plotMesh      = 0;    % Plot mesh, turbine locations, and print grid offset values. Default: false.

%% Script core functions
[Wp,sol,sys] = InitWFSim(Wp,modelOptions,verboseOptions.plotMesh); % Initialize WFSim model
sol.VinfInput = sol_VinfInput;
sol.update = sol_update;

if length(Vinf) > 1 && length(Wp.site.u_Inf) == 1 %Vinf is set by the user
    VinfStr = strrep(sprintf('Vinf%2.1fTo%2.1f',Vinf(1),Vinf(end)),'.','dot'); % ramp [V0 V1] or time series
elseif length(Wp.site.u_Inf) == 1
    VinfStr = strrep(sprintf('Vinf%2.1f',Wp.site.u_Inf),'.','dot'); % information for plot
else
    VinfStr = strrep(sprintf('Vinf%2.1fNTM',mean(Wp.site.u_Inf)),'.','dot'); % information for plot
end
sol.kUinfChange = 0; %counter since change of Uinf %ToDo ADi remove this
sol.nkUinfChange = 0; %counter since this effected the grid points
sol.diffV = 0;

% Information for storage
if controller == 0
    folderName = sprintf('%s_OL_Ct',VinfStr);
    Rprint = ControlSetStr;
elseif controller == 1 || measured == 1
    folderName = sprintf('%s_CLmeasWind',VinfStr);
    Rprint = strrep(strrep(sprintf('R%2.1e',R),'.','dot'),'-','_neg');
else
    folderName = sprintf('%s_CL_K%02d_P%d',VinfStr,KoopmanStates,PolyLiftingFunction);%'Vinf%d_diffComb_states%d'
    Rprint = strrep(strrep(sprintf('R%2.1e',R),'.','dot'),'-','_neg');
end

dirFig = fullfile(DataOut,folderName);
if ~exist(dirFig,'dir')
    mkdir(dirFig);
end

if refstairs >= 1
    figName = sprintf('Stairs%d%s',refstairs,Rprint);
elseif OptYaw && ~measured
    figName = sprintf('%s_KoopAIC%d',Rprint,KoopAIC);
else
    figName = sprintf('%s',Rprint);
end

VinfStrFig = strrep(VinfStr,'To','');
fileName = [sprintf('%s_%s',VinfStrFig,ControlSetStr),'.mat'];
matFileExists = exist(fullfile(dirFig,fileName),'file') == 2;

if controller == 2 && ~measured && ~KoopAIC && OptYaw ~= 1 % Koopman AIC with yaw 0 or static yaw
    figName = sprintf('%s_KoopAIC0_OptYaw%d',Rprint,OptYaw);
end
if KoopFeat == 1
    figName = [figName,'_KoopFeat'];
elseif KoopFeat > 1
    figName = sprintf('%s_KoopFeat%d',figName,KoopFeat);
end
if AICnoV1
    figName = [figName,'_AICnoV1'];
    if ~isempty(ad) % settings of the online update other than in MPCinitv
        adFull = adaptSettings(ad);
        figName = sprintf('%s_nu%03d_h%d', figName, round(100*adFull.nu), adFull.h);
    end
end
filenamepng = matlab.lang.makeValidName(figName);
fullFigName = fullfile(dirFig,[filenamepng,'.fig']);

figFileExists = exist(fullFigName,'file') == 2;

plotDataFile = fullfile(dirFig,[filenamepng,sprintf('_Update%d',sol_update),'_plotData.mat']);
if controller == 2 && exist(plotDataFile,'file') == 2
    disp(['Plot data exist: ',plotDataFile]);
    tmp = load(plotDataFile,'pd');
    tmp.pd = addKoopFeatEstimate(tmp.pd, fileparts(DataKoopman)); % estimates of the input-only model
    plotWFSimTimeseriesTex(tmp.pd, dirFig, [filenamepng,sprintf('_Update%d',sol_update)]);
    sol_array = []; JR = tmp.pd.JR; JQ = tmp.pd.JQ;
    return;
end

if matFileExists && figFileExists %&& controller == 0 && 0%gt = get(gca, 'title');
    disp(['Mat file exists: ',fileName]);
    disp(['Fig file exists: ',filenamepng]);
    open(fullFigName);
    sol_array = [];JR = NaN; JQ = NaN;
    return;
end

% Initialize variables and figure specific to this script
CPUTime   = nan(NN,1);% Create empty matrix to save CPU timings
if verboseOptions.Animate > 0 % Create empty figure if Animation is on
    hfig = figure('color',[0 166/255 214/255],'units','normalized',...
        'outerposition',[0 0 1 1],'ToolBar','none','visible', 'on');
end

if controller > 0 && refstairs <= 0 % for all controllers
    NN = min(length(S.AGCdata(:,2)),NN);
elseif refstairs == 10
    %NN = 660;%
end

mpc.OptYaw = OptYaw;
mpc.KoopAIC = KoopAIC;
mpc.Rvec = Rvec;
mpc.KoopFeat = KoopFeat;
mpc.AICnoV1 = AICnoV1;
if ~isempty(ad), mpc.ad = adaptSettings(ad); end

mpc.RU = RU;

if controller == 1 %for original MPC controller
    mpc = MPCinit(sol,Wp,NN);
elseif controller == 2
    sol.mpc.OptYaw = 0;
    [mpc,K] = MPCinitv(sol,Wp,NN,refstairs,DataKoopman,KoopmanStates,R,PolyLiftingFunction,S,mpc);
    mpc.controller = controller;
end

% Performing forward time propagations
disp(['Performing ' num2str(NN) ' simulations..']);

if length(Wp.site.u_Inf) > 1
    Vinf = Wp.site.u_Inf;
    Vinfy = Wp.site.v_Inf;
end
Wp.site.u_Inf = Vinf(1);
if exist('Vinfy','var'), Wp.site.v_Inf  =  Vinfy(1); end

while sol.k < NN
    
    
    tic; % Start stopwatch
    
    % Determine control setting at current time by interpolation of time series
    turbInput = struct('t',sol.time);
    
    if ~controller || sol.k<=0 % sol.k>=0 % set to ">=0" for open-loop |or| "<=" for closed-loop -> changed to be able to run this as a function
        for i = 1:Wp.turbine.N
            % calculate the thrust coefficent and yaw angle with the help
            % of interpolation method from the given loop up table.
            turbInput.CT_prime(i,1) = interp1(turbInputSet.t,turbInputSet.CT_prime(i,:),sol.time,turbInputSet.interpMethod);
            turbInput.phi(i,1)      = interp1(turbInputSet.t,turbInputSet.phi(i,:),sol.time,turbInputSet.interpMethod);
        end
    else
        switch controller % decided controller will run under the close-loop configuration.
            case 1
                % uses YALMIP toolbox and its optimizer function
                [turbInput.CT_prime(:,1),turbInput.phi(:,1),mpc] = NMPCcontroller(sol,Wp,mpc);
            case 2
                % MATLAB quadprog (QP cost fun: min f(x) = 1/2*x'Hx+f'x s.t. Ax <= b
                % measured = 0: Koopman model for wind
                [turbInput.CT_prime(:,1),turbInput.phi(:,1),mpc] = qLPVMPCcontroller_dU_quadprog(sol,Wp,mpc,K,measured);
            otherwise
                % not yet implemented -> this should stop
                error('Controller switch has to be set to 0,1 or 2')
        end
    end
    
    % Propagate the WFSim model to the next timestep.
    if sol.k >= kTimestep && length(Vinf) == 2 && Wp.site.u_Inf < Vinf(2) %&& length(Vinf) <= 2
        kVinf = kVinf + 1;
        Wp.site.u_Inf = Vinf(1) + kVinf*0.01; %Vinf(2); %
    elseif length(Vinf) > 2
        Wp.site.u_Inf = Vinf(max(1,sol.k)); %Vinf(1) + kVinf*0.01;
        if exist('Vinfy','var') % time series without lateral component: v_Inf unchanged
            Wp.site.v_Inf = Vinfy(max(1,sol.k));
        end
    end
    
    [sol,sys]      = WFSim_timestepping(sol,sys,Wp,turbInput,modelOptions); % forward timestep: x_k+1 = f(x_k)
    CPUTime(sol.k) = toc;   % Stop stopwatch
    if controller > 1
        sol.xprev = mpc.xprev;
        sol.error_norm = mpc.error_norm;
        if sol.k >1
            sol.P1 = mpc.P1;
        else
            sol.P1 = 0;
        end
    end
    % Save sol to cell array
    sol_array(sol.k) = sol; %#ok<AGROW> this is not a huge delay
    
    % Print progress, if necessary
    if verboseOptions.printProgress
        disp(['Simulated t(' num2str(sol.k) ') = ' num2str(sol.time) ...
            ' s. CPU: ' num2str(CPUTime(sol.k)*1e3,3) ' ms.']);
    end
    
    % Plot animations, if necessary
    if verboseOptions.Animate > 0
        if ~rem(sol.k,verboseOptions.Animate)
            hfig = WFSim_animation(Wp,sol,hfig);
        end
    end
end

disp(' ')
disp(['Completed ' num2str(NN) ' forward simulations. Average CPU time: ' num2str(mean(CPUTime)*10^3,3) ' ms.']);
t_display = toc;
fprintf('Lemke qLPVMPC solved in %2.3f \n',t_display);
profile off

%% plot figures and save data
[time, CT_prime, Phi, Power, ~, V, JR] = saveDataEDMD(Wp,sol_array,controller,figName,dirFig,240);

if controller == 0
    time1 = time(time >=500);
    [JR,JQ] = plotWFSimTimeseries(Wp,sol_array,controller,Power,CT_prime,Phi,V,JR,time,measured,figName,dirFig);
elseif controller == 1
    [JR,JQ] = plotWFSimTimeseries(Wp,sol_array,controller,Power,CT_prime,Phi,V,JR,time,measured,figName,dirFig,mpc);
else
    time1 = time(time >= 500 & time <= 1200);
    if refstairs == 20 % adaptive AIC test case: whole simulation
        time1 = time(time >= 1);
    end
    [JR,JQ] = plotWFSimTimeseriesLatex2(Wp,sol_array,controller,Power,CT_prime,Phi,V,JR,time1,measured,figName,dirFig,mpc);

    % Save the plot data of this run (controller and model in the file name)
    kk = length(Power);
    utemp = reshape(cell2mat(arrayfun(@(x)(x.u),sol_array,'UniformOutput',false)),Wp.mesh.Nx,Wp.mesh.Ny,kk);
    nTxP = round(Wp.turbine.Crx(1)/Wp.mesh.Lx*Wp.mesh.Nx);
    nTyP = round(Wp.turbine.Cry(1)/Wp.mesh.Ly*Wp.mesh.Ny);
    X_est = cell2mat(arrayfun(@(x)(x.xprev),sol_array,'UniformOutput',false));
    pd = struct('time',time1,'Power',Power,'CT_prime',CT_prime,'Phi',Phi,'V',V, ...
        'V_est',X_est(1:2,:),'u1',squeeze(utemp(nTxP,nTyP,:))','uInf',squeeze(utemp(1,nTyP,:))', ...
        'Pref',mpc.Pref(1:kk),'Pgreedy',mpc.Pgreedy,'JQ',JQ,'JR',JR,'R',R, ...
        'controller',controller,'measured',measured,'OptYaw',OptYaw,'KoopAIC',KoopAIC, ...
        'KoopFeat',KoopFeat,'KoopmanStates',KoopmanStates,'ControlSetStr',ControlSetStr);
    if AICnoV1 % adaptive AIC: update flags, title, no greedy power (8 m/s)
        P1 = cell2mat(arrayfun(@(x)(x.P1),sol_array,'UniformOutput',false));
        pd.update = min(max(P1-1,0),1);
        pd.modeStr = 'Non-adaptive AIC';
        if sol_update, pd.modeStr = 'Adaptive AIC'; end
        pd.showGreedy = false;
        pd.showYaw = false; % thrust control only
        pd.AICnoV1 = AICnoV1;
        pd.sol_update = sol_update;
        if ~isempty(ad) % settings of the online update in the title
            pd.ad = mpc.ad;
            pd.modeStr = sprintf('AIC (\\nu = %.2f m/s, h = %d)', mpc.ad.nu, mpc.ad.h); % short: fits the width
        end
    end
    pd = addKoopFeatEstimate(pd, fileparts(DataKoopman)); % estimates of the input-only model
    save(plotDataFile,'pd');
    plotWFSimTimeseriesTex(pd, dirFig, [filenamepng,sprintf('_Update%d',sol_update)]);
end