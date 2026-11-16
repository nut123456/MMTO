function Cantilever3D_MMTO_READY(nelx,nely,nelz,volfrac,rmin,nl,cgtol,cgmax,filename,Evec)
% V3.4 AL material-volume control (subcell fractions) + truthful rendering.
% Multi-material 3-D cantilever topology optimization: compliance minimization with a VOLUME constraint
% CLUSTERED-MATERIAL version: differentiable correlation-length-inspired CMF
% gives coherent, clearly separated material domains inside the optimization.
% - Single design field x in [0,1]
% - SVMM projection -> hierarchical partition-of-unity interpolation
% - MG-preconditioned CG (MGCG) is the optimization solver
% - Solver benchmark at the final design: CG, Jacobi-PCG, IC-PCG, MGCG (optional direct)
% - VFLS-guided ADAPTIVE SUBCELL integration is used INSIDE every FE iteration
%   for BOTH stiffness and physical volume, while the global structured FE DOFs are unchanged
% - Normalized Heaviside projections enforce h(0)=0 and h(1)=1 so that
%   the material partition and the volume constraint use the same physical representation
% - Neighborhood-supported minimum-feature control (MFC) is embedded in the analysis field
% - Clustered material field (CMF) is embedded in FE/volume/sensitivity analysis
%   so thin/isolated members are suppressed before adaptive subcell FE integration
% - TRUE non-smoothed optimized material-domain visualization and 6-face connectivity diagnostics
%
% READY-TO-RUN DEFAULT (press Run or call with no arguments):
%   mesh = 48 x 24 x 24
%   volfrac = 0.12 (full cantilever design domain)
%   3 materials E=[0.25 0.70 1.00]
% - Constraint: mean design volume <= volfrac; NO mass/resource constraint
%
% Equivalent explicit call:
% Cantilever3D_VOL_MGCG_VFLS_Adaptive_MFC_CLUSTERED_FAST(48,24,24,0.12,2.5,1,3,1e-6,1000, ...
%     'Cantilever3D_VOL_VFLSadaptive_MFC_CLUSTERED_48x24x24.mat',[0.25 0.70 1.00]);
if nargin ~= 0 && nargin ~= 10
    error('Use either no inputs (default 48x24x24 3-D Cantilever case) or all 10 input arguments.');
end
if nargin == 0
    nelx     = 48;
    nely     = 24;
    nelz     = 24;
    volfrac  = 0.12;
    rmin     = 2.5;
    nl       = 3;
    cgtol    = 1e-3;
    cgmax    = 800;
    filename = sprintf('Cantilever3D_V31_INTERFACE_CONNECTIVITY_%dx%dx%d.mat',nelx,nely,nelz);
    %Evec     = [0.50 1.00];
    Evec     = [0.25 0.70 1.00];
    % Evec     = [0.25 0.50 0.70 0.85 1.00];
end
% MG hierarchy requires integer coarsening at every level
for ll = 1:nl-1
    div = 2^ll;
    assert(mod(nelx,div)==0 && mod(nely,div)==0 && mod(nelz,div)==0, ...
        'nelx, nely, nelz must be divisible by 2^(nl-1) for this MG hierarchy.');
end
tic;
%% ===== (A) MATERIAL / CONTINUATION =====
Emin = 1e-9;                 % void stiffness (normalized)
nu   = 0.30;
% continuation: [start, max, start_iter, period, step]
penalCn = {1,  4, 20, 5, 0.12};    % FAST-V2: smoother p continuation, avoids large compliance spikes
betaCn  = {1, 32, 25, 5, 0.25};    % FAST-V2: smoother beta continuation, reaches beta=32 progressively
penal = penalCn{1};
beta  = betaCn{1};
% sanity check materials
mnum = numel(Evec);
assert(mnum>=2, 'Evec must contain at least two material stiffness values.');
% IMPORTANT: there is deliberately NO density/mass vector Dvec in this version.
% The only global resource constraint is the design-volume fraction volfrac.
%% ===== (B) FE PREP =====
KE = Ke3D(nu);
nelem = nelx*nely*nelz;
nx = nelx+1; ny = nely+1; nz = nelz+1;
ndof = 3*nx*ny*nz;
% node numbering & edof
nodenrs(1:ny,1:nz,1:nx) = reshape(1:ny*nz*nx,ny,nz,nx);
edofVec(1:nelem,1) = reshape(3*nodenrs(1:ny-1,1:nz-1,1:nx-1)+1,nelem,1);
edofMat(1:nelem,1:24) = repmat(edofVec(1:nelem),1,24) + ...
    repmat([0 1 2 3*ny*nz+[0 1 2 -3 -2 -1] -3 -2 -1 ...
            3*ny+[0 1 2] 3*ny*(nz+1)+[0 1 2 -3 -2 -1] 3*ny+[-3 -2 -1]],nelem,1);
iK = reshape(kron(edofMat,ones(24,1))',576*nelem,1);
jK = reshape(kron(edofMat,ones(1,24))',576*nelem,1);
% -------------------------------------------------------------------------
% Element corner nodes + element-to-node averaging operator.
% These are used by the VFLS-guided adaptive subcell integration.
% IMPORTANT: the global FE mesh/DOFs remain unchanged; only the constitutive
% integration inside interface elements is locally refined.
% Node order for scalar trilinear interpolation:
%   1:(y-,z-,x-)  2:(y+,z-,x-)  3:(y+,z+,x-)  4:(y-,z+,x-)
%   5:(y-,z-,x+)  6:(y+,z-,x+)  7:(y+,z+,x+)  8:(y-,z+,x+)
% -------------------------------------------------------------------------
elemNodes = [ ...
    reshape(nodenrs(1:nely,   1:nelz,   1:nelx),[],1), ...
    reshape(nodenrs(2:nely+1, 1:nelz,   1:nelx),[],1), ...
    reshape(nodenrs(2:nely+1, 2:nelz+1, 1:nelx),[],1), ...
    reshape(nodenrs(1:nely,   2:nelz+1, 1:nelx),[],1), ...
    reshape(nodenrs(1:nely,   1:nelz,   2:nelx+1),[],1), ...
    reshape(nodenrs(2:nely+1, 1:nelz,   2:nelx+1),[],1), ...
    reshape(nodenrs(2:nely+1, 2:nelz+1, 2:nelx+1),[],1), ...
    reshape(nodenrs(1:nely,   2:nelz+1, 2:nelx+1),[],1)];
nodeElemInc = sparse(elemNodes(:),repmat((1:nelem)',8,1),1,ny*nz*nx,nelem);
nodeElemDeg = full(sum(nodeElemInc,2));
nodeElemDeg(nodeElemDeg==0) = 1;
Wne = spdiags(1./nodeElemDeg,0,ny*nz*nx,ny*nz*nx)*nodeElemInc;
% Multigrid prolongation operators
Pu=cell(nl-1,1);
for l = 1:nl-1
    Pu{l,1} = prepcoarse(nelz/2^(l-1),nely/2^(l-1),nelx/2^(l-1));
end
%% ===== 3-D CANTILEVER: LOADS / SUPPORTS / PASSIVE SETS =====
U = zeros(ndof,1);
% Element numbering: (y,z,x)

% Fixed support: complete left face, x = 1
fixedNodes = double(nodenrs(:,:,1));
fixeddofs = unique([3*fixedNodes(:)-2; ... % ux
                    3*fixedNodes(:)-1; ... % uy
                    3*fixedNodes(:)]);      % uz
% Sine-distributed downward load on bottom edge of free/right face:
% x = nelx+1, z = 1, y = 1...nely+1
loadNodes = double(nodenrs(:,1,nelx+1));
loadDofs  = 3*loadNodes(:);
loadVals  = -sin((0:nely)'/nely*pi);
F = sparse(loadDofs,1,loadVals,ndof,1);
F(fixeddofs) = 0;
% Entire cantilever domain is designable
passiveS = false(nelem,1);
passiveV = false(nelem,1);
passiveSolid = passiveS;
passiveVoid  = passiveV;
activeDesign = ~(passiveVoid | passiveSolid);
nActive = sum(activeDesign);
% Null-space elimination for fixed DOFs
N = ones(ndof,1);
N(fixeddofs) = 0;
Null = spdiags(N,0,ndof,ndof);
% Display coordinates for supports and distributed load
[iyFix,izFix] = ndgrid(1:ny,1:nz);
supportXYZ = [ones(numel(iyFix),1), iyFix(:), izFix(:)];
loadPlotX = (nelx+1)*ones(size(loadNodes));
loadPlotY = (1:ny)';
loadPlotZ = ones(size(loadNodes));
loadPlotW = loadVals;
fprintf('\n============================================================\n');
fprintf(' 3-D CANTILEVER MMTO + VFLS-ADAPTIVE + MINIMUM-FEATURE CONTROL + MGCG\n');
fprintf(' Mesh                 : %d x %d x %d\n',nelx,nely,nelz);
fprintf(' Elements             : %d\n',nelem);
fprintf(' DOFs                 : %d\n',ndof);
fprintf(' Fixed face nodes     : %d\n',numel(fixedNodes));
fprintf(' Sine-load edge nodes : %d\n',numel(loadNodes));
fprintf(' Total Fz             : %.6e\n',sum(loadVals));
fprintf(' Active design elems  : %d\n',nActive);
fprintf(' Passive solid elems  : %d\n',nnz(passiveSolid));
fprintf(' Passive void elems   : %d\n',nnz(passiveVoid));
fprintf(' Volume fraction      : %.4f (FULL CANTILEVER DESIGN DOMAIN)\n',volfrac);
fprintf('============================================================\n\n');
%% ===== (C0) NEIGHBORHOOD-SUPPORTED MINIMUM-FEATURE CONTROL =====
% A second, wider neighborhood filter is used to measure local support.
% The filtered design xTilde is transformed to an analysis field
%
%   q_e      = (G*xTilde)_e / Gs_e
%   s_e      = H_N(q_e; eta_N,beta_N)
%   xMFC_e   = xTilde_e * [(1-alpha)+alpha*s_e]
%
% where alpha is continued from 0 to maxAlpha. Thin/isolated members have
% low neighborhood support q_e and are therefore attenuated. Stiffness,
% physical volume, adaptive-subcell classification, and their gradients all
% use xMFC, so this is an optimization/analysis control rather than display
% smoothing. The global structured FE mesh and MGCG hierarchy are unchanged.
mfcCfg.enabled       = true;
mfcCfg.radius        = 2.50;    % support radius in parent-element units
mfcCfg.etaNeighbor   = 0.060;   % neighborhood-support threshold
mfcCfg.betaNeighbor  = 24.0;    % smooth support projection sharpness
mfcCfg.startIter     = 20;      % FAST: activate MFC earlier
mfcCfg.fullIter      = 90;      % FAST: reaches maxAlpha earlier
mfcCfg.maxAlpha      = 0.90;    % 0=no MFC, 1=full support gating
mfcCfg.reportEvery   = 10;
fprintf(' Minimum-feature control : ON, radius=%.2f, etaN=%.3f, alphaMax=%.2f\n', ...
    mfcCfg.radius,mfcCfg.etaNeighbor,mfcCfg.maxAlpha);
%% ===== (C1) CORRELATION-LENGTH-INSPIRED CLUSTERED MATERIAL FIELD =====
% The cited MFSE work obtains coherent material regions using smooth spatial
% material fields with a controllable correlation length. This implementation
% keeps the present SVMM/MMTO representation, but introduces a differentiable
% clustered material field inside FE analysis and the sensitivity chain:
%   qM = (Gmat*xMFC)./GmatS
%   T(qM) = sum_r w_r*c_r,   w_r ~ exp[-kappa*(qM-c_r)^2]
%   xCMF = (1-alphaM)*xMFC + alphaM*T(qM).
% A larger cmfCfg.radius produces larger, cleaner MAT domains.
cmfCfg.enabled        = true;
cmfCfg.radius         = 3.0;  % raise to 4.0 only if more RAM is available
cmfCfg.startIter      = 10;
cmfCfg.fullIter       = 90;
cmfCfg.maxAlpha       = 0.92;
cmfCfg.kappaStart     = 8.0;
cmfCfg.kappaMax       = 35.0;
cmfCfg.kappaFullIter  = 120;
cmfCfg.reportEvery    = 10;
fprintf(' Clustered material field : ON, radius=%.2f cells, alphaMax=%.2f, kappaMax=%.1f\n', ...
    cmfCfg.radius,cmfCfg.maxAlpha,cmfCfg.kappaMax);
% Material-volume balancing. With only a total volume constraint, pure
% compliance minimization would prefer the stiffest phase everywhere. The
% cited three-material example prescribes each material amount. Here the
% default targets are equal and sum to volfrac; change targetFrac as desired.
matBalCfg.enabled       = true;
matBalCfg.targetFrac    = (volfrac/mnum)*ones(1,mnum);
matBalCfg.startIter     = 15;
matBalCfg.fullIter      = 90;
matBalCfg.maxWeight     = 0.35;
fprintf(' Material target fractions : %s (sum = %.4f)\n', ...
    mat2str(matBalCfg.targetFrac,4),sum(matBalCfg.targetFrac));
% V3.4: augmented-Lagrangian targets, updated after each accepted OC step.
% Relative residual g_j=(V_j-target_j)/target_j. Total volume remains
% enforced separately by the existing adaptive-subcell OC corrector.
% Because the M phase volumes sum to Vtotal, enforce only M-1 independent
% targets; the final material is implied by total volume equality.
alCfg.enabled = true;
alCfg.startIter = 45;
alCfg.updateEvery = 5;
alCfg.rho = 0.75;
alCfg.rhoGrowth = 1.30;
alCfg.rhoMax = 120;
alCfg.lambdaMax = 250;
alCfg.targetRelTol = 0.025;  % 2.5% tolerance on individual material fractions
alCfg.mu = zeros(1,mnum-1);
alCfg.previousError = inf;
alCfg.stallUpdates = 0;
assert(all(matBalCfg.targetFrac > 0) && ...
       abs(sum(matBalCfg.targetFrac)-volfrac) < 1e-10, ...
       'Individual material targets must sum to volfrac.');
fprintf(' V3.4 AL material volumes: independent constraints=%d, relative tolerance=%.1f%%\n', ...
    mnum-1,100*alCfg.targetRelTol);

%% ===== V3.3: TRUE OPTIMIZATION INTERFACE / LOCAL CONNECTIVITY =====
% All regularizers act on physical hierarchical phase fields psiV, never
% on a rendered mask. Their analytic gradients pass through CMF, MFC,
% and filtering via the same transpose chain as compliance.
regCfg.enabled = true;
regCfg.startIter = 35;
regCfg.fullIter  = 125;
regCfg.interfaceWeight = 0.045; % dimensionless augmented-compliance weight
regCfg.connectWeight   = 0.060;
regCfg.epsNeighbor = 1e-3;      % differentiable local solid support
regCfg.connRadius = 1;          % direct 6-face neighbors; local surrogate only
regCfg.gradientCheck = true;
regCfg.thinWeight = 0.025;        % two-hop same-material neighborhood surrogate
regCfg.balanceSubcell = true;    % integrate materials at same sites as stiffness
regCfg.thinStartIter = 65;    % finite-difference check at initial iteration
% Symmetric 6-face adjacency: element ordering (y,z,x)
[Aface,nFaces] = mmto_face_adjacency(nely,nelz,nelx);
faceDegree = full(sum(Aface,2));
faceDegreeSafe = max(faceDegree,1);
% Two-hop face graph: intended to reduce thin and isolated material stripes.
A2 = spones(Aface*Aface);
A2 = A2 - spdiags(diag(A2),0,nelem,nelem);
A2 = spones(A2+Aface);
degree2 = max(full(sum(A2,2)),1);
fprintf(' V3.3 local interface/connection surrogate: %d undirected faces\n',nFaces);
fprintf(' Interface penalty=%.4f, local support penalty=%.4f\n', ...
    regCfg.interfaceWeight,regCfg.connectWeight);
%% ===== (C) DENSITY FILTER =====
[iH,jH,sH] = deal(zeros(nelem*(2*(ceil(rmin)-1)+1)^3,1));
k = 0;
for i1 = 1:nelx
  for k1 = 1:nelz
    for j1 = 1:nely
      e1 = (i1-1)*nely*nelz + (k1-1)*nely + j1;
      for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
        for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),nelz)
          for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
            e2 = (i2-1)*nely*nelz + (k2-1)*nely + j2;
            k=k+1; iH(k)=e1; jH(k)=e2;
            sH(k)=max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
          end
        end
      end
    end
  end
end
H  = sparse(iH(1:k),jH(1:k),sH(1:k));
Hs = sum(H,2);
% Wider neighborhood operator used only by minimum-feature control.
[G,Gs] = build_scalar_filter_3d(nelx,nely,nelz,mfcCfg.radius);
% Separate material-correlation operator for coherent MAT domains.
[Gmat,GmatS] = build_scalar_filter_3d(nelx,nely,nelz,cmfCfg.radius);
%% ===== (D) INIT =====
% SVMM thresholds
tgrid = (1:mnum)/(mnum+1);
% For three materials: centers = [0, 0.375, 0.625, 0.875].
cmfCenters = zeros(1,mnum+1);
cmfCenters(1) = 0.0;
for jj = 1:mnum-1
    cmfCenters(jj+1) = 0.5*(tgrid(jj)+tgrid(jj+1));
end
cmfCenters(end) = 0.5*(tgrid(end)+1.0);
cmfCfg.centers = cmfCenters;
fprintf(' CMF phase centers         : %s\n',mat2str(cmfCfg.centers,3));
% Initialize the design so that the NORMALIZED physical solid indicator
% rho=h_1 is already equal to volfrac for a uniform field.  Using x=volfrac
% directly is not volume-consistent when the Heaviside threshold eta_1 is
% nonzero (for 3 materials eta_1=0.25).
x0 = inverse_normalized_heaviside(volfrac,tgrid(1),beta);
x = x0*ones(nelem,1);
x(passiveVoid)  = 1e-10;
x(passiveSolid) = 1.0;

loop = 0; change = inf; 
volCurrent = inf;
Stats = zeros(300,4);
% ===== History for convergence plot =====
maxLoop = 500;
HistComp       = nan(maxLoop,1);
HistVolTotal   = nan(maxLoop,1);
HistMatVolFrac = nan(maxLoop,mnum);
HistAdaptCount = nan(maxLoop,3);   % number of elements using 1^3, 2^3, 4^3 subcells
HistMFCAlpha   = nan(maxLoop,1);
HistMFCSupport = nan(maxLoop,1);   % mean neighborhood support projection
HistCMFAlpha   = nan(maxLoop,1);
HistCMFKappa   = nan(maxLoop,1);
HistMatPurity  = nan(maxLoop,1);
HistInterface = nan(maxLoop,1);
HistConnection = nan(maxLoop,1);
HistRegWeight = nan(maxLoop,2);
HistALResidual = nan(maxLoop,mnum);
HistALPenalty = nan(maxLoop,1);
HistALMultipliers = nan(maxLoop,mnum-1);
% Robust stopping controls.  The old condition change>1e-2 could stop after
% Iteration 1 when the OC move limit was exactly 0.01, so the 2^3/4^3
% adaptive stages were never reached.
maxOptIter = 420;
minOptIter = 120;       % FAST: enough iterations to enter selected 4^3 stage
changeTol  = 5.0e-3;
volTolAbs  = 5.0e-3;
% FAST early-stop based on relative compliance stagnation
stableCount = 0;
stableNeed  = 12;
relCompTol  = 2.0e-4;
prevCompForStop = inf;
% %% ===== GIF setup (before MAIN LOOP) =====
gifname = sprintf('Cantilever3D_MM_%dx%dx%d.gif', nelx, nely, nelz);
if exist(gifname, 'file'), delete(gifname); end   % ลบไฟล์เก่าถ้ามี
% เตรียม figure ล่วงหน้า
figure(2);
clf;
set(gcf, 'Color', 'w', 'Renderer', 'opengl', 'Position', [100 100 960 540]);
%% ===== VIDEO setup (before MAIN LOOP) =====
vname = sprintf('Cantilever3D_MM_%dx%dx%d.mp4', nelx, nely, nelz);
if exist(vname,'file'), delete(vname); end
v = VideoWriter(vname, 'MPEG-4');   % ใช้ H.264 codec
v.Quality   = 100;                   % ความคมชัด (0–100)
v.FrameRate = 4;                    % 🔹ช้าลง (3–6 = เหมาะสุดสำหรับ topology animation)
open(v);
% สร้าง figure สำหรับแสดงผลคงที่
figure(2); clf;
set(gcf,'Color','w','Renderer','opengl','Position',[100 100 960 540]);
%% ===== Save figures folder =====
figDir = ['Results_' datestr(now,'yyyy_mm_dd_HH_MM_SS')];
if ~exist(figDir,'dir')
    mkdir(figDir);
end
%% ===== CSV optimization-history setup =====
% Saves one row per optimization iteration. The file is closed after every
% row so that the history remains available even if MATLAB is interrupted.
iterCSV = fullfile(figDir,'Optimization_History.csv');
csvFID = fopen(iterCSV,'w');
if csvFID < 0
    error('Cannot create optimization-history CSV: %s',iterCSV);
end
fprintf(csvFID,['Iteration,Cfull,Vphys,Vtarget,Lambda,Change,MGIterations,MGRelRes,' ...
                'Penal,Beta,MFCAlpha,CMFAlpha,CMFKappa,Purity,MBWeight,MBError,' ...
                'Sub1,Sub2,Sub4,MeanMFCSupport,InterfacePenalty,LocalConnectivityPenalty,InterfaceWeight,ConnectivityWeight']);
for jj = 1:mnum
    fprintf(csvFID,',MaterialVol%d',jj);
end
fprintf(csvFID,'\n');
fclose(csvFID);
fprintf(' Optimization history CSV : %s\n',iterCSV);
%% ===== VFLS-GUIDED ADAPTIVE SUBCELL ANALYSIS =====
% This section changes the ANALYSIS, not only the picture.
%
% Let xTilde be the filtered design field and t_k the hierarchical material
% thresholds.  Each material interface is represented implicitly by
%
%       Phi_k = xTilde - t_k = 0.
%
% Elements crossed by Phi_k=0 are integrated with 4x4x4 subcell samples;
% elements close to an interface use 2x2x2 samples; regular bulk elements
% retain one constitutive sample.  The global Hex8 FE mesh, DOFs, geometric
% multigrid hierarchy, and MGCG solver therefore remain unchanged.
%
% The subcell integration evaluates the SAME smooth Heaviside/material law
% at local trilinear samples and supplies a consistent chain-rule gradient
% back to x through nodal averaging and the density filter.
adaptCfg.enabled          = true;
adaptCfg.level2Band       = 0.070;   % conservative near-interface band
adaptCfg.level4Band       = 0.020;   % narrow high-resolution band
adaptCfg.spanLevel2       = 0.140;   % moderate nodal-variation trigger
adaptCfg.spanLevel4       = 0.180;   % high-resolution variation trigger
adaptCfg.highSubdiv       = 4;       % final-stage value; continuation below
adaptCfg.midSubdiv        = 2;
adaptCfg.bulkSubdiv       = 1;
adaptCfg.reportEvery      = 10;
% Adaptive-resolution continuation:
%   Iter   1-40 : bulk 1^3, interface 2^3
%   Iter  41-90 : bulk 1^3, near/cut 2^3
%   Iter  91-end: bulk 1^3, near 2^3, selected cut/high-gradient 4^3
fprintf(' VFLS adaptive subcells : ON with 1^3 -> 2^3 -> selected 4^3 continuation\n');
fprintf(' Physical volume         : SAME adaptive subcells as stiffness\n');
fprintf(' Minimum feature         : neighborhood-supported field used by FE + volume\n');
fprintf(' Stop rule               : min %d iters, max %d, stagnation + |dV|<=%.1e\n', ...
    minOptIter,maxOptIter,volTolAbs);
%% ===== VFLS-INSPIRED DISPLAY BOUNDARY REGULARIZATION =====
% The optimizer/FE fields are NOT modified.  For visualization only, each
% material field is converted to an implicit level-set Phi = V-lvl, then:
%   (1) trilinear grid refinement,
%   (2) VFLS-style diffusion of Phi (curvature regularization),
%   (3) signed-distance reinitialization using the 3-D upwind scheme,
%   (4) optional volume-preserving zero-level correction,
%   (5) very light Taubin smoothing of the extracted triangle mesh.
% This follows the boundary-regularization ideas of the 3-D VFLS method,
% while preserving the MMTO + MGCG optimization itself.
renderCfg.isoUpsample       = 2;      % 2 for live plots; use 3 for final/high-quality
renderCfg.vflsDiffusePasses = 2;      % diffusion/reinitialization cycles
renderCfg.vflsDiffusion     = 0.035;  % explicit Laplacian coefficient (< 1/6 is conservative)
renderCfg.vflsReinitIter    = 14;     % pseudo-time reinitialization steps/pass
renderCfg.vflsReinitDt      = 0.10;   % CFL-like pseudo-time step
renderCfg.preserveVolume    = true;   % approximately preserve V>=lvl fraction
renderCfg.taubinIter        = 2;      % only a light final mesh cleanup
renderCfg.lambda            = 0.35;
renderCfg.mu                = -0.36;
renderCfg.freezeBoundary    = true;   % keep open/support-side mesh boundaries fixed
renderCfg.liveEvery         = 2;      % render combined views every N iterations
renderCfg.saveEvery         = 10;     % save PNGs every N iterations
%% ===== (E) MAIN LOOP =====
while loop < maxOptIter && ...
      (loop < minOptIter || change > changeTol || ...
       ~isfinite(volCurrent) || abs(volCurrent-volfrac) > volTolAbs)
  loop = loop + 1;
  % Parameters actually used to analyze/update this iteration.  These are
  % retained for a final system that is exactly consistent with the last OC step.
  penalUsed = penal;
  betaUsed  = beta;
  % ---- E.1 Filter + Heaviside projection
  assert(isvector(x) && numel(x)==nelem, ...
      'Design vector x was corrupted: expected %d entries, got size %s.', ...
      nelem, mat2str(size(x)));
  x = x(:);
  xTilde = (H*x)./Hs;
  % ---- E.1b Minimum-feature analysis field
  mfcIter = mfc_cfg_for_iter(mfcCfg,loop);
  [xMFC,mfcData] = minimum_feature_field(xTilde,G,Gs,mfcIter);
  mfcIterLast = mfcIter;
  % ---- E.1c Correlation-length material clustering
  cmfIter = cmf_cfg_for_iter(cmfCfg,loop);
  [xCMF,cmfData] = clustered_material_field(xMFC,Gmat,GmatS,cmfIter);
  xCMF(passiveVoid)  = 0.0;
  xCMF(passiveSolid) = 1.0;
  cmfIterLast = cmfIter;
  [h,dh] = heaviside_smooth(xCMF,tgrid,betaUsed);   % ne×m
  % ---- E.2 Multi-material fields + analysis-consistent adaptive subcells
  % Stiffness, physical volume, and VFLS interface classification all use
  % the SAME clustered material field xCMF.
  [psiV,dpsiVCell] = MatIntModel(h,dh,1);          % normalized material-volume fields
  %[psiE,~] = MatIntModel(h,dh,penalUsed);          % compatibility/output field
  % V3.3: evaluate material fractions at the SAME adaptive constitutive
  % samples as FE stiffness and physical occupancy; use exact local pullback.
  adaptIter = adaptive_cfg_for_iter(adaptCfg,loop);
  [Ee,Ve,adaptData] = vfls_adaptive_subcell_response( ...
      xCMF,elemNodes,Wne,tgrid,betaUsed,penalUsed,Evec,Emin,adaptIter);
  adaptIterLast = adaptIter;
  Psub = adaptData.Psi;   % nelem x mnum; sum(Psub,2) == Ve
  matBalW = material_balance_weight(matBalCfg,loop);
  matFracNow = sum(Psub(activeDesign,:),1)/nActive;
  matDev = (matFracNow-matBalCfg.targetFrac)./max(matBalCfg.targetFrac,1e-12);
  matMismatch = sum(matDev.^2);
  % AL is dimensionless and applied to M-1 independent phase volumes.
  % Shifted Lagrangian: sum mu*g + rho/2*sum(g^2).
  % The final phase is implied by sum(Psi)=Ve, not separately penalized.
  alReady = alCfg.enabled && loop >= alCfg.startIter;
  if alReady
      alResidual = matDev(1:mnum-1);
      alValue = sum(alCfg.mu.*alResidual) + ...
                0.5*alCfg.rho*sum(alResidual.^2);
  else
      alResidual = zeros(1,mnum-1);
      alValue = 0;
  end
  gBalP = zeros(nelem,mnum);
  for jj=1:mnum
      coeff = matBalW*2*matDev(jj)/max(matBalCfg.targetFrac(jj),1e-12);
      if alReady && jj<mnum
          coeff = coeff + (alCfg.mu(jj)+alCfg.rho*alResidual(jj)) / ...
                    max(matBalCfg.targetFrac(jj),1e-12);
      end
      gBalP(activeDesign,jj) = coeff/nActive;
  end
  % Interface, one-hop same-material support and two-hop thin-domain control
  % are all applied to physical subcell-integrated material fractions.
  [pInterface,pConnect,gRegPsiI,gRegPsiC] = ...
      mmto_v31_penalties(Psub,Aface,faceDegreeSafe,activeDesign,regCfg.epsNeighbor);
  [~,pThin,~,gThin] = ...
      mmto_v31_penalties(Psub,A2,degree2,activeDesign,regCfg.epsNeighbor);
  rw = mmto_v31_reg_weights(regCfg,loop);
  thinW = regCfg.thinWeight * max(0,min(1,(loop-regCfg.thinStartIter)/ ...
      max(1,regCfg.fullIter-regCfg.thinStartIter)));
  if loop==1 && regCfg.gradientCheck
      mmto_v31_gradient_check(xCMF,tgrid,betaUsed,Aface,faceDegreeSafe,activeDesign,regCfg);
  end
  % Pull back derivative from subcell phase fractions to scalar CMF.
  gRegP = gBalP + rw(1)*gRegPsiI + rw(2)*gRegPsiC + thinW*gThin;
  dRegScalar = sum(gRegP.*adaptData.dPsi_direct,2);
  gSubNode = zeros(nelem,8);
  for jj=1:mnum
      gSubNode = gSubNode + adaptData.dPsi_dNode{jj} .* gRegP(:,jj);
  end
  gRegNode = accumarray(elemNodes(:),gSubNode(:),[ny*nz*nx 1],@sum,0);
  dRegScalar = dRegScalar + Wne'*gRegNode;
  dRegScalar(~activeDesign)=0;
  % ---- E.3 FE assembly & MGCG solve
  sK = reshape(KE(:)*Ee',576*nelem,1);
  K1 = sparse(iK,jK,sK);
  K1 = Null'*K1*Null - (Null-speye(ndof));
  K = cell(nl,1); K{1}=K1;
  for l=1:nl-1, K{l+1}=Pu{l}'*(K{l}*Pu{l}); end
  Lfac = chol(K{nl},'lower'); Ufac=Lfac';
  [mgcgi,mgcgr,U] = mgcg(K{1},F,U,K,Lfac,Ufac,Pu,nl,1,cgtol,cgmax);
  % ---- E.4 Objective + CONSISTENT adaptive-subcell compliance gradient
  ce   = sum(U(edofMat)*KE.*U(edofMat),2);         % unit-E element energy
  comp = sum(Ee.*ce);                 % physical compliance of the full cantilever design domain
  compFull = comp; % full cantilever compliance
  penaltySum = matBalW*matMismatch + alValue + rw(1)*pInterface + ...
               rw(2)*pConnect + thinW*pThin;
  % exp meriting positive factor; analytic product-rule gradient.
  augFactor = exp(min(30,max(-20,penaltySum)));
   compAug = comp*augFactor; % augmented cantilever merit function used by OC gradient
  compAugFull = compAug;
  % Compliance derivative with respect to clustered analysis field.
  dCdxCMF = -ce .* adaptData.dE_direct;
  gElemNodeC = (-ce) .* adaptData.dE_dNode;        % nelem x 8
  gNodeC = accumarray(elemNodes(:),gElemNodeC(:),[ny*nz*nx 1],@sum,0);
  dCdxCMF = dCdxCMF + Wne' * gNodeC;
  dCdxCMF(passiveVoid | passiveSolid) = 0;
  % Product rule for Caug = C*(1 + wM*Pmat).
  dCdxCMF = augFactor*(dCdxCMF + comp*dRegScalar);
  % Exact pullback CMF -> MFC -> density filter -> design variable.
  dCdxMFC = clustered_material_pullback(dCdxCMF,Gmat,GmatS,cmfData);
  dCdxTilde = minimum_feature_pullback(dCdxMFC,xTilde,G,Gs,mfcData);
  dc = H' * (dCdxTilde ./ Hs);
  % ---- E.5 PHYSICAL VOLUME from the SAME adaptive subcells
  % Ve is the element solid occupancy obtained from the normalized first
  % Heaviside threshold.  Thus sum_j psi_j = rho and stiffness/volume share
  % exactly the same implicit-interface representation.
  volCurrent = mean(Ve(activeDesign));
  dVdxCMF = zeros(nelem,1);
  dVdxCMF(activeDesign) = adaptData.dV_direct(activeDesign)/nActive;
  gElemNodeV = zeros(nelem,8);
  gElemNodeV(activeDesign,:) = adaptData.dV_dNode(activeDesign,:)/nActive;
  gNodeV = accumarray(elemNodes(:),gElemNodeV(:),[ny*nz*nx 1],@sum,0);
  dVdxCMF = dVdxCMF + Wne' * gNodeV;
  dVdxCMF(passiveVoid | passiveSolid) = 0;
  dVdxMFC = clustered_material_pullback(dVdxCMF,Gmat,GmatS,cmfData);
  dVdxTilde = minimum_feature_pullback(dVdxMFC,xTilde,G,Gs,mfcData);
  dv = H' * (dVdxTilde ./ Hs);
  dv(passiveVoid | passiveSolid) = 0;
  f0val = compAugFull;
  fval  = volCurrent/volfrac - 1;
  % ---- E.6 OC update using the adaptive physical-volume constraint

  % FAST adaptive move limit: large early steps, smaller late steps
  if loop <= 40
      move = 0.050;
  elseif loop <= 90
      move = 0.030;
  elseif loop <= 150
      move = 0.015;
  else
      move = 0.008;
  end
  innerit = 0;
  % Classification-consistent OC corrector.  The adaptive level map is frozen
  % during each scalar lambda bisection, then refreshed at the accepted
  % candidate.  If the refreshed map changes the physical volume noticeably,
  % the bisection is repeated with the refreshed classification.
  fixedLevel = adaptData.level;
  maxClassPass = 3;
  ocVolTol = 2.5e-5;
  B = -dc ./ max(dv,1e-16);
  B(~activeDesign) = 0;
  B = max(B,1e-12);
  for classPass = 1:maxClassPass
      l1 = 0;
      l2 = 1e9;
      while (l2-l1)/(l1+l2+eps) > 1e-6
          innerit = innerit + 1;
          lmid = 0.5*(l1+l2);
          xnew = x .* sqrt(B/lmid);
          xnew = max(1e-9, max(x-move, min(1.0, min(x+move, xnew))));
          xnew(passiveVoid)  = 1e-9;
          xnew(passiveSolid) = 1.0;
          volTrial = adaptive_volume_from_design( ...
              xnew,H,Hs,G,Gs,mfcIter,Gmat,GmatS,cmfIter,Wne,elemNodes,fixedLevel,tgrid,betaUsed,activeDesign,passiveSolid,passiveVoid);
          if volTrial > volfrac
              l1 = lmid;
          else
              l2 = lmid;
          end
      end
      % Refresh the VFLS-guided classification at the accepted candidate.
      xTnew = (H*xnew)./Hs;
      [xMFCnew,~] = minimum_feature_field(xTnew,G,Gs,mfcIter);
      [xCMFnew,~] = clustered_material_field(xMFCnew,Gmat,GmatS,cmfIter);
      xCMFnew(passiveVoid)  = 0.0;
      xCMFnew(passiveSolid) = 1.0;
      refreshedLevel = adaptive_level_from_xtilde( ...
          xCMFnew,Wne,elemNodes,tgrid,adaptIter);
      volRefresh = adaptive_volume_from_design( ...
          xnew,H,Hs,G,Gs,mfcIter,Gmat,GmatS,cmfIter,Wne,elemNodes,refreshedLevel,tgrid,betaUsed,activeDesign,passiveSolid,passiveVoid);
      if isequal(refreshedLevel,fixedLevel) || abs(volRefresh-volfrac) <= ocVolTol
          %fixedLevel = refreshedLevel;
          break;
      end
      fixedLevel = refreshedLevel;
  end
  % No post-hoc density scaling is used: feasibility is enforced inside the
  % nonlinear OC bisection with refreshed adaptive classification.
  lam = lmid;  % OC multiplier retained for reporting
  change = full(max(abs(xnew-x)));
  x = xnew;



  % ---- E.7 Continuation
  penal = cont_update(penal,penalCn,loop);
  beta  = cont_update(beta, betaCn, loop);
  % ---- E.8 Log & Plot (multi-material iso)
  nSub1 = nnz(adaptData.level==1);
  nSub2 = nnz(adaptData.level==2);
  nSub4 = nnz(adaptData.level==4);
  meanSupport = mean(mfcData.support(activeDesign));
  rhoAct = sum(psiV(activeDesign,:),2);
  phaseAct = psiV(activeDesign,:)./max(rhoAct,1e-14);
  solidAct = rhoAct > 0.10;
  if any(solidAct)
      matPurity = mean(max(phaseAct(solidAct,:),[],2));
  else
      matPurity = 1.0;
  end
  mbErr = max(abs(matDev));
  fprintf([' It.:%4i  Cfull:%10.4e  Vphys:%7.4f/%7.4f  lam:%9.2e  ch:%6.3f' ...
           '  MG:%3i rr:%3.2e  p:%3.1f b:%3.1f  MFC:a=%.2f  CMF:a=%.2f k=%.1f' ...
           '  Pur=%.3f MB:w=%.2f e=%.2f IF=%.4g CN=%.4g Sub:[%d %d %d]\n'], ...
      loop,compFull,volCurrent,volfrac,lam,change,mgcgi,mgcgr,penalUsed,betaUsed, ...
      mfcIter.alpha,cmfIter.alpha,cmfIter.kappa,matPurity,matBalW,mbErr, ...
      nSub1,nSub2,nSub4);
  Stats(loop,:) = [f0val fval mgcgi mgcgr];
  % ===== Store convergence/adaptive history =====
  HistComp(loop)       = compFull;
  HistVolTotal(loop)   = volCurrent;
  HistAdaptCount(loop,:) = [nSub1 nSub2 nSub4];
  HistMFCAlpha(loop)     = mfcIter.alpha;
  HistMFCSupport(loop)   = meanSupport;
  HistCMFAlpha(loop)     = cmfIter.alpha;
  HistCMFKappa(loop)     = cmfIter.kappa;
  HistMatPurity(loop)    = matPurity;
  HistInterface(loop) = pInterface;
  HistConnection(loop) = pConnect;
  HistRegWeight(loop,:) = rw;
  for jj = 1:mnum
      HistMatVolFrac(loop,jj) = sum(Psub(activeDesign,jj))/nActive;
  end
  % ===== Save this iteration to CSV =====
  % Columns reproduce the command-window convergence line and also include
  % mean MFC support and the individual material volume fractions.
  csvFID = fopen(iterCSV,'a');
  if csvFID < 0
      error('Cannot append optimization-history CSV: %s',iterCSV);
  end
  fprintf(csvFID,['%d,%.15g,%.15g,%.15g,%.15g,%.15g,%d,%.15g,' ...
                  '%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,' ...
                  '%d,%d,%d,%.15g'], ...
      loop,compFull,volCurrent,volfrac,lam,change,mgcgi,mgcgr, ...
      penalUsed,betaUsed,mfcIter.alpha,cmfIter.alpha,cmfIter.kappa, ...
      matPurity,matBalW,mbErr,nSub1,nSub2,nSub4,meanSupport);
  for jj = 1:mnum
      fprintf(csvFID,',%.15g',HistMatVolFrac(loop,jj));
  end
  fprintf(csvFID,',%.15g,%.15g,%.15g,%.15g\n',pInterface,pConnect,rw(1),rw(2));
  fclose(csvFID);
  % ===== FAST convergence test =====
  if isfinite(prevCompForStop)
      relCompChange = abs(compFull-prevCompForStop)/max(abs(prevCompForStop),1e-12);
      if relCompChange <= relCompTol && abs(volCurrent-volfrac) <= volTolAbs && ...
              max(abs(matDev)) <= alCfg.targetRelTol
          stableCount = stableCount + 1;
      else
          stableCount = 0;
      end
  end
  HistALResidual(loop,:) = matDev;
  HistALPenalty(loop) = alCfg.rho;
  HistALMultipliers(loop,:) = alCfg.mu;
  if alReady && mod(loop,alCfg.updateEvery)==0
      % Multiplier ascent for signed equality residuals; no post-hoc
      % material relabeling or density correction.
      alCfg.mu = max(-alCfg.lambdaMax,min(alCfg.lambdaMax, ...
          alCfg.mu + alCfg.rho*alResidual));
      err = max(abs(alResidual));
      if err > 0.90*alCfg.previousError
          alCfg.stallUpdates = alCfg.stallUpdates+1;
      else
          alCfg.stallUpdates = 0;
      end
      if alCfg.stallUpdates>=3 && err>alCfg.targetRelTol
          alCfg.rho = min(alCfg.rhoMax,alCfg.rho*alCfg.rhoGrowth);
          alCfg.stallUpdates = 0;
      end
      alCfg.previousError=err;
      fprintf(' V3.4 AL | maxRelMatErr=%.4f | rho=%.3g | mu=%s\n', ...
          max(abs(matDev)),alCfg.rho,mat2str(alCfg.mu,3));
  end
  prevCompForStop = compFull;
  if loop >= minOptIter && stableCount >= stableNeed && ...
          max(abs(matDev)) <= alCfg.targetRelTol
      fprintf(' FAST-V2 STOP: stable for %d iterations (rel dC <= %.1e, |dV| <= %.1e).\n', ...
          stableNeed,relCompTol,volTolAbs);
      break;
  end
nm = size(psiV,2);
%% ===== MATERIAL COLORS ACCORDING TO TABLE 1 =====
% Material 1: Red     E/Es = 1.00
% Material 2: Orange  E/Es = 0.85
% Material 3: Blue    E/Es = 0.70
% Material 4: Purple  E/Es = 0.50
% Material 5: Green   E/Es = 0.25

%% ===== CASE 1: MATERIALS 4, 1 =====
% base = [0.50 0.00 0.65;   % Material 4: Purple
%         0.85 0.00 0.00];  % Material 1: Red

% %% ===== CASE 2: MATERIALS 5, 3, 1 =====
base = [0.00 0.65 0.00;   % Material 5: Green
        0.00 0.35 1.00;   % Material 3: Blue
        0.85 0.00 0.00];  % Material 1: Red

%% ===== CASE 3: MATERIALS 5, 4, 3, 2, 1 =====
% base = [0.00 0.65 0.00;   % Material 5: Green
%         0.50 0.00 0.65;   % Material 4: Purple
%         0.00 0.35 1.00;   % Material 3: Blue
%         1.00 0.50 0.00;   % Material 2: Orange
%         0.85 0.00 0.00];  % Material 1: Red
%%
if nm <= size(base,1)
    colors = base(1:nm,:);
else
    colors = lines(nm);
end
% ============================================================
% HONEST MATERIAL VISUALIZATION -- no artificial smoothing or dilation
% Element occupancy is decided by the original optimization phase fractions.
% Thresholds are explicit and identical for every material and iteration.
% The hard material IDs are EXCLUSIVELY assigned (one owner per element).
% ============================================================
rhoThr = 0.05;
ownThr = 0.10;
% NOTE: ownThr uses absolute psi_j, as in original code.
% Binary surfaces show the actual selected element-cell regions.
[psiMax,matID] = max(psiV,[],2);
rhoMat = sum(psiV,2);
showEl = activeDesign & (rhoMat > rhoThr) & (psiMax > ownThr);
matID(~showEl) = 0;
matID(passiveSolid) = nm+1; % passive solid shown gray when present
matID(passiveVoid) = 0;

% Diagnose connectivity from the 3-D voxel mask, NOT from rendered triangles.
if mod(loop,10)==0 || loop==1
    diagTable = mmto_true_connectivity(matID,nely,nelz,nelx,nm,loop);
    diagFile = fullfile(figDir,'True_Connectivity_History.csv');
    if loop==1 || ~exist(diagFile,'file')
        writetable(diagTable,diagFile);
    else
        writetable(diagTable,diagFile,'WriteMode','append','WriteVariableNames',false);
    end
    fprintf(' TRUE connectivity (6-face): global components=%d, detached voxels=%d\n', ...
        diagTable.Components(1),diagTable.DetachedCells(1));
    for jj=1:nm
        fprintf('   MAT%d: components=%d, detached=%d, cells=%d, contacts(other MAT)=%d\n', ...
            jj,diagTable.Components(jj+1),diagTable.DetachedCells(jj+1), ...
            diagTable.OccupiedCells(jj+1),diagTable.CrossMaterialFaces(jj+1));
    end
end

if mod(loop,renderCfg.liveEvery)==0 || loop==1
    figure(1); clf;
    set(gcf,'Color','w','Renderer','opengl','Position',[80 80 1450 680]);
    ax=axes; hold(ax,'on');
    mmto_draw_true_voxels(ax,matID,nely,nelz,nelx,nm,colors);
    scatter3(ax,supportXYZ(:,1),supportXYZ(:,2),supportXYZ(:,3), ...
        18,[0.75 0.05 0.05],'s','filled','HandleVisibility','off');
    quiver3(ax,loadPlotX,loadPlotY,loadPlotZ, ...
        zeros(size(loadPlotX)),zeros(size(loadPlotX)),loadPlotW,0, ...
        'Color',[0.10 0.35 0.95],'LineWidth',0.75,'HandleVisibility','off');
    mmto_true_axes(ax,nelx,nely,nelz);
    title(ax,sprintf('TRUE thresholded material domains | Iter %d | rhoThr=%.2f, ownThr=%.2f', ...
        loop,rhoThr,ownThr),'FontName','Times New Roman','FontSize',13);
    hleg=gobjects(nm,1);
    for j=1:nm
        hleg(j)=plot3(ax,NaN,NaN,NaN,'s','MarkerSize',9, ...
            'MarkerFaceColor',colors(j,:),'MarkerEdgeColor','none', ...
            'DisplayName',sprintf('Material %d',j));
    end
    legend(ax,hleg,'Location','northeastoutside'); legend boxoff;
    drawnow limitrate;

    figure(2); clf;
    set(gcf,'Color','w','Renderer','opengl','Position',[60 80 1500 680]);
    ax=axes; hold(ax,'on');
    mmto_draw_true_voxels(ax,matID,nely,nelz,nelx,nm,colors);
    mmto_true_axes(ax,nelx,nely,nelz);
    title(ax,sprintf('MMTO thresholded voxel domains (no display smoothing), Iter %d',loop), ...
        'FontName','Times New Roman','FontSize',14);
    drawnow limitrate;
    if exist('v','var')
        try, writeVideo(v,getframe(gcf)); catch, end
    end
end
if mod(loop,renderCfg.saveEvery)==0
    exportgraphics(figure(1),fullfile(figDir,sprintf('Cantilever3D_TRUE_Iter_%04d.png',loop)), ...
        'Resolution',300);
    exportgraphics(figure(2),fullfile(figDir,sprintf('Cantilever3D_TRUE_Publication_Iter_%04d.png',loop)), ...
        'Resolution',300);
end

% Separated material domains use the EXACT SAME hard ownership map as Figure 1.
if mod(loop,10)==0
    figure(10); clf;
    set(gcf,'Color','w','Renderer','opengl','Position',[40 100 1700 500]);
    for j=1:nm
        ax=subplot(1,nm,j); hold(ax,'on');
        mmto_draw_true_single(ax,matID==j,nely,nelz,nelx,colors(j,:));
        mmto_true_axes(ax,nelx,nely,nelz);
        title(ax,sprintf('MAT %d | physical %.4f | owned %d',j, ...
            sum(psiV(activeDesign,j))/nActive,nnz(matID==j)), ...
            'FontName','Times New Roman','FontSize',11);
    end
    sgtitle(sprintf('TRUE separated material voxels | Iter %d | no smoothing',loop));
    drawnow;
    exportgraphics(gcf,fullfile(figDir,sprintf('Separated_TRUE_Iter_%04d.png',loop)), ...
        'Resolution',300);
end

% ============================================================
% Figure 20: Convergence graph
% Left Y-axis  : Compliance
% Right Y-axis : Total design volume fraction + material volume fractions
% ============================================================
if mod(loop,1)==0 || loop==1
    it = 1:loop;
    figure(20); clf;
    set(gcf,'Color','w','Position',[120 120 900 560]);
    yyaxis left
    plot(it, HistComp(it), '-ok', ...
        'LineWidth',2, ...
        'MarkerSize',3, ...
        'MarkerFaceColor','w');
    ylabel('Compliance','FontName','Times New Roman','FontSize',14);
    xlabel('Iteration number','FontName','Times New Roman','FontSize',14);
    grid on; box on;
    yyaxis right
    hold on;
    plot(it, HistVolTotal(it), '-s', ...
        'Color',[0.25 0.25 0.25], ...
        'LineWidth',2, ...
        'MarkerSize',3, ...
        'MarkerFaceColor','w');
    for jj = 1:mnum
        plot(it, HistMatVolFrac(it,jj), ':', ...
            'Color',colors(jj,:), ...
            'LineWidth',2);
    end
    ylabel('Volume fraction','FontName','Times New Roman','FontSize',16);
    ax = gca;
    ax.FontName = 'Times New Roman';
    ax.FontSize = 16;
    ax.LineWidth = 1.0;
    xlim([1 max(loop,2)]);
    %yyaxis right
    % volume-axis autoscaling
    yyaxis right
    volTotalVec = HistVolTotal(it);
    matVolFrac  = HistMatVolFrac(it,:);
    volAll = [volTotalVec(:); matVolFrac(:)];
    volAll = volAll(isfinite(volAll));
    if isempty(volAll)
        maxVolPlot = max(0.35,1.15*volfrac);
    else
        maxVolPlot = max(volAll);
    end
    ylim([0 max(0.35,1.15*maxVolPlot)]);
    yyaxis left
    minC = min(HistComp(it),[],'omitnan');
    maxC = max(HistComp(it),[],'omitnan');
    if isfinite(minC) && isfinite(maxC) && maxC > minC
        ylim([0.95*minC 1.05*maxC]);
    end
    legText = [{'Compliance','Design Volume Fraction'}, ...
               arrayfun(@(jj)sprintf('Material Volume %d',jj),1:mnum,'UniformOutput',false)];
    legend(legText,'Location','northeast','FontName','Times New Roman','FontSize',16);
    title(sprintf('Convergence history (Iter %d)',loop), ...
        'FontName','Times New Roman','FontSize',16,'FontWeight','bold');
    drawnow limitrate;
end

end % ===== end while loop =====
if exist('v','var')
    try
        close(v);
    catch
    end
end
%% ===== (F) FINAL SYSTEM + FAIR LINEAR-SOLVER BENCHMARK =====
% The optimization itself uses MGCG.  For the paper/reviewer comparison,
% all solvers below are tested on EXACTLY the same final K and F, with:
%   * zero initial guess,
%   * identical tolerance and maximum iteration count,
%   * setup and solve times reported separately.
%
% Baselines:
%   1) CG (no preconditioner)
%   2) PCG-Jacobi
%   3) PCG-ICHOL (non-multigrid preconditioner)
%   4) MGCG (the proposed/used solver)
%   5) Sparse direct Cholesky (optional; OFF by default for 3-D memory safety)
fprintf('\n============================================================\n');
fprintf(' FINAL LINEAR-SOLVER BENCHMARK\n');
fprintf(' Same K, F, x0=0, tol=%.1e, maxit=%d\n',cgtol,cgmax);
fprintf('============================================================\n');
% Rebuild the SAME adaptive-subcell stiffness AND volume representation at the final design.
xTildeFinal = (H*x)./Hs;
if exist('mfcIterLast','var')
    mfcFinalCfg = mfcIterLast;
else
    mfcFinalCfg = mfc_cfg_for_iter(mfcCfg,max(loop,1));
end
[xMFCFinal,mfcFinal] = minimum_feature_field(xTildeFinal,G,Gs,mfcFinalCfg);
if exist('cmfIterLast','var')
    cmfFinalCfg = cmfIterLast;
else
    cmfFinalCfg = cmf_cfg_for_iter(cmfCfg,max(loop,1));
end
[xCMFFinal] = clustered_material_field(xMFCFinal,Gmat,GmatS,cmfFinalCfg);
xCMFFinal(passiveVoid)  = 0.0;
xCMFFinal(passiveSolid) = 1.0;
xPhys = xCMFFinal;
[hFinal,dhFinal] = heaviside_smooth(xCMFFinal,tgrid,betaUsed);
[psiVFinal,~] = MatIntModel(hFinal,dhFinal,1);
[psiEFinal,~] = MatIntModel(hFinal,dhFinal,penalUsed);  % output compatibility
if exist('adaptIterLast','var')
    adaptFinalCfg = adaptIterLast;
else
    adaptFinalCfg = adaptive_cfg_for_iter(adaptCfg,max(loop,1));
end
[EeFinal,VeFinal,adaptFinal] = vfls_adaptive_subcell_response( ...
    xCMFFinal,elemNodes,Wne,tgrid,betaUsed,penalUsed,Evec,Emin,adaptFinalCfg);
AdaptiveLevelFinal = adaptFinal.level;
sKfinal = reshape(KE(:)*EeFinal',576*nelem,1);
Kfinal  = sparse(iK,jK,sKfinal);
Kfinal  = Null'*Kfinal*Null - (Null-speye(ndof));
benchOpts = struct();
benchOpts.caseName       = sprintf('Cantilever3D_%dx%dx%d_V%.3f',nelx,nely,nelz,volfrac);
benchOpts.runDirect      = false;   % set true only if enough RAM is available
benchOpts.directMaxDOF   = 50000;   % safety threshold for direct factorization
benchOpts.icholDropTol   = 1e-3;
benchOpts.icholDiagComp  = 1e-3;
benchOpts.verbose        = true;
SolverBenchmark = benchmark_spd_solvers( ...
    Kfinal,F,Pu,nl,cgtol,cgmax,benchOpts);
benchCSV = fullfile(figDir,'Solver_Benchmark.csv');
writetable(SolverBenchmark,benchCSV);
fprintf('\nSolver benchmark saved to:\n  %s\n',benchCSV);
% Publication-ready solver comparison figure.
try
    solverNames = categorical(SolverBenchmark.Solver);
    figure(30); clf;
    set(gcf,'Color','w','Position',[150 120 980 520]);
    yyaxis left
    bar(solverNames,SolverBenchmark.TotalTime_s);
    ylabel('Total time (s)','FontName','Times New Roman','FontSize',13);
    yyaxis right
    plot(solverNames,SolverBenchmark.Iterations,'-o','LineWidth',1.8,'MarkerSize',6);
    ylabel('Iterations','FontName','Times New Roman','FontSize',13);
    grid on; box on;
    title('Linear-solver comparison on the final FE system', ...
        'FontName','Times New Roman','FontSize',14,'FontWeight','bold');
    set(gca,'FontName','Times New Roman','FontSize',12);
    exportgraphics(gcf,fullfile(figDir,'Solver_Benchmark.png'),'Resolution',500);
catch MEplot
    %warning('Solver benchmark plot skipped: %s',MEplot.message);
end
% Final adaptive physical-volume check (NO mass/resource constraint).
finalVolume = mean(VeFinal(activeDesign));
fprintf('\nFINAL VOLUME = %.8f / %.8f\n',finalVolume,volfrac);
fprintf('NO MASS/RESOURCE CONSTRAINT IS USED; volfrac is evaluated on the full active 3-D Cantilever design domain.\n');
fprintf('FINAL MFC alpha = %.3f, mean support = %.4f, support radius = %.2f cells\n', ...
    mfcFinalCfg.alpha,mean(mfcFinal.support(activeDesign)),mfcCfg.radius);
% FINAL diagnostic uses the REBUILT FINAL physical material field.
[psiMaxFinal,matIDFinal] = max(psiVFinal,[],2);
rhoFinal = sum(psiVFinal,2);
matIDFinal(~(activeDesign & rhoFinal>0.05 & psiMaxFinal>0.10))=0;
matIDFinal(passiveSolid)=mnum+1;
matIDFinal(passiveVoid)=0;
FinalConnectivity = mmto_true_connectivity(matIDFinal,nely,nelz,nelx,mnum,loop);
writetable(FinalConnectivity,fullfile(figDir,'FINAL_True_Connectivity.csv'));
fprintf('FINAL true 6-face structural components: %d (detached cells: %d)\n', ...
    FinalConnectivity.Components(1),FinalConnectivity.DetachedCells(1));
finalMatFrac = sum(adaptFinal.Psi(activeDesign,:),1)/nActive;
fprintf('FINAL adaptive-subcell material fractions = %s (sum=%.6f, V=%.6f)\n',mat2str(finalMatFrac,5),sum(finalMatFrac),finalVolume);
fprintf('TARGET material fractions = %s\n',mat2str(matBalCfg.targetFrac,5));
finalRelativeMaterialError=(finalMatFrac-matBalCfg.targetFrac)./matBalCfg.targetFrac;
fprintf('FINAL material relative errors = %s (tol %.3f)\n', ...
   mat2str(finalRelativeMaterialError,5),alCfg.targetRelTol);
if max(abs(finalRelativeMaterialError)) > alCfg.targetRelTol
   warning('V3.4 material target tolerance NOT met: max relative error %.3g. Increase iterations / tune AL.', ...
       max(abs(finalRelativeMaterialError)));
end
matCheck = table((1:mnum)',matBalCfg.targetFrac(:),finalMatFrac(:), ...
    finalRelativeMaterialError(:),'VariableNames', ...
    {'Material','Target','Actual','RelativeError'});
writetable(matCheck,fullfile(figDir,'FINAL_Material_Constraints.csv'));
alHist = table((1:loop)',HistALPenalty(1:loop), ...
    max(abs(HistALResidual(1:loop,:)),[],2), ...
    'VariableNames',{'Iteration','ALPenalty','MaxRelativeMaterialError'});
for jj=1:mnum
    alHist.(sprintf('Material%dResidual',jj))=HistALResidual(1:loop,jj);
end
writetable(alHist,fullfile(figDir,'Material_Constraint_History.csv'));
% Final adaptive-subcell statistics
subdivList = [1;2;4];
elemCount  = [nnz(AdaptiveLevelFinal==1); ...
              nnz(AdaptiveLevelFinal==2); ...
              nnz(AdaptiveLevelFinal==4)];
elemFrac   = elemCount/nelem;
sampleCount = elemCount .* (subdivList.^3);
AdaptiveSubcellStats = table(subdivList,elemCount,elemFrac,sampleCount, ...
    'VariableNames',{'SubdivPerAxis','Elements','ElementFraction','ConstitutiveSamples'});
fprintf('\nFINAL VFLS-GUIDED ADAPTIVE SUBCELL DISTRIBUTION\n');
disp(AdaptiveSubcellStats);
writetable(AdaptiveSubcellStats,fullfile(figDir,'Adaptive_Subcell_Stats.csv'));
% Save fields consistent with the final design used in the benchmark.
psiV = psiVFinal;
psiE = psiEFinal;
%% ===== Save final xPhys result to MAT file =====
if ~exist(figDir,'dir')
    mkdir(figDir);
end
save(fullfile(figDir, filename), ...
     'xPhys','x','psiV','psiE','EeFinal','VeFinal','AdaptiveLevelFinal', ...
     'alCfg','HistALResidual','HistALPenalty','HistALMultipliers','finalRelativeMaterialError','HistComp','HistVolTotal','HistMatVolFrac','HistAdaptCount','HistMFCAlpha','HistMFCSupport','HistCMFAlpha','HistCMFKappa','HistMatPurity','HistInterface','HistConnection','HistRegWeight','regCfg', ...
     'AdaptiveSubcellStats','SolverBenchmark','FinalConnectivity','matIDFinal','A2','degree2', ...
     'nelx','nely','nelz','volfrac','Evec','adaptCfg','adaptFinalCfg','mfcCfg','mfcFinalCfg','cmfCfg','cmfFinalCfg','matBalCfg','xTildeFinal','xMFCFinal','xCMFFinal', ...
     'penalUsed','betaUsed','minOptIter','maxOptIter','changeTol','volTolAbs', ...
     'passiveS','passiveV','passiveSolid','passiveVoid','activeDesign','fixedNodes','loadNodes','supportXYZ','fixeddofs', ...
     'loadVals','loadPlotX','loadPlotY','loadPlotZ','loadPlotW','-v7.3');
fprintf('\nSaved final xPhys and MMTO results to: %s\n', ...
        fullfile(figDir, filename));
fprintf('Saved optimization history CSV to: %s\n',iterCSV);
end  % ===== end main function =====

function mmto_true_axes(ax,nelx,nely,nelz)
    view(ax,[138 20]);
    axis(ax,[0.5 nelx+0.5 0.5 nely+0.5 0.5 nelz+0.5]);
    axis(ax,'equal'); axis(ax,'manual'); axis(ax,'off');
    daspect(ax,[1 1 1]);
    lighting(ax,'gouraud'); camlight(ax,'headlight'); camlight(ax,'right');
end

function mmto_draw_true_voxels(ax,matID,nely,nelz,nelx,nm,colors)
    for j=1:nm
        mmto_draw_true_single(ax,matID==j,nely,nelz,nelx,colors(j,:));
    end
    if any(matID==nm+1)
        mmto_draw_true_single(ax,matID==nm+1,nely,nelz,nelx,[0.48 0.48 0.48]);
    end
end

function mmto_draw_true_single(ax,mask,nely,nelz,nelx,color)
    if ~any(mask), return; end
    % Cells are located at integer centers 1,...,n. Zero padding places
    % the outer surface at the true element boundary (half integers).
    V=permute(reshape(logical(mask),[nely,nelz,nelx]),[3 1 2]);
    P=zeros(nelx+2,nely+2,nelz+2);
    P(2:end-1,2:end-1,2:end-1)=double(V);
    fv=isosurface(0:nelx+1,0:nely+1,0:nelz+1,permute(P,[2 1 3]),0.5);
    % MATLAB isosurface grid x->columns (dimension 2), y->rows (dimension 1).
    if isempty(fv.faces), return; end
    h=patch(ax,'Faces',fv.faces,'Vertices',fv.vertices, ...
        'FaceColor',color,'EdgeColor','none','FaceAlpha',1, ...
        'HandleVisibility','off');
    set(h,'FaceLighting','gouraud','AmbientStrength',0.45,'DiffuseStrength',0.55);
end

function T=mmto_true_connectivity(matID,ny,nz,nx,nm,iter)
    % Face-adjacent (6-connected) voxel components; no toolbox required.
    phase=reshape(matID,[ny,nz,nx]);

    count=zeros(nm+1,1); ncomp=count; detached=count; contacts=count;
    allSolid=phase>0;
    for k=1:nm+1
        if k==1, A=allSolid; else, A=(phase==k-1); end
        count(k)=nnz(A);
        [ncomp(k),largest]=mmto_true_component_count(A);
        detached(k)=count(k)-largest;
        if k>1
            j=k-1;
            contacts(k)=nnz(phase(1:end-1,:,:) ==j & phase(2:end,:,:) >0 & phase(2:end,:,:) ~=j) + ...
                nnz(phase(2:end,:,:) ==j & phase(1:end-1,:,:) >0 & phase(1:end-1,:,:) ~=j) + ...
                nnz(phase(:,1:end-1,:) ==j & phase(:,2:end,:) >0 & phase(:,2:end,:) ~=j) + ...
                nnz(phase(:,2:end,:) ==j & phase(:,1:end-1,:) >0 & phase(:,1:end-1,:) ~=j) + ...
                nnz(phase(:,:,1:end-1) ==j & phase(:,:,2:end) >0 & phase(:,:,2:end) ~=j) + ...
                nnz(phase(:,:,2:end) ==j & phase(:,:,1:end-1) >0 & phase(:,:,1:end-1) ~=j);
        end
    end
    % Row 1 is all materials combined; rows 2:end are MAT1,...,MATm.
    Iteration=repmat(iter,nm+1,1);
    Domain=string(["ALL";compose('MAT%d',(1:nm)')]);
    OccupiedCells=count; Components=ncomp; DetachedCells=detached;
    CrossMaterialFaces=contacts;
    T=table(Iteration,Domain,OccupiedCells,Components,DetachedCells,CrossMaterialFaces);
end

function [components,largest]=mmto_true_component_count(A)
    % Pure MATLAB breadth-first search, 6-face adjacency in (y,z,x) order.
    sz=size(A); if numel(sz)<3, sz(3)=1; end
    ny=sz(1); nz=sz(2); nx=sz(3);
    occupied=find(A(:)); components=0; largest=0;
    if isempty(occupied), return; end
    remaining=logical(A(:)); queue=zeros(numel(occupied),1);
    for seed=occupied(:)'
        if ~remaining(seed), continue; end
        components=components+1; head=1; tail=1; queue(1)=seed;
        remaining(seed)=false;
        while head<=tail
            id=queue(head); head=head+1;
            iy=mod(id-1,ny)+1;
            iz=mod(floor((id-1)/ny),nz)+1;
            ix=floor((id-1)/(ny*nz))+1;
            if iy>1,             cand=id-1;       if remaining(cand), tail=tail+1;queue(tail)=cand;remaining(cand)=false;end,end
            if iy<ny,            cand=id+1;       if remaining(cand), tail=tail+1;queue(tail)=cand;remaining(cand)=false;end,end
            if iz>1,             cand=id-ny;      if remaining(cand), tail=tail+1;queue(tail)=cand;remaining(cand)=false;end,end
            if iz<nz,            cand=id+ny;      if remaining(cand), tail=tail+1;queue(tail)=cand;remaining(cand)=false;end,end
            if ix>1,             cand=id-ny*nz;   if remaining(cand), tail=tail+1;queue(tail)=cand;remaining(cand)=false;end,end
            if ix<nx,            cand=id+ny*nz;   if remaining(cand), tail=tail+1;queue(tail)=cand;remaining(cand)=false;end,end
        end
        largest=max(largest,tail);
    end
end

%% ======== V3.3 differentiable optimization regularization ========
function [A,nFaces]=mmto_face_adjacency(ny,nz,nx)
% 6-face graph with no wrap across (y,z,x) index boundaries.
    id=reshape(1:ny*nz*nx,ny,nz,nx);
    u=[reshape(id(1:end-1,:,:),[],1); reshape(id(:,1:end-1,:),[],1); reshape(id(:,:,1:end-1),[],1)];
    v=[reshape(id(2:end,:,:),[],1); reshape(id(:,2:end,:),[],1); reshape(id(:,:,2:end),[],1)];
    nFaces=numel(u); n=ny*nz*nx;
    A=sparse([u;v],[v;u],1,n,n);
end
function rw=mmto_v31_reg_weights(cfg,iter)
    if ~cfg.enabled || iter<=cfg.startIter, a=0;
    elseif iter>=cfg.fullIter, a=1;
    else
        t=(iter-cfg.startIter)/max(1,cfg.fullIter-cfg.startIter);
        a=t*t*(3-2*t);
    end
    rw=a*[cfg.interfaceWeight,cfg.connectWeight];
end
function [pI,pC,gI,gC]=mmto_v31_penalties(P,A,degree,active,epsN)
% P: actual (not sharpened) phase occupancies, ne x material count.
% INTERFACE: graph Dirichlet energy on phase occupancies. This includes
% some solid/void boundary penalty (perimeter regularization).
% CONNECTIVITY SURROGATE: penalize local phase occupancy lacking same-phase
% neighbors, with neighborhood weight normalized by nearby solid amount.
% It is NOT a global graph connectivity guarantee.
    P=max(0,double(P)); [ne,nm]=size(P);
    active=logical(active(:)); D=sum(active); D=max(D,1);
    PA=A*P;
    deg=degree(:);
    pI=sum(sum(P.*(deg.*P-PA)))/(D*nm);
    gI=2*(deg.*P-PA)/(D*nm);
    gI(~active,:)=0;
    % Local void-safe co-occupancy normalization:
    % q_j=(A*psi_j)/(eps + A*sum_k psi_k).
    % derivative includes both numerator and shared solid denominator.
    rho=sum(P,2); nbr=A*rho;
    den=epsN+nbr; q=PA./den;
    delta=1-q;
    pC=sum(sum(P.*delta.^2))/(D*nm);
    gC=delta.^2/(D*nm);
    % dP/d(q_j)=-2*P_j*(1-q_j)/(D*nm)
    t=-2*P.*delta/(D*nm);
    % d q_j = (A dP_j)/den - (PA_j/den^2)*(A d rho)
    gC=gC + A'*(t./den);
    cross=sum(t.*PA,2)./(den.^2);
    common=A'*cross;
    gC=gC-repmat(common,1,nm);
    gC(~active,:)=0;
end
function mmto_v31_gradient_check(x,tgrid,beta,A,degree,active,cfg)
% Check the complete scalar->hierarchical phase->regularizer derivative.
    x=x(:); [h,dh]=heaviside_smooth(x,tgrid,beta);
    [P,dCells]=MatIntModel(h,dh,1);
    dP=zeros(size(P));
    for k=1:numel(dCells), dP=dP+dCells{k}; end
    [~,~,gI,gC]=mmto_v31_penalties(P,A,degree,active,cfg.epsNeighbor);
    g=sum((gI+gC).*dP,2);
    ids=find(active & x>1e-4 & x<1-1e-4);
    if isempty(ids), return; end
    ids=ids(unique(round(linspace(1,numel(ids),min(4,numel(ids))))));
    hStep=1e-6; err=0;
    for id=ids(:)'
        xx=x; xx(id)=xx(id)+hStep;
        [hh,dd]=heaviside_smooth(xx,tgrid,beta);
        PP=MatIntModel(hh,dd,1);
        [i1,c1]=mmto_v31_penalties(PP,A,degree,active,cfg.epsNeighbor);
        xx(id)=x(id)-hStep;
        [hh,dd]=heaviside_smooth(xx,tgrid,beta);
        PP=MatIntModel(hh,dd,1);
        [i2,c2]=mmto_v31_penalties(PP,A,degree,active,cfg.epsNeighbor);
        numerical=((i1+c1)-(i2+c2))/(2*hStep);
        err=max(err,abs(numerical-g(id))/max(1e-8,abs(numerical)+abs(g(id))));
    end
    fprintf(' V3.3 regularization gradient relative discrepancy = %.3e\n',err);
    if err>5e-3, error('V3.3 regularization gradient validation failed (%.3e).',err); end
end

%% ======== local helper functions ========
function x0 = inverse_normalized_heaviside(rhoTarget,eta,beta)
% Uniform design value whose normalized Heaviside occupancy equals rhoTarget.
    rhoTarget = min(1-1e-12,max(1e-12,double(rhoTarget)));
    a = tanh(beta*eta);
    b = tanh(beta*(1-eta));
    den = a+b;
    z = rhoTarget*den-a;
    z = min(1-1e-12,max(-1+1e-12,z));
    x0 = eta + atanh(z)/beta;
    x0 = min(1,max(0,x0));
end
function y = cont_update(v,cn,iter)
  y = min(v + v*(iter>=cn{3})*(mod(iter-cn{3},cn{4})==0)*cn{5}, cn{2});
end
function [h,dh] = heaviside_smooth(xT, tgrid, beta)
% Normalized threshold projection: h(0)=0 and h(1)=1 for every threshold.
% This keeps the void/material partition consistent with the volume measure.
  ne = numel(xT); m = numel(tgrid);
  h  = zeros(ne,m); dh = zeros(ne,m);
  xT = min(1,max(0,double(xT(:))));
  for j=1:m
    eta = tgrid(j);
    a = tanh(beta*eta);
    b = tanh(beta*(1-eta));
    den = max(a+b,1e-14);
    s = beta*(xT-eta);
    th = tanh(s);
    h(:,j)  = (a + th)/den;
    dh(:,j) = beta*(1-th.^2)/den;
    h(:,j) = min(1,max(0,h(:,j)));
  end
end
function [psi,dpsi] = MatIntModel(h,dh,pen)
% Exact product rule. dpsi{k}(:,j) = d psi_j / d x via threshold k.
% The prior routine inadvertently multiplied non-differentiated factors
% (and could expand arrays for k>j+1); that corrupted material gradients.
    [ne,m] = size(h);
    hp = h.^pen;
    dhp = pen*max(h,1e-14).^(pen-1).*dh;
    psi=zeros(ne,m);
    dpsi=cell(m,1);
    for k=1:m, dpsi{k}=zeros(ne,m); end
    for j=1:m
        if j<m
            fac=[hp(:,1:j), 1-hp(:,j+1)];
        else
            fac=hp;
        end
        psi(:,j)=prod(fac,2);
        for k=1:size(fac,2)
            others=ones(ne,1);
            for t=1:size(fac,2)
                if t~=k, others=others.*fac(:,t); end
            end
            sg=1;
            if j<m && k==j+1, sg=-1; end
            dpsi{k}(:,j)=sg*dhp(:,k).*others;
        end
    end
end

function [G,Gs] = build_scalar_filter_3d(nelx,nely,nelz,radius)
% Weighted neighborhood operator for minimum-feature support measurement.
% Same structured element ordering as the main density filter.
    ne = nelx*nely*nelz;
    rr = max(1,double(radius));
    span = ceil(rr)-1;
    maxN = ne*(2*span+1)^3;
    iG = zeros(maxN,1); jG = zeros(maxN,1); sG = zeros(maxN,1);
    c = 0;
    for ix = 1:nelx
      for iz = 1:nelz
        for iy = 1:nely
          e1 = (ix-1)*nely*nelz + (iz-1)*nely + iy;
          for jx = max(ix-span,1):min(ix+span,nelx)
            for jz = max(iz-span,1):min(iz+span,nelz)
              for jy = max(iy-span,1):min(iy+span,nely)
                d = sqrt((ix-jx)^2+(iy-jy)^2+(iz-jz)^2);
                w = max(0,rr-d);
                if w > 0
                  e2 = (jx-1)*nely*nelz + (jz-1)*nely + jy;
                  c = c+1; iG(c)=e1; jG(c)=e2; sG(c)=w;
                end
              end
            end
          end
        end
      end
    end
    G = sparse(iG(1:c),jG(1:c),sG(1:c),ne,ne);
    Gs = full(sum(G,2));
    Gs(Gs<=0) = 1;
end
function cfgI = mfc_cfg_for_iter(cfg,iter)
% Smooth continuation prevents the minimum-feature control from dictating
% the initial global load path. It becomes fully active only after the
% topology has substantially formed.
    cfgI = cfg;
    if ~(isfield(cfg,'enabled') && cfg.enabled)
        cfgI.alpha = 0;
        return;
    end
    if iter <= cfg.startIter
        a = 0;
    elseif iter >= cfg.fullIter
        a = cfg.maxAlpha;
    else
        a = cfg.maxAlpha*(iter-cfg.startIter)/max(1,cfg.fullIter-cfg.startIter);
    end
    cfgI.alpha = min(cfg.maxAlpha,max(0,a));
end
function [xMFC,data] = minimum_feature_field(xTilde,G,Gs,cfg)
% Neighborhood-supported density control.
% q is a wider local average. Unsupported thin/isolated density is smoothly
% attenuated through a normalized support projection. No hard deletion is
% performed, keeping derivatives continuous for OC.
    xTilde = min(1,max(0,double(xTilde(:))));
    q = (G*xTilde)./Gs;
    alpha = 0;
    if isfield(cfg,'alpha'), alpha = cfg.alpha; end
    if ~(isfield(cfg,'enabled') && cfg.enabled) || alpha <= 0
        support = ones(size(xTilde));
        ds = zeros(size(xTilde));
        retention = ones(size(xTilde));
        xMFC = xTilde;
    else
        [support,ds] = normalized_support_projection(q,cfg.etaNeighbor,cfg.betaNeighbor);
        retention = (1-alpha) + alpha*support;
        xMFC = xTilde .* retention;
    end
    xMFC = min(1,max(0,xMFC));
    data = struct('q',q,'support',support,'ds',ds, ...
                  'retention',retention,'alpha',alpha);
end
function gTilde = minimum_feature_pullback(gMFC,xTilde,G,Gs,data)
% Chain-rule transpose for xMFC = xTilde .* r(q), q=(G*xTilde)./Gs.
    gMFC = double(gMFC(:));
    xTilde = double(xTilde(:));
    if data.alpha <= 0
        gTilde = gMFC;
        return;
    end
    localTerm = data.retention .* gMFC;
    qTerm = (xTilde .* data.alpha .* data.ds .* gMFC) ./ Gs;
    gTilde = localTerm + G' * qTerm;
end
function [s,ds] = normalized_support_projection(q,eta,beta)
% Smooth normalized 0-1 support gate with exact endpoint consistency.
    q = min(1,max(0,double(q(:))));
    a = tanh(beta*eta);
    b = tanh(beta*(1-eta));
    den = max(a+b,1e-14);
    th = tanh(beta*(q-eta));
    s = (a+th)/den;
    ds = beta*(1-th.^2)/den;
    s = min(1,max(0,s));
end
function w = material_balance_weight(cfg,iter)
% Continuation for material-volume balancing; leaves early topology formation free.
    if ~(isfield(cfg,'enabled') && cfg.enabled) || iter <= cfg.startIter
        w = 0;
    elseif iter >= cfg.fullIter
        w = cfg.maxWeight;
    else
        w = cfg.maxWeight*(iter-cfg.startIter)/max(1,cfg.fullIter-cfg.startIter);
    end
end
function cfgI = cmf_cfg_for_iter(cfg,iter)
% Smooth continuation for clustered material-domain formation.
    cfgI = cfg;
    if ~(isfield(cfg,'enabled') && cfg.enabled)
        cfgI.alpha = 0; cfgI.kappa = cfg.kappaStart; return;
    end
    if iter <= cfg.startIter
        a = 0;
    elseif iter >= cfg.fullIter
        a = cfg.maxAlpha;
    else
        a = cfg.maxAlpha*(iter-cfg.startIter)/max(1,cfg.fullIter-cfg.startIter);
    end
    if iter <= cfg.startIter
        k = cfg.kappaStart;
    elseif iter >= cfg.kappaFullIter
        k = cfg.kappaMax;
    else
        tt = (iter-cfg.startIter)/max(1,cfg.kappaFullIter-cfg.startIter);
        k = cfg.kappaStart + tt*(cfg.kappaMax-cfg.kappaStart);
    end
    cfgI.alpha = min(cfg.maxAlpha,max(0,a));
    cfgI.kappa = min(cfg.kappaMax,max(cfg.kappaStart,k));
end
function [xC,data] = clustered_material_field(xIn,Gm,Gms,cfg)
% Differentiable correlation-length-inspired material clustering map.
    xIn = min(1,max(0,double(xIn(:))));
    q = (Gm*xIn)./Gms;
    alpha = 0; if isfield(cfg,'alpha'), alpha = cfg.alpha; end
    kappa = cfg.kappaStart; if isfield(cfg,'kappa'), kappa = cfg.kappa; end
    centers = double(cfg.centers(:)');
    if ~(isfield(cfg,'enabled') && cfg.enabled) || alpha <= 0
        xC = xIn;
        data = struct('q',q,'target',q,'dtarget',ones(size(q)), ...
                      'alpha',0,'kappa',kappa,'centers',centers);
        return;
    end
    D = q - centers;
    Z = -kappa*(D.^2);
    Z = Z - max(Z,[],2);
    W = exp(Z); W = W./max(sum(W,2),1e-14);
    target = W*centers';
    c2mean = W*(centers.^2)';
    dtarget = 2*kappa*max(0,c2mean-target.^2);
    xC = (1-alpha)*xIn + alpha*target;
    xC = min(1,max(0,xC));
    data = struct('q',q,'target',target,'dtarget',dtarget, ...
                  'alpha',alpha,'kappa',kappa,'centers',centers);
end
function gIn = clustered_material_pullback(gC,Gm,Gms,data)
% Exact transpose chain rule for xC=(1-a)x + a*T(Gm*x/Gms).
    gC = double(gC(:));
    if data.alpha <= 0, gIn = gC; return; end
    gIn = (1-data.alpha)*gC + ...
          Gm' * ((data.alpha*data.dtarget.*gC)./Gms);
end
function v = material_display_field(psiV,j,gamma)
% Display-only soft winner field for crisp non-overlapping MAT colors.
    P = max(0,double(psiV));

    rho = sum(P,2);
    phase = zeros(size(P));
    ids = rho > 1e-14;
    phase(ids,:) = P(ids,:)./rho(ids);
    phase = phase.^max(1,double(gamma));
    den = sum(phase,2); den(den<1e-14)=1;
    phase = phase./den;
    v = rho.*phase(:,j);
end
function cfgI = adaptive_cfg_for_iter(cfg,iter)
% Conservative adaptive-resolution continuation to avoid subcell-scale
% branches before the global load path has stabilized.
    cfgI = cfg;
    if iter <= 40
        cfgI.midSubdiv  = 1;
        cfgI.highSubdiv = 2;
    elseif iter <= 90
        cfgI.midSubdiv  = 2;
        cfgI.highSubdiv = 2;
    else
        cfgI.midSubdiv  = 2;
        cfgI.highSubdiv = 4;
    end
end
function level = adaptive_level_from_xtilde(xTilde,Wne,elemNodes,tgrid,cfg)
% Classify parent Hex8 elements for 1^3 / 2^3 / 4^3 constitutive sampling
% from the VFLS implicit material interfaces Phi_k=xTilde-t_k.
    xTilde = min(1,max(0,double(xTilde(:))));
    ne = numel(xTilde);
    level = ones(ne,1);
    if ~(isfield(cfg,'enabled') && cfg.enabled), return; end
    xNode = Wne*xTilde;
    qN = xNode(elemNodes);
    qMin  = min(qN,[],2);
    qMax  = max(qN,[],2);
    qMean = mean(qN,2);
    qSpan = qMax-qMin;
    centerDist = inf(ne,1);
    crossed = false(ne,1);
    for kk = 1:numel(tgrid)
        tk = tgrid(kk);
        centerDist = min(centerDist,abs(qMean-tk));
        crossed = crossed | ((qMin <= tk) & (qMax >= tk));
    end
    midMask = crossed | (centerDist <= cfg.level2Band) | ...
              (qSpan >= cfg.spanLevel2);
    level(midMask) = cfg.midSubdiv;
    highMask = crossed & ((centerDist <= cfg.level4Band) | ...
                          (qSpan >= cfg.spanLevel4));
    level(highMask) = cfg.highSubdiv;
end
function [Ee,Ve,data] = vfls_adaptive_subcell_response( ...
    xTilde,elemNodes,Wne,tgrid,beta,penal,Evec,Emin,cfg)
% Analysis-consistent VFLS-guided adaptive subcell response.
%
% Ee : effective element modulus used in K.
% Ve : effective physical solid occupancy used in the volume constraint.
% Both are integrated over exactly the same 1^3/2^3/4^3 subcells.
%
% The physical solid indicator is the normalized first threshold h_1, and
% the hierarchical material fractions satisfy sum_j psi_j = h_1. Therefore
% the volume constraint is consistent with the actual material partition.
    xTilde = min(1,max(0,double(xTilde(:))));
    ne = numel(xTilde);
    xNode = Wne*xTilde;
    qN = xNode(elemNodes);                    % ne x 8
    level = adaptive_level_from_xtilde(xTilde,Wne,elemNodes,tgrid,cfg);
    Ee          = zeros(ne,1);
    Ve          = zeros(ne,1);
    dE_direct   = zeros(ne,1);
    dV_direct   = zeros(ne,1);
    dE_dNode    = zeros(ne,8);
    dV_dNode    = zeros(ne,8);
    m=numel(Evec);
    Psi=zeros(ne,m);
    dPsi_direct=zeros(ne,m);
    dPsi_dNode=cell(m,1);
    for jj=1:m, dPsi_dNode{jj}=zeros(ne,8); end
    levelsToRun = unique(level(:))';
    for q = levelsToRun
        ids = find(level==q);
        if isempty(ids), continue; end
        if q <= 1
            [Etmp,dEtmp] = scalar_material_E_dE( ...
                xTilde(ids),tgrid,beta,penal,Evec,Emin);
            [Rtmp,dRtmp] = scalar_solid_rho_drho( ...
                xTilde(ids),tgrid(1),beta);
            Ee(ids)        = Etmp;
            Ve(ids)        = Rtmp;
            dE_direct(ids) = dEtmp;
            dV_direct(ids) = dRtmp;
            [hh,dd]=heaviside_smooth(xTilde(ids),tgrid,beta);
            [Ptmp,dCells]=MatIntModel(hh,dd,1);
            Psi(ids,:)=Ptmp;
            for kk=1:m, dPsi_direct(ids,:)=dPsi_direct(ids,:)+dCells{kk}; end
        else
            Ns = subcell_shape_matrix(q);     % ns x 8
            xs = qN(ids,:)*Ns.';              % nGroup x ns
            [Es,dEs] = scalar_material_E_dE( ...
                xs(:),tgrid,beta,penal,Evec,Emin);
            [Rs,dRs] = scalar_solid_rho_drho( ...
                xs(:),tgrid(1),beta);
            Es   = reshape(Es,size(xs));
            dEs  = reshape(dEs,size(xs));
            Rs   = reshape(Rs,size(xs));
            dRs  = reshape(dRs,size(xs));
            nsmp = size(Ns,1);
            Ee(ids) = mean(Es,2);
            Ve(ids) = mean(Rs,2);
            dE_dNode(ids,:) = (dEs/nsmp)*Ns;
            dV_dNode(ids,:) = (dRs/nsmp)*Ns;
            [hh,dd]=heaviside_smooth(xs(:),tgrid,beta);
            [Ps,dCells]=MatIntModel(hh,dd,1);
            dPs=zeros(size(Ps));
            for kk=1:m, dPs=dPs+dCells{kk}; end
            for jj=1:m
                Pj=reshape(Ps(:,jj),size(xs));
                Gj=reshape(dPs(:,jj),size(xs));
                Psi(ids,jj)=mean(Pj,2);
                dPsi_dNode{jj}(ids,:)=(Gj/nsmp)*Ns;
            end
        end
    end
    data = struct();
    data.level       = level;
    data.dE_direct   = dE_direct;
    data.dV_direct   = dV_direct;
    data.dE_dNode    = dE_dNode;
    data.dV_dNode    = dV_dNode;
    data.xNode       = xNode;
    data.Psi = Psi;
    data.dPsi_direct = dPsi_direct;
    data.dPsi_dNode = dPsi_dNode;
    assert(max(abs(sum(Psi,2)-Ve))<1e-7,'Subcell material volume mismatch');
end
function vol = adaptive_volume_from_design( ...
    xDesign,H,Hs,G,Gs,mfcCfg,Gmat,GmatS,cmfCfg,Wne,elemNodes,level,tgrid,beta,activeDesign,passiveSolid,passiveVoid)
% Physical volume inside OC using the SAME MFC + CMF mapping as FE analysis.
    xT = (H*xDesign(:))./Hs;
    xT = min(1,max(0,xT));
    [xM,~] = minimum_feature_field(xT,G,Gs,mfcCfg);
    [xA,~] = clustered_material_field(xM,Gmat,GmatS,cmfCfg);
    ne = numel(xA);
    if nargin < 16 || isempty(passiveSolid), passiveSolid = false(ne,1); end
    if nargin < 17 || isempty(passiveVoid),  passiveVoid  = false(ne,1); end
    xA(passiveSolid) = 1.0;
    xA(passiveVoid)  = 0.0;
    Ve = zeros(ne,1);
    qNeeded = unique(level(:))';
    needNodes = any(qNeeded > 1);
    if needNodes
        xNode = Wne*xA;
        qN = xNode(elemNodes);
    end
    for q = qNeeded
        ids = find(level==q);
        if isempty(ids), continue; end
        if q <= 1
            Ve(ids) = scalar_solid_rho_drho(xA(ids),tgrid(1),beta);
        else
            Ns = subcell_shape_matrix(q);
            xs = qN(ids,:)*Ns.';
            Rs = scalar_solid_rho_drho(xs(:),tgrid(1),beta);
            Rs = reshape(Rs,size(xs));
            Ve(ids) = mean(Rs,2);
        end
    end
    vol = mean(Ve(activeDesign));
end
function [rho,drho] = scalar_solid_rho_drho(xv,eta,beta)
% Normalized smooth solid indicator with exact endpoint consistency.
    xv = min(1,max(0,double(xv(:))));
    a = tanh(beta*eta);
    b = tanh(beta*(1-eta));
    den = max(a+b,1e-14);
    th = tanh(beta*(xv-eta));
    rho  = (a + th)/den;
    drho = beta*(1-th.^2)/den;
    rho = min(1,max(0,rho));
end
function [E,dEdx] = scalar_material_E_dE(xv,tgrid,beta,penal,Evec,Emin)
% Hierarchical partition-of-unity stiffness law using normalized threshold
% projections.  h_k(0)=0 and h_k(1)=1 for every material threshold.
    xv = min(1,max(0,double(xv(:))));
    Evec = Evec(:);
    m = numel(Evec);
    E    = Emin*ones(size(xv));
    dEdx = zeros(size(xv));
    P  = ones(size(xv));
    dP = zeros(size(xv));
    Eprev = Emin;
    for kk = 1:m
        eta = tgrid(kk);
        a = tanh(beta*eta);
        b = tanh(beta*(1-eta));
        den = max(a+b,1e-14);
        s  = beta*(xv-eta);
        th = tanh(s);
        h  = (a + th)/den;
        dh = beta*(1-th.^2)/den;
        h = min(1,max(0,h));
        hp = h.^penal;
        if penal == 1
            dhp = dh;
        else
            dhp = penal*(max(h,1e-15).^(penal-1)).*dh;
        end
        Pold  = P;
        dPold = dP;
        P  = Pold.*hp;
        dP = dPold.*hp + Pold.*dhp;
        dEmat = Evec(kk)-Eprev;
        E     = E    + dEmat*P;
        dEdx  = dEdx + dEmat*dP;
        Eprev = Evec(kk);
    end
end
function Ns = subcell_shape_matrix(q)
% SUBCELL_SHAPE_MATRIX
% Trilinear Hex8 shape functions evaluated at centers of q^3 uniform
% subcells in natural coordinates [-1,1]^3.
%
% Node ordering matches elemNodes in the main routine:
% 1(-,-,-), 2(+,-,-), 3(+,+,-), 4(-,+,-),
% 5(-,-,+), 6(+,-,+), 7(+,+,+), 8(-,+,+).
    q = max(1,round(q));
    a = -1 + ((1:q)-0.5)*(2/q);
    [r,s,t] = ndgrid(a,a,a);
    r = r(:); s = s(:); t = t(:);
    Ns = [ ...
        (1-r).*(1-s).*(1-t), ...
        (1+r).*(1-s).*(1-t), ...
        (1+r).*(1+s).*(1-t), ...
        (1-r).*(1+s).*(1-t), ...
        (1-r).*(1-s).*(1+t), ...
        (1+r).*(1-s).*(1+t), ...
        (1+r).*(1+s).*(1+t), ...
        (1-r).*(1+s).*(1+t)]/8;
end
function fv = vfls_iso_mesh(V,lvl,cfg)
% VFLS_ISO_MESH  VFLS-inspired smooth implicit surface for rendering.
%
% The MMTO design variables are NOT changed.  The continuous material field
% V is converted to an implicit function Phi = V-lvl.  The zero contour is
% regularized using the same two numerical ideas used by the velocity-field
% level-set formulation: Laplacian diffusion and signed-distance
% reinitialization with first-order 3-D upwind gradients.  The final zero
% level set is triangulated by ISOSURFACE and receives only a light Taubin
% mesh cleanup.  No SMOOTH3 is used.
    V = double(V);
    sz = size(V);
    if numel(sz) < 3 || any(sz < 2) || ~any(isfinite(V(:)))
        fv = struct('faces',zeros(0,3),'vertices',zeros(0,3));
        return;
    end
    % ---- 1) Refine the scalar grid for surface extraction
    q = max(1,round(cfg.isoUpsample));
    if q > 1
        gx = 1:sz(1); gy = 1:sz(2); gz = 1:sz(3);
        gxq = linspace(1,sz(1),(sz(1)-1)*q+1);
        gyq = linspace(1,sz(2),(sz(2)-1)*q+1);
        gzq = linspace(1,sz(3),(sz(3)-1)*q+1);
        [GX,GY,GZ]    = ndgrid(gx,gy,gz);
        [GXq,GYq,GZq] = ndgrid(gxq,gyq,gzq);
        Vq = interpn(GX,GY,GZ,V,GXq,GYq,GZq,'linear');
    else
        gxq = 1:sz(1); gyq = 1:sz(2); gzq = 1:sz(3);
        [GXq,GYq,GZq] = ndgrid(gxq,gyq,gzq);
        Vq = V;
    end
    if ~any(Vq(:) > lvl) || ~any(Vq(:) < lvl)
        fv = struct('faces',zeros(0,3),'vertices',zeros(0,3));
        return;
    end
    % ---- 2) Density/material field -> implicit level-set
    % Positive Phi = displayed material; negative Phi = void/outside.
    Phi = Vq - lvl;
    targetFrac = nnz(Phi >= 0)/numel(Phi);
    % ---- 3) VFLS-style diffusion + signed-distance reinitialization
    nPass = max(0,round(cfg.vflsDiffusePasses));
    for ip = 1:nPass
        [~,~,LapPhi,~] = vfls_phi_diff(Phi);
        Phi = Phi + cfg.vflsDiffusion*LapPhi;
        Phi = vfls_reinitialize(Phi,cfg.vflsReinitIter,cfg.vflsReinitDt);
    end
    % ---- 4) Preserve the original displayed volume fraction approximately.
    % Diffusion can move the zero contour slightly.  A constant level-set
    % shift restores the number of positive grid samples without changing
    % local boundary smoothness.
    if isfield(cfg,'preserveVolume') && cfg.preserveVolume && ...
            targetFrac > 0 && targetFrac < 1
        vals = sort(Phi(:),'ascend');
        n = numel(vals);
        kth = max(1,min(n,round((1-targetFrac)*n)));
        phiShift = vals(kth);
        Phi = Phi - phiShift;
    end
    % ---- 5) Extract the zero level set
    fv = isosurface(GXq,GYq,GZq,Phi,0);
    if isempty(fv.faces) || size(fv.vertices,1) < 4
        return;
    end
    % ---- 6) Light Taubin cleanup only (not the primary smoother)
    nTaubin = max(0,round(cfg.taubinIter));
    if nTaubin == 0
        return;
    end
    Fm = fv.faces;
    P  = fv.vertices;
    nv = size(P,1);
    I = [Fm(:,1); Fm(:,2); Fm(:,3); Fm(:,2); Fm(:,3); Fm(:,1)];
    J = [Fm(:,2); Fm(:,3); Fm(:,1); Fm(:,1); Fm(:,2); Fm(:,3)];
    A = spones(sparse(I,J,1,nv,nv));
    deg = full(sum(A,2));
    deg(deg==0) = 1;
    W = spdiags(1./deg,0,nv,nv)*A;
    freeze = false(nv,1);
    if isfield(cfg,'freezeBoundary') && cfg.freezeBoundary
        Edges = sort([Fm(:,[1 2]); Fm(:,[2 3]); Fm(:,[3 1])],2);
        [Eu,~,ic] = unique(Edges,'rows');
        cnt = accumarray(ic,1);
        boundaryEdges = Eu(cnt==1,:);
        if ~isempty(boundaryEdges)
            freeze(unique(boundaryEdges(:))) = true;
        end
    end
    P0 = P;
    lamS = cfg.lambda;
    muS  = cfg.mu;
    for itS = 1:nTaubin
        D = W*P - P;
        P1 = P + lamS*D;
        if any(freeze), P1(freeze,:) = P0(freeze,:); end
        D = W*P1 - P1;
        P2 = P1 + muS*D;
        if any(freeze), P2(freeze,:) = P0(freeze,:); end
        P = P2;
    end
    fv.vertices = P;
end
function Phi = vfls_reinitialize(Phi,reNum,dtRe)
% Reinitialize Phi toward a signed-distance function using the same
% first-order upwind structure as the 3-D VFLS MATLAB implementation.
    reNum = max(0,round(reNum));
    if reNum == 0, return; end
    [~,~,~,S] = vfls_phi_diff(Phi);
    smax = max(abs(S(:)));
    if ~isfinite(smax) || smax < 1e-12, return; end
    dtEff = min(abs(dtRe),0.45)/smax;
    for ir = 1:reNum
        [GradP,GradM,~,~] = vfls_phi_diff(Phi);
        Phi = Phi - dtEff .* (max(S,0).*GradP + min(S,0).*GradM - S);
    end
end
function [GradP,GradM,Grad2,S] = vfls_phi_diff(Phi)
% 3-D first-order upwind gradients, Laplacian and regularized sign(Phi).
% Mirrored ghost values impose a zero-normal-gradient boundary treatment.
    % x-like direction = dimension 2
    DxL = Phi - cat(2,2*Phi(:,1,:)-Phi(:,2,:),Phi(:,1:end-1,:));
    DxR = cat(2,Phi(:,2:end,:),2*Phi(:,end,:)-Phi(:,end-1,:)) - Phi;
    % y-like direction = dimension 1
    DyL = Phi - cat(1,2*Phi(1,:,:)-Phi(2,:,:),Phi(1:end-1,:,:));
    DyR = cat(1,Phi(2:end,:,:),2*Phi(end,:,:)-Phi(end-1,:,:)) - Phi;
    % z direction = dimension 3
    DzL = Phi - cat(3,2*Phi(:,:,1)-Phi(:,:,2),Phi(:,:,1:end-1));
    DzR = cat(3,Phi(:,:,2:end),2*Phi(:,:,end)-Phi(:,:,end-1)) - Phi;
    GradP = sqrt(max(DxL,0).^2 + min(DxR,0).^2 + ...
                 max(DyL,0).^2 + min(DyR,0).^2 + ...
                 max(DzL,0).^2 + min(DzR,0).^2);
    GradM = sqrt(min(DxL,0).^2 + max(DxR,0).^2 + ...
                 min(DyL,0).^2 + max(DyR,0).^2 + ...
                 min(DzL,0).^2 + max(DzR,0).^2);
    Grad2 = (DxR-DxL) + (DyR-DyL) + (DzR-DzL);
    gx = 0.5*(DxL+DxR);
    gy = 0.5*(DyL+DyR);
    gz = 0.5*(DzL+DzR);
    S = Phi ./ sqrt(Phi.^2 + gx.^2 + gy.^2 + gz.^2 + eps);
end
%% ===== Low-level functions =====
function KE = Ke3D(nu)
C = [2/9 1/18 1/24 1/36 1/48 5/72 1/3 1/6 1/12];
A11 = [-C(1) -C(3) -C(3) C(2) C(3) C(3); -C(3) -C(1) -C(3) -C(3) -C(4) -C(5);...
    -C(3) -C(3) -C(1) -C(3) -C(5) -C(4); C(2) -C(3) -C(3) -C(1) C(3) C(3);...
    C(3) -C(4) -C(5) C(3) -C(1) -C(3); C(3) -C(5) -C(4) C(3) -C(3) -C(1)];
B11 = [C(7) 0 0 0 -C(8) -C(8); 0 C(7) 0 C(8) 0 0; 0 0 C(7) C(8) 0 0;...
    0 C(8) C(8) C(7) 0 0; -C(8) 0 0 0 C(7) 0; -C(8) 0 0 0 0 C(7)];
A22 = [-C(1) -C(3) C(3) C(2) C(3) -C(3); -C(3) -C(1) C(3) -C(3) -C(4) C(5);...
    C(3) C(3) -C(1) C(3) C(5) -C(4); C(2) -C(3) C(3) -C(1) C(3) -C(3);...
    C(3) -C(4) C(5) C(3) -C(1) C(3); -C(3) C(5) -C(4) -C(3) C(3) -C(1)];
B22 = [C(7) 0 0 0 -C(8) C(8); 0 C(7) 0 C(8) 0 0; 0 0 C(7) -C(8) 0 0;...
    0 C(8) -C(8) C(7) 0 0; -C(8) 0 0 0 C(7) 0; C(8) 0 0 0 0 C(7)];
A12 = [C(6) C(3) C(5) -C(4) -C(3) -C(5); C(3) C(6) C(5) C(3) C(2) C(3);...
    -C(5) -C(5) C(4) -C(5) -C(3) -C(4); -C(4) C(3) C(5) C(6) -C(3) -C(5);...
    -C(3) C(2) C(3) -C(3) C(6) C(5); C(5) -C(3) -C(4) C(5) -C(5) C(4)];
B12 = [-C(9) 0 -C(9) 0 C(8) 0; 0 -C(9) -C(9) -C(8) 0 -C(8); C(9) C(9) -C(9) 0 C(8) 0;...
    0 -C(8) 0 -C(9) 0 C(9); C(8) 0 -C(8) 0 -C(9) -C(9); 0 C(8) 0 -C(9) C(9) -C(9)];
A13 = [-C(4) -C(5) -C(3) C(6) C(5) C(3); -C(5) -C(4) -C(3) -C(5) C(4) -C(5);...
    C(3) C(3) C(2) C(3) C(5) C(6); C(6) -C(5) -C(3) -C(4) C(5) C(3);...
    C(5) C(4) -C(5) C(5) -C(4) -C(3); -C(3) C(5) C(6) -C(3) C(3) C(2)];
B13 = [0 0 C(8) -C(9) -C(9) 0; 0 0 C(8) C(9) -C(9) C(9); -C(8) -C(8) 0 0 -C(9) -C(9);...
    -C(9) C(9) 0 0 0 -C(8); -C(9) -C(9) C(9) 0 0 C(8); 0 -C(9) -C(9) C(8) -C(8) 0];
A14 = [C(2) C(5) C(5) C(4) -C(5) -C(5); C(5) C(2) C(5) C(5) C(6) C(3);...
    C(5) C(5) C(2) C(5) C(3) C(6); C(4) C(5) C(5) C(2) -C(5) -C(5);...
    -C(5) C(6) C(3) -C(5) C(2) C(5); -C(5) C(3) C(6) -C(5) C(5) C(2)];
B14 = [-C(9) 0 0 -C(9) C(9) C(9); 0 -C(9) 0 -C(9) -C(9) 0; 0 0 -C(9) -C(9) 0 -C(9);...
    -C(9) -C(9) -C(9) -C(9) 0 0; C(9) -C(9) 0 0 -C(9) 0; C(9) 0 -C(9) 0 0 -C(9)];
A23 = [C(2) C(5) -C(5) C(4) -C(5) C(5); C(5) C(2) -C(5) C(5) C(6) -C(3);...
    -C(5) -C(5) C(2) -C(5) -C(3) C(6); C(4) C(5) -C(5) C(2) -C(5) C(5);...
    -C(5) C(6) -C(3) -C(5) C(2) -C(5); C(5) -C(3) C(6) C(5) -C(5) C(2)];
B23 = [-C(9) 0 0 -C(9) C(9) -C(9); 0 -C(9) 0 -C(9) -C(9) 0; 0 0 -C(9) C(9) 0 -C(9);...
    -C(9) -C(9) C(9) -C(9) 0 0; C(9) -C(9) 0 0 -C(9) 0; -C(9) 0 -C(9) 0 0 -C(9)];
KE = 1/(1+nu)/(2*nu-1)*([A11 A12 A13 A14; A12' A22 A23 A13'; A13' A23' A22 A12'; A14' A13 A12 A11] +...
    nu*[B11 B12 B13 B14; B12' B22 B23 B13'; B13' B23' B22 B12'; B14' B13 B12 B11]);
end
function [Pu] = prepcoarse(nex,ney,nez)
% Assemble state variable prolongation (structured trilinear)
maxnum = nex*ney*nez*20;
iP = zeros(maxnum,1); jP = zeros(maxnum,1); sP = zeros(maxnum,1);
nexc = nex/2; neyc = ney/2; nezc = nez/2;
vals = [1,0.5,0.25,0.125];
cc = 0;
for nx = 1:nexc+1
  for ny = 1:neyc+1
    for nz = 1:nezc+1
      col = (nx-1)*(neyc+1)+ny+(nz-1)*(neyc+1)*(nexc+1);
      nx1 = nx*2 - 1; ny1 = ny*2 - 1; nz1 = nz*2 - 1;
      for k = max(nx1-1,1):min(nx1+1,nex+1)
        for l = max(ny1-1,1):min(ny1+1,ney+1)
          for h = max(nz1-1,1):min(nz1+1,nez+1)
            row = (k-1)*(ney+1)+l+(h-1)*(nex+1)*(ney+1);
            ind = 1+((nx1-k)^2+(ny1-l)^2+(nz1-h)^2);
            cc=cc+1; iP(cc)=3*row-2; jP(cc)=3*col-2; sP(cc)=vals(ind);
            cc=cc+1; iP(cc)=3*row-1; jP(cc)=3*col-1; sP(cc)=vals(ind);
            cc=cc+1; iP(cc)=3*row;   jP(cc)=3*col;   sP(cc)=vals(ind);
          end
        end
      end
    end
  end
end
Pu = sparse(iP(1:cc),jP(1:cc),sP(1:cc));
end
function T = benchmark_spd_solvers(A,b,Pu,nl,tol,maxit,opts)
% BENCHMARK_SPD_SOLVERS  Fair comparison on the same SPD system A*u=b.
% Returns a table used by the final publication-style benchmark plot.
    if nargin < 7, opts = struct(); end
    if ~isfield(opts,'runDirect'),     opts.runDirect = false; end
    if ~isfield(opts,'directMaxDOF'),  opts.directMaxDOF = 50000; end
    if ~isfield(opts,'icholDropTol'),  opts.icholDropTol = 1e-3; end
    if ~isfield(opts,'icholDiagComp'), opts.icholDiagComp = 1e-3; end
    if ~isfield(opts,'verbose'),       opts.verbose = true; end
    n = size(A,1);
    b = b(:);
    x0 = zeros(n,1);
    nMax = 5;
    Solver = strings(nMax,1);
    SetupTime_s = nan(nMax,1);
    SolveTime_s = nan(nMax,1);
    TotalTime_s = nan(nMax,1);
    Iterations = nan(nMax,1);
    RelRes = nan(nMax,1);
    Flag = nan(nMax,1);
    krow = 0;
    % 1) Unpreconditioned CG
    krow = krow + 1;
    Solver(krow) = "CG";
    t0 = tic; SetupTime_s(krow) = toc(t0);
    t0 = tic;
    [~,Flag(krow),RelRes(krow),Iterations(krow)] = pcg(A,b,tol,maxit,[],[],x0);
    SolveTime_s(krow) = toc(t0);
    TotalTime_s(krow) = SetupTime_s(krow)+SolveTime_s(krow);
    % 2) Jacobi-PCG
    krow = krow + 1;
    Solver(krow) = "PCG-Jacobi";
    t0 = tic;
    d = full(diag(A));
    d(abs(d)<eps) = 1;
    Mj = spdiags(d,0,n,n);
    SetupTime_s(krow) = toc(t0);
    t0 = tic;
    [~,Flag(krow),RelRes(krow),Iterations(krow)] = pcg(A,b,tol,maxit,Mj,[],x0);
    SolveTime_s(krow) = toc(t0);
    TotalTime_s(krow) = SetupTime_s(krow)+SolveTime_s(krow);
    % 3) Incomplete-Cholesky PCG. Keep the benchmark alive if ICHOL fails.
    krow = krow + 1;
    Solver(krow) = "PCG-ICHOL";
    try
        t0 = tic;
        ichOpts = struct('type','ict','droptol',opts.icholDropTol, ...
                         'diagcomp',opts.icholDiagComp);
        Lich = ichol(A,ichOpts);
        SetupTime_s(krow) = toc(t0);
        t0 = tic;
        [~,Flag(krow),RelRes(krow),Iterations(krow)] = pcg(A,b,tol,maxit,Lich,Lich',x0);
        SolveTime_s(krow) = toc(t0);
        TotalTime_s(krow) = SetupTime_s(krow)+SolveTime_s(krow);
    catch ME
        Flag(krow) = -99;
        if opts.verbose
            %warning('ICHOL benchmark skipped: %s',ME.message);
        end
    end
    % 4) MGCG using the same geometric hierarchy as the optimizer.
    krow = krow + 1;
    Solver(krow) = "MGCG";
    try
        t0 = tic;
        Apre = cell(nl,1);
        Apre{1} = A;
        for lev = 1:nl-1
            Apre{lev+1} = Pu{lev}'*(Apre{lev}*Pu{lev});
        end
        Lc = chol(Apre{nl},'lower');
        Uc = Lc';
        SetupTime_s(krow) = toc(t0);
        t0 = tic;
        [Iterations(krow),RelRes(krow),~] = mgcg(A,b,x0,Apre,Lc,Uc,Pu,nl,1,tol,maxit);
        SolveTime_s(krow) = toc(t0);
        TotalTime_s(krow) = SetupTime_s(krow)+SolveTime_s(krow);
        Flag(krow) = double(RelRes(krow) >= tol);
    catch ME
        Flag(krow) = -99;
        if opts.verbose
            % warning('MGCG benchmark failed: %s',ME.message);
        end
    end
    % 5) Optional sparse direct Cholesky for small systems only.
    if opts.runDirect && n <= opts.directMaxDOF
        krow = krow + 1;
        Solver(krow) = "Direct-Cholesky";
        try
            t0 = tic;
            Ld = chol(A,'lower');
            SetupTime_s(krow) = toc(t0);
            t0 = tic;
            xd = Ld'\(Ld\b);
            SolveTime_s(krow) = toc(t0);
            TotalTime_s(krow) = SetupTime_s(krow)+SolveTime_s(krow);
            RelRes(krow) = norm(b-A*xd)/max(norm(b),eps);
            Iterations(krow) = 1;
            Flag(krow) = 0;
        catch ME
            Flag(krow) = -99;
            if opts.verbose
                % warning('Direct Cholesky benchmark failed: %s',ME.message);
            end
        end
    end
    Solver = Solver(1:krow);
    SetupTime_s = SetupTime_s(1:krow);
    SolveTime_s = SolveTime_s(1:krow);
    TotalTime_s = TotalTime_s(1:krow);
    Iterations = Iterations(1:krow);
    RelRes = RelRes(1:krow);
    Flag = Flag(1:krow);
    T = table(Solver,SetupTime_s,SolveTime_s,TotalTime_s,Iterations,RelRes,Flag);
    if opts.verbose
        fprintf('\nLinear solver benchmark:\n');
        disp(T);
    end
end
function [i,relres,u] = mgcg(A,b,u,Apre,Lfac,Ufac,Pu,nl,nswp,tol,maxiter)
r = b - A*u; res0 = norm(b);
omega = 0.6;
invD = cell(nl-1,1);
for l = 1:nl-1
  invD{l,1}= 1./spdiags(Apre{l,1},0);
end
for i = 1:1e6
  z = VCycle(Apre,r,Lfac,Ufac,Pu,1,nl,invD,omega,nswp);
  rz = r'*z;
  if i==1, p = z; else, beta = rz/rzold; p = beta*p + z; end
  Ap = A*p; alpha = rz/(p'*Ap);
  u = u + alpha*p; r = r - alpha*Ap; rzold = rz;
  relres = norm(r)/max(res0,1e-16);
  if (relres < tol || i >= maxiter), break; end
end
end
function z = VCycle(A,r,Lfac,Ufac,Pu,l,nl,invD,omega,nswp)
z = 0*r;
if (l == nl)
  z = Ufac \ (Lfac \ r);
else
  z = smthdmpjac(z,A{l,1},r,invD{l,1},omega,nswp);
  d = r - A{l,1}*z;
  dh2 = Pu{l,1}'*d;
  vh2 = VCycle(A,dh2,Lfac,Ufac,Pu,l+1,nl,invD,omega,nswp);
  z = z + Pu{l,1}*vh2;
  z = smthdmpjac(z,A{l,1},r,invD{l,1},omega,nswp);
end
end
function u = smthdmpjac(u,A,b,invD,omega,nswp)
for i = 1:nswp
  r = b - A*u;
  u = u + omega*(invD.*r);
end
end
