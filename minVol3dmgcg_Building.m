function minVol3dmgcg_Building(nelx,nely,nelz,volfrac,rmin,ft,nl,cgtol,cgmax,filename,Evec,Dvec)
% Multi-material MGCG bridge TO: Mass minimization with compliance constraint
% - Single design field x \in [0,1]
% - SVMM projection -> hierarchical partition-of-unity interpolation
% - MG preconditioned CG solver


% Call example: Vol=0.20

% minVol3dmgcg_bridge_Thype_II(140,20,80,0.12,sqrt(3),1,3, 1e-6,1000, 'bridge_2mm.mat', [0.5 1.0], [0.4 1]);

% E = [0.5 0.7 1.0];
% D = [0.4 0.6 1];
% minVol3dmgcg_bridge_Thype_II(140,20,80,0.12,sqrt(3),1,3, 1e-6,1000, 'bridge_3mm.mat', [0.25 0.7 1.0], [0.2 0.6 1]);

% E = [0.25 0.50 0.70 0.85 1.00];
% D = [0.20 0.40 0.60 0.80 1.00];
% minVol3dmgcg_Building(60,60,60*6,0.2,sqrt(3),1,3, 1e-6,150, 'bridge_5mm.mat', [0.25 0.50 0.70 0.85 1.00], [0.20 0.40 0.60 0.80 1.00]);
tic;

%% ===== (A) MATERIAL / CONTINUATION =====
Emin = 1e-9;                 % void stiffness (normalized)
nu   = 0.30;
% continuation: [start, max, start_iter, period, step]
penalCn = {1,  5, 20, 5, 0.25};    % SIMP penalization on stiffness blending
betaCn  = {1, 60, 20, 5, 0.10};    % Heaviside sharpness on filtered x
penal = penalCn{1}; 
beta  = betaCn{1};

% sanity check materials
mnum = numel(Evec);
assert(mnum==numel(Dvec), 'Evec and Dvec must have the same length (number of materials).');

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

% Multigrid prolongation operators
Pu=cell(nl-1,1);
for l = 1:nl-1
    Pu{l,1} = prepcoarse(nelz/2^(l-1),nely/2^(l-1),nelx/2^(l-1));
end

%% ============================================================
% Z-Topos building slabs model
% yellow horizontal slabs = passive solid
% internal organic part   = design domain
%% ============================================================

U = zeros(ndof,1);

%% ===== Element numbering =====
elmNrs = int32(reshape(1:nelx*nely*nelz,nely,nelz,nelx));

passiveSolid3D = false(nely,nelz,nelx);

%% ===== Passive horizontal slabs =====
nSlab     = 7;                      % จำนวนแผ่น slab
slabThick = max(1,round(0.04*nelz));
zSlab     = round(linspace(2,nelz-1,nSlab));

for k = 1:nSlab
    z1 = max(1,zSlab(k)-floor(slabThick/2));
    z2 = min(nelz,zSlab(k)+ceil(slabThick/2));

    passiveSolid3D(:,z1:z2,:) = true;
end

%% ===== Design domain: central organic core only =====
[Y,Z,X] = ndgrid(1:nely,1:nelz,1:nelx);

xc = (nelx+1)/2;
yc = (nely+1)/2;

xn = (X-xc)/(0.42*nelx);
yn = (Y-yc)/(0.42*nely);
zn = Z/nelz;

% variable organic radius along height
rCore = 0.58 ...
      - 0.22*sin(3*pi*zn).^2 ...
      + 0.10*cos(6*pi*zn);

coreRegion = xn.^2 + yn.^2 <= rCore.^2;

% outside core = passive void, except slab
passiveVoid3D = ~coreRegion & ~passiveSolid3D;

%% ===== Large vertical openings inside core =====
nHole = 3;
for h = 1:nHole
    ang = 2*pi*(h-1)/nHole;
    hx = xc + 0.18*nelx*cos(ang);
    hy = yc + 0.18*nely*sin(ang);

    hole = ((X-hx)/(0.07*nelx)).^2 + ...
           ((Y-hy)/(0.11*nely)).^2 <= 1;

    passiveVoid3D = passiveVoid3D | (hole & coreRegion & ~passiveSolid3D);
end

%% ===== Convert passive regions to vector =====
passiveSolid = passiveSolid3D(:);
passiveVoid  = passiveVoid3D(:);

% %% ===== Load: distributed load on top slab =====
% topNodes = nodenrs(1:ny,nz,1:nx);
% F = sparse(3*topNodes(:),1,-1/numel(topNodes),ndof,1);
%% ===== Load: distributed load on every slab =====
F = sparse(ndof,1);

totalLoad = -1;   % แรงรวมทั้งหมด

for k = 1:nSlab

    zMid = zSlab(k);

    % nodes บนระดับ slab นั้น
    slabNodes = nodenrs(1:ny,zMid,1:nx);

    % กระจายแรงเท่ากันทุก node และทุก slab
    F = F + sparse( ...
        3*slabNodes(:), ...
        1, ...
        totalLoad/(nSlab*numel(slabNodes)), ...
        ndof,1);

end
%% ===== Supports: fixed bottom slab / base plane =====
baseNodes = nodenrs(1:ny,1,1:nx);

fix_x = 3*baseNodes(:)-2;
fix_y = 3*baseNodes(:)-1;
fix_z = 3*baseNodes(:);

fixeddofs = unique([fix_x; fix_y; fix_z]);

%% ===== Null space elimination =====
N = ones(ndof,1);
N(fixeddofs) = 0;
Null = spdiags(N,0,ndof,ndof);
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
%% ===== (D) INIT =====
x = volfrac*ones(nelem,1);
x(passiveVoid)  = 1e-10;
x(passiveSolid) = 1.0;
xPhys = x;

loop = 0; change = 1; lam = 1e10;
Stats = zeros(300,4);
% ===== History for convergence plot =====
maxLoop = 300;
HistComp      = nan(maxLoop,1);
HistMassTotal = nan(maxLoop,1);
HistMassFrac  = nan(maxLoop,mnum);

% SVMM thresholds
tgrid = (1:mnum)/(mnum+1);
% %% ===== GIF setup (before MAIN LOOP) =====
gifname = sprintf('Bridge_MM_%dx%dx%d.gif', nelx, nely, nelz);
if exist(gifname, 'file'), delete(gifname); end   % ลบไฟล์เก่าถ้ามี
gif_started = false;    
% เตรียม figure ล่วงหน้า
figure(2);
clf;
set(gcf, 'Color', 'w', 'Renderer', 'opengl', 'Position', [100 100 960 540]);
%% ===== VIDEO setup (before MAIN LOOP) =====
vname = sprintf('Bridge_MM_%dx%dx%d.mp4', nelx, nely, nelz);
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
%% ===== (E) MAIN LOOP =====
while (change > 1e-2 && loop < 250)
  loop = loop + 1;

  % ---- E.1 Filter + Heaviside projection
  xTilde = (H*x)./Hs;                          
  [h,dh] = heaviside_smooth(xTilde,tgrid,beta);   % ne×m

  % ---- E.2 Multi-material interpolation
  [psiV,dpsiV] = MatIntModel(h,dh,1);            % for mass/volume
  [psiE,dpsiE] = MatIntModel(h,dh,penal);        % for stiffness
  De = psiV*Dvec(:);                              % ne×1
  Ee = Emin + psiE*(Evec(:)-Emin);                % ne×1


  % ---- E.3 FE assembly & MGCG solve
  sK = reshape(KE(:)*Ee',576*nelem,1);
  K1 = sparse(iK,jK,sK);
  K1 = Null'*K1*Null - (Null-speye(ndof));
  K = cell(nl,1); K{1}=K1;
  for l=1:nl-1, K{l+1}=Pu{l}'*(K{l}*Pu{l}); end
  Lfac = chol(K{nl},'lower'); Ufac=Lfac';
  [mgcgi,mgcgr,U] = mgcg(K{1},F,U,K,Lfac,Ufac,Pu,nl,1,cgtol,cgmax);

  % ---- E.4 Objective & constraint
  ce   = sum(U(edofMat)*KE.*U(edofMat),2);       % ne×1
  comp = sum(Ee.*ce);
  M    = sum(De)/nelem;

  volCurrent = mean(xPhys(~passiveVoid));
  % % sensitivities dC/dx, dM/dx via chain rule and filter
        % ===== compliance sensitivity dC/dx =====
  dC = zeros(nelem,1);
    
  for j = 1:mnum
        dEe_dxT_j = dpsiE{j}*(Evec(:)-Emin);
        dC = dC + H * (-(dEe_dxT_j).*ce ./ Hs);
  end
    % ===== volume sensitivity dV/dx =====
    dv = ones(nelem,1)/sum(~passiveVoid);
    dv(passiveVoid)  = 0;
    dv(passiveSolid) = 0;
    
    % objective = compliance minimization
    f0val = comp;
    fval  = volCurrent/volfrac - 1;
    
    dc = dC;
    dv = ones(nelem,1)/sum(~passiveVoid);  % volume sensitivity
    dv(passiveVoid) = 0;
    dv(passiveSolid) = 0;
  % ---- E.5 Sensitivity filter (optional)
  if ft==1
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
  elseif ft==2
    dc(:) = H*(dc(:)./Hs);  
    dv(:) = H*(dv(:)./Hs);
  end

  % ---- E.6 OC update
%  l1 = 0; l2 = lam*2; move = 0.2; innerit = 0;
      % ---- E.6 OC update for volume constraint ----
    activeDesign = ~(passiveVoid | passiveSolid);
    nActive = sum(activeDesign);
    
    l1 = 0;
    l2 = 1e9;
    move = 0.01;
    innerit = 0;

  while (l2-l1)/(l1+l2+eps) > 1e-6
    innerit = innerit + 1;
    lmid = 0.5*(l1+l2);

    B = -dc ./ max(dv,1e-16);
    B(~activeDesign) = 0;
    B = max(B,1e-12);

    xnew = x .* sqrt(B/lmid);

    xnew = max(1e-9, max(x-move, min(1.0, min(x+move, xnew))));

    % enforce passive regions
    xnew(passiveVoid)  = 1e-9;
    xnew(passiveSolid) = 1.0;

    % volume check only active design domain
    volTrial = mean(xnew(activeDesign));

    if volTrial > volfrac
        l1 = lmid;
    else
        l2 = lmid;
    end
  end
        % hard projection to target volume if still high
  volTrial = mean(xnew(activeDesign));
  if volTrial > volfrac
        xnew(activeDesign) = xnew(activeDesign) * volfrac / max(volTrial,1e-12);
  end
    
  xnew(activeDesign) = min(1,max(1e-9,xnew(activeDesign)));
  xnew(passiveVoid)  = 1e-9;
  xnew(passiveSolid) = 1.0;

  change = full(max(abs(xnew-x)));  
  x = xnew;
  if ft==1, xPhys = x; else, xPhys = (H*x)./Hs; end

  % ---- E.7 Continuation
  penal = cont_update(penal,penalCn,loop);
  beta  = cont_update(beta, betaCn, loop);

  % ---- E.8 Log & Plot (multi-material iso)
  fprintf(' It.:%4i  Mass:%8.4f  Comp:%10.4e  lam:%9.2e  ch:%6.3f  it:%3i  rr:%3.2e  p:%3.1f  b:%3.1f\n',...
      loop,f0val,comp,lam,change,mgcgi,mgcgr,penal,beta);
  Stats(loop,:) = [f0val fval mgcgi mgcgr];

      % ===== Store convergence history =====
  HistComp(loop)      = comp;
  HistMassTotal(loop) = sum(De)/nelem;
  for jj = 1:mnum
      HistMassFrac(loop,jj) = sum(psiV(:,jj))/nelem;
  end

nm = size(psiV,2);

% % %case 1
%base = [0.60 0.00 0.80];  % Material 5: purple
% % %case 1
% base = [0.00 0.35 1.00;   % Material 4: blue
%         0.60 0.00 0.80];  % Material 5: purple
%case 2
% base = [0.00 0.65 0.00;   % Material 1: green
%         0.85 0.00 0.00;   % Material 3: red
%         0.60 0.00 0.80];  % Material 5: purple
% %case 3
base = [0.00 0.65 0.00;   % Material 1: green
        0.00 0.35 1.00;   % Material 2: blue
        0.85 0.00 0.00;   % Material 3: red
        1.00 0.50 0.00;   % Material 4: orange
        0.60 0.00 0.80];  % Material 5: purple

if nm <= size(base,1)
    colors = base(1:nm,:);
else
    colors = lines(nm); % fallback
end

% ============================================================
% Figure 1: Half / current computational domain, clean material plot
% ============================================================
if mod(loop,1)==0
    figure(1); clf; set(gcf,'Color','w');
    ax = axes; hold(ax,'on');

    drew_any = false;

    for j = 1:nm
        V = reshape(psiV(:,j), [nely, nelz, nelx]);   % y-z-x
        V = permute(V, [3 1 2]);              % x-y-z

        % smooth only the displayed binary material field
        V = smooth3(V,'box',5);

        lvl = 0.35;
        if nnz(V > lvl) < 50
            lvl = 0.25;
        end

        if any(V(:) > lvl)
            patch(isosurface(V,lvl), ...
                'FaceColor', colors(j,:), ...
                'EdgeColor','none', ...
                'FaceAlpha',0.96, ...
                'HandleVisibility','off');

            patch(isocaps(V,lvl), ...
                'FaceColor', colors(j,:), ...
                'EdgeColor','none', ...
                'FaceAlpha',0.25, ...
                'HandleVisibility','off');

            drew_any = true;
        end
    end

    % fallback: plot xPhys if material assignment is empty
    if ~drew_any
        Vx = reshape(xPhys, [nely, nelz, nelx]);
        Vx = permute(Vx, [3 1 2]);
        Vx = smooth3(Vx,'box',5);

        lvl = 0.3;
        if any(Vx(:) > lvl)
            patch(isosurface(Vx,lvl), ...
                'FaceColor',[0.3 0.7 1], ...
                'EdgeColor','none', ...
                'FaceAlpha',0.95, ...
                'HandleVisibility','off');
        end
    end

    view(ax,[140 15]);
    axis(ax,'equal'); axis(ax,'off'); daspect(ax,[1 1 1]);
    lighting(ax,'gouraud'); material(ax,'shiny');
    camlight(ax,'headlight'); camlight(ax,'right');

    title(ax, sprintf('MMTO (Iter %d)',loop), ...
        'FontName','Times New Roman','FontSize',14,'FontWeight','bold');

    % legend
    hleg = gobjects(nm,1);
    for j = 1:nm
        hleg(j) = plot3(ax,NaN,NaN,NaN,'s', ...
            'MarkerSize',10, ...
            'MarkerFaceColor',colors(j,:), ...
            'MarkerEdgeColor','none', ...
            'DisplayName',sprintf('Material %d',j));
    end
    legend(ax,hleg,'Location','northeastoutside'); legend boxoff;

    drawnow limitrate;
end
if mod(loop,10)==0

    exportgraphics(figure(1),...
        fullfile(figDir,...
        sprintf('Fig1_Iter_%04d.png',loop)),...
        'Resolution',500);

end
% ============================================================
% Figure 10: Separated hard materials with auto-zoom
% ============================================================
if mod(loop,10)==0

    figure(10); clf;
    set(gcf,'Color','w','Renderer','opengl',...
        'Position',[100 100 1400 420]);

    [psiMax, matID] = max(psiV, [], 2);

    rhoShow = xPhys > 0.30;
    ownThr  = 0.25;
    showEl  = rhoShow & (psiMax > ownThr);

    matID(~showEl) = 0;

    for j = 1:nm
        ax = subplot(1,nm,j); hold(ax,'on');

        Vj = double(matID == j);
        Vj3 = reshape(Vj,[nely,nelz,nelx]);
        Vj3 = permute(Vj3,[3 1 2]);      % x-y-z
        Vj3 = smooth3(Vj3,'box',5);

        lvl = 0.35;
        if nnz(Vj3 > lvl) < 50
            lvl = 0.25;
        end

        if any(Vj3(:) > lvl)

            p = patch(ax,isosurface(Vj3,lvl));
            isonormals(Vj3,p);
            set(p,'FaceColor',colors(j,:),...
                  'EdgeColor','none',...
                  'FaceAlpha',0.98);

            pc = patch(ax,isocaps(Vj3,lvl));
            set(pc,'FaceColor',colors(j,:),...
                   'EdgeColor','none',...
                   'FaceAlpha',0.20);

        else
            text(ax,0.5,0.5,'No region',...
                'Units','normalized',...
                'HorizontalAlignment','center',...
                'FontWeight','bold');
        end

        view(ax,[140 20]);
        axis(ax,'equal'); axis(ax,'off');
        daspect(ax,[1 1 1]);

        camlight(ax,'headlight');
        camlight(ax,'right');
        lighting(ax,'gouraud');
        material(ax,'shiny');

        Mi = sum(matID == j)/nelem;
        title(ax,sprintf('Material %d (M=%.3f)',j,Mi),...
            'FontWeight','bold','FontSize',10);
    end

    sgtitle(sprintf('Separated Hard Materials (Iter %d)',loop),...
        'FontWeight','bold','FontSize',15);

    drawnow;

    exportgraphics(gcf,...
        fullfile(figDir,sprintf('Fig10_Iter_%04d.png',loop)),...
        'Resolution',500);
end
% ============================================================
% Figure 2: Double-symmetry full bridge, clean material plot
% ============================================================
if mod(loop,1)==0
    disp('--- Rendering clean double-symmetry full multi-material bridge ---');

    figure(2); clf;
    set(gcf,'Color','w','Renderer','opengl','Position',[100 100 1200 540]);
    ax = axes; hold(ax,'on');

    for j = 1:nm
        Vq = reshape(psiV(:,j), [nely, nelz, nelx]);   % y-z-x
        Vq = permute(Vq, [3 1 2]);              % x-y-z

        % Double symmetry expansion
        isovals = zeros(2*nelx, 2*nely, nelz);

        isovals(nelx+1:2*nelx, 1:nely, 1:nelz) = Vq;
        isovals(1:nelx, 1:nely, 1:nelz) = ...
            isovals(2*nelx:-1:nelx+1, 1:nely, 1:nelz);

        isovals(1:nelx, nely+1:2*nely, 1:nelz) = ...
            isovals(1:nelx, nely:-1:1, 1:nelz);

        isovals(nelx+1:2*nelx, nely+1:2*nely, 1:nelz) = ...
            isovals(nelx+1:2*nelx, nely:-1:1, 1:nelz);

        isovals = smooth3(isovals,'box',5);

        lvl = 0.35;
        if nnz(isovals > lvl) < 50
            lvl = 0.25;
        end

        if any(isovals(:) > lvl)
            patch(isosurface(isovals,lvl), ...
                'FaceColor', colors(j,:), ...
                'EdgeColor','none', ...
                'FaceAlpha',0.96, ...
                'HandleVisibility','off');

            patch(isocaps(isovals,lvl), ...
                'FaceColor', colors(j,:), ...
                'EdgeColor','none', ...
                'FaceAlpha',0.25, ...
                'HandleVisibility','off');
        end
    end

    view(ax,[140 15]);
    axis(ax,[1 2*nelx 1 2*nely 1 nelz]);
    axis(ax,'manual'); axis(ax,'equal'); axis(ax,'off');
    daspect(ax,[1 1 1]);

    lighting(ax,'gouraud'); material(ax,'shiny');
    camlight(ax,'headlight'); camlight(ax,'right');

    title(ax, sprintf('MMTO (Iter %d)', loop), ...
          'FontName','Times New Roman','FontSize',14,'FontWeight','bold');

    % clean legend
    hleg = gobjects(nm,1);
    for j = 1:nm
        hleg(j) = plot3(ax,NaN,NaN,NaN,'s', ...
            'MarkerSize',10, ...
            'MarkerFaceColor',colors(j,:), ...
            'MarkerEdgeColor','none', ...
            'DisplayName',sprintf('Material %d',j));
    end
    legend(ax,hleg,'Location','northeastoutside');
    legend boxoff;

    drawnow limitrate;

    % Optional: write video frame if VideoWriter v exists and is open
    if exist('v','var')
        try
            writeVideo(v, getframe(gcf));
        catch
            % ignore video writing error to keep optimization running
        end
    end
end
% ===== Save image every 10 iterations =====
if mod(loop,10)==0

    fname = fullfile(figDir,...
        sprintf('Bridge_Iter_%04d.png',loop));

    exportgraphics(gcf,fname,...
        'Resolution',500);

end
% ============================================================
% Figure 20: Convergence graph
% Left Y-axis  : Compliance
% Right Y-axis : Total mass fraction + material mass fractions
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

    plot(it, HistMassTotal(it), '-s', ...
        'Color',[0.25 0.25 0.25], ...
        'LineWidth',2, ...
        'MarkerSize',3, ...
        'MarkerFaceColor','w');

    for jj = 1:mnum
        plot(it, HistMassFrac(it,jj), ':', ...
            'Color',colors(jj,:), ...
            'LineWidth',2);
    end

    ylabel('Mass fraction','FontName','Times New Roman','FontSize',16);

    ax = gca;
    ax.FontName = 'Times New Roman';
    ax.FontSize = 16;
    ax.LineWidth = 1.0;

    xlim([1 max(loop,2)]);

    %yyaxis right
    %maxMass = max([HistMassTotal(it); HistMassFrac(it,:)],[],'all','omitnan');
    yyaxis right
    
    massTotalVec = HistMassTotal(it);
    massFracMat  = HistMassFrac(it,:);
    
    massAll = [massTotalVec(:); massFracMat(:)];
    massAll = massAll(isfinite(massAll));
    
    if isempty(massAll)
        maxMass = 0.35;
    else
        maxMass = max(massAll);
    end
    
    ylim([0 max(0.35,1.15*maxMass)]);

    yyaxis left
    minC = min(HistComp(it),[],'omitnan');
    maxC = max(HistComp(it),[],'omitnan');
    if isfinite(minC) && isfinite(maxC) && maxC > minC
        ylim([0.95*minC 1.05*maxC]);
    end

    legText = [{'Compliance','Total Mass Fraction'}, ...
               arrayfun(@(jj)sprintf('Mass Fraction %d',jj),1:mnum,'UniformOutput',false)];

    legend(legText,'Location','northeast','FontName','Times New Roman','FontSize',16);

    title(sprintf('Convergence history (Iter %d)',loop), ...
        'FontName','Times New Roman','FontSize',16,'FontWeight','bold');

    drawnow limitrate;
end


end % ===== end while loop =====
%% ===== Save final xPhys result to MAT file =====
if ~exist(figDir,'dir')
    mkdir(figDir);
end

save(fullfile(figDir, filename), ...
     'xPhys','x','psiV','psiE','HistComp','HistMassTotal','HistMassFrac', ...
     'nelx','nely','nelz','volfrac','Evec','Dvec','-v7.3');

fprintf('\nSaved final xPhys and MMTO results to: %s\n', ...
        fullfile(figDir, filename));
%% ======== nested helpers ========
function y = cont_update(v,cn,iter)
  y = min(v + v*(iter>=cn{3})*(mod(iter-cn{3},cn{4})==0)*cn{5}, cn{2});
end

function [h,dh] = heaviside_smooth(xT, tgrid, beta)
% xT: ne×1 -> h,dh: ne×m
  ne = numel(xT); m = numel(tgrid);
  h  = zeros(ne,m); dh = zeros(ne,m);
  for j=1:m
    s = (xT - tgrid(j))*beta;
    th = tanh(s);
    h(:,j)  = 0.5*(1 + th);          % (0,1)
    dh(:,j) = 0.5*beta*(1 - th.^2);  % d h_j / d xT
  end
end

function [psi,dpsi] = MatIntModel(h,dh,pen)
% Hierarchical product model (partition-of-unity)
  [ne,m] = size(h);
  hp  = h.^pen;
  dhp = pen*(h.^(pen-1)).*dh;
  psi = zeros(ne,m);
  dpsi = cell(m,1); for j=1:m, dpsi{j} = zeros(ne,m); end

  for i=1:m
    if i < m
      hm = [hp(:,1:i), 1 - hp(:,i+1)];      % ne×(i+1)
      psi(:,i) = prod(hm,2);
      for j=1:m
        dhm = hm;
        if j >= i+2
          dhm(:,j) = 0;
        elseif j == i+1
          dhm(:,j) = -dhp(:,j);
        else
          dhm(:,j) = dhp(:,j);
        end
        dpsi{j}(:,i) = prod(dhm,2);
      end
    else % i == m
      hm = hp; 
      psi(:,i) = prod(hm,2);
      for j=1:m
        dhm = hm; dhm(:,j) = dhp(:,j);
        dpsi{j}(:,i) = prod(dhm,2);
      end
    end
  end
end

end  % ===== end main =====

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