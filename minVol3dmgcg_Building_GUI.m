function minVol3dmgcg_Building_GUI
clc;

fig = uifigure('Name','MMTO Building GUI','Position',[60 40 1350 760]);
fig.Color = [0.97 0.97 0.97];

mainGrid = uigridlayout(fig,[1 2]);
mainGrid.ColumnWidth = {360,'1x'};
mainGrid.Padding = [8 8 8 8];

%% ===== LEFT CONTROL PANEL =====
left = uipanel(mainGrid,'Title','Input Parameters','FontWeight','bold');
left.Layout.Row = 1;
left.Layout.Column = 1;

grid = uigridlayout(left,[18 2]);
grid.RowHeight = repmat({29},1,18);
grid.ColumnWidth = {135,'1x'};
grid.Padding = [12 12 12 12];
grid.RowSpacing = 7;

uilabel(grid,'Text','nelx');
nelxBox = uieditfield(grid,'numeric','Value',40,'RoundFractionalValues','on');

uilabel(grid,'Text','nely');
nelyBox = uieditfield(grid,'numeric','Value',40,'RoundFractionalValues','on');

uilabel(grid,'Text','nelz');
nelzBox = uieditfield(grid,'numeric','Value',60,'RoundFractionalValues','on');

uilabel(grid,'Text','Number of slabs');
nSlabBox = uieditfield(grid,'numeric','Value',6,'Limits',[1 Inf],'RoundFractionalValues','on');

uilabel(grid,'Text','Volume fraction');
volBox = uieditfield(grid,'numeric','Value',0.20,'Limits',[0.001 1]);

uilabel(grid,'Text','rmin');
rminBox = uieditfield(grid,'numeric','Value',sqrt(3),'Limits',[0.1 Inf]);

uilabel(grid,'Text','Filter type');
ftDrop = uidropdown(grid,...
    'Items',{'1: sensitivity filter','2: density filter'},...
    'ItemsData',[1 2],...
    'Value',1);

uilabel(grid,'Text','MG levels nl');
nlBox = uieditfield(grid,'numeric','Value',3,'RoundFractionalValues','on');

uilabel(grid,'Text','CG tolerance');
cgtolBox = uieditfield(grid,'numeric','Value',1e-6);

uilabel(grid,'Text','CG max iter');
cgmaxBox = uieditfield(grid,'numeric','Value',50,'RoundFractionalValues','on');

uilabel(grid,'Text','Output MAT file');
fileBox = uieditfield(grid,'text','Value','building_mmto_result.mat');

uilabel(grid,'Text','Evec');
EBox = uieditfield(grid,'text','Value','[0.25 0.50 0.70 0.85 1]');

uilabel(grid,'Text','Dvec');
DBox = uieditfield(grid,'text','Value','[0.20 0.40 0.60 0.80 1]');

uilabel(grid,'Text','Preset materials');
presetDrop = uidropdown(grid,...
    'Items',{'5 materials','3 materials','2 materials'},...
    'Value','5 materials',...
    'ValueChangedFcn',@applyPreset);

note = uitextarea(grid,'Editable','off');
note.Layout.Row = [15 16];
note.Layout.Column = [1 2];
note.Value = { ...
    'Tip: nelx, nely, nelz should be divisible by 2^(nl-1).', ...
    'Number of slabs controls passive horizontal floors.'};

runBtn = uibutton(grid,'push',...
    'Text','Run Optimization',...
    'FontWeight','bold',...
    'BackgroundColor',[0.2 0.6 0.2],...
    'FontColor',[1 1 1],...
    'ButtonPushedFcn',@runOptimization);
runBtn.Layout.Row = 17;
runBtn.Layout.Column = 1;

exitBtn = uibutton(grid,'push',...
    'Text','Exit',...
    'FontWeight','bold',...
    'BackgroundColor',[0.8 0.1 0.1],...
    'FontColor',[1 1 1],...
    'ButtonPushedFcn',@exitProgram);
exitBtn.Layout.Row = 17;
exitBtn.Layout.Column = 2;

status = uilabel(grid,...
    'Text','Ready',...
    'FontWeight','bold');
status.Layout.Row = 18;
status.Layout.Column = [1 2];

%% ===== RIGHT DISPLAY PANEL =====
right = uigridlayout(mainGrid,[2 1]);
right.Layout.Row = 1;
right.Layout.Column = 2;
right.RowHeight = {'2.3x','1.2x'};
right.RowSpacing = 10;

p1 = uipanel(right,...
    'Title','Figure 1: MMTO Topology',...
    'FontWeight','bold');
p1.Layout.Row = 1;

topoGrid = uigridlayout(p1,[1 1]);
topoGrid.Padding = [5 5 5 5];

axTopo = uiaxes(topoGrid);
view(axTopo,[140 15]);
axis(axTopo,'equal');
axis(axTopo,'off');
title(axTopo,'Topology');

p2 = uipanel(right,...
    'Title','Figure 2: Convergence History',...
    'FontWeight','bold');
p2.Layout.Row = 2;

convGrid = uigridlayout(p2,[1 1]);
convGrid.Padding = [5 5 5 5];

axConv = uiaxes(convGrid);
xlabel(axConv,'Iteration number');
ylabel(axConv,'Compliance');
title(axConv,'Convergence History');

%% ===== Callback functions =====
    function applyPreset(~,~)
        switch presetDrop.Value
            case '5 materials'
                EBox.Value = '[0.25 0.50 0.70 0.85 1]';
                DBox.Value = '[0.20 0.40 0.60 0.80 1]';
            case '3 materials'
                EBox.Value = '[0.25 0.70 1]';
                DBox.Value = '[0.20 0.60 1]';
            case '2 materials'
                EBox.Value = '[0.50 1]';
                DBox.Value = '[0.40 1]';
        end
    end

    function runOptimization(~,~)
        try
            runBtn.Enable = 'off';
            exitBtn.Enable = 'off';
            status.Text = 'Running...';
            drawnow;

            nelx = round(nelxBox.Value);
            nely = round(nelyBox.Value);
            nelz = round(nelzBox.Value);
            nSlab = round(nSlabBox.Value);

            volfrac = volBox.Value;
            rmin = rminBox.Value;
            ft = ftDrop.Value;
            nl = round(nlBox.Value);
            cgtol = cgtolBox.Value;
            cgmax = round(cgmaxBox.Value);
            filename = strtrim(fileBox.Value);

            Evec = str2num(EBox.Value); %#ok<ST2NM>
            Dvec = str2num(DBox.Value); %#ok<ST2NM>

            validateattributes(Evec,{'numeric'},{'vector','nonempty','positive'});
            validateattributes(Dvec,{'numeric'},{'vector','nonempty','positive'});

            if numel(Evec) ~= numel(Dvec)
                error('Evec and Dvec must have the same number of materials.');
            end

            if nSlab < 1 || nSlab > nelz
                error('Number of slabs must be between 1 and nelz.');
            end

            if isempty(filename) || ~endsWith(filename,'.mat')
                error('Output filename must end with .mat');
            end

            if mod(nelx,2^(nl-1)) ~= 0 || mod(nely,2^(nl-1)) ~= 0 || mod(nelz,2^(nl-1)) ~= 0
                error('nelx, nely, and nelz must be divisible by 2^(nl-1).');
            end

            cla(axTopo,'reset');
            cla(axConv,'reset');

            minVol3dmgcg_Building( ...
                nelx,nely,nelz,nSlab,volfrac,rmin,ft,nl,cgtol,cgmax,...
                filename,Evec,Dvec,axTopo,axConv,status);

            status.Text = 'Finished. Check Results_* folder.';
            uialert(fig,'Optimization finished. Check Results_* folder.','Done');

        catch ME
            status.Text = ['Error: ' ME.message];
            uialert(fig,ME.message,'Run error');
        end

        runBtn.Enable = 'on';
        exitBtn.Enable = 'on';
    end

    function exitProgram(~,~)
        selection = uiconfirm(fig,...
            'Do you want to exit MMTO GUI?',...
            'Exit Confirmation',...
            'Options',{'Yes','No'},...
            'DefaultOption',2,...
            'CancelOption',2);

        if strcmp(selection,'Yes')
            delete(fig);
        end
    end
end