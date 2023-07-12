%% SB_NWB_import_main_v2
% Second iteration of the import script that computes the following:
% - Behavioral metrics
% - Spike sorting metrics
% - SU selectivity metrics (SB)
% - SU selectivity metrics (SC)
%
% SB refers to sternberg main task, and SC to sternberg screen task
% throughout.
%
% Michael Kyzar 4/5/2023

% NOTE: Being used for main analysis. For QA metrics, see QA_graphs.mat
% May need the most recent version of MATLAB (2023a) to run the anova
% stats. 
clear; clc; close all
%% Parameters

fs = filesep;

% Operation Flags: Should either be  'SCREENING' (1), 'STERNBERG' (2),or 'BOTH' (3)
taskFlag = 3;

% importRange = 1:2; % subject IDs for dataset. Empty defaults to all available subjects.
% importRange = 6:19; % K2017 dataset
% importRange = 1:20; % K2020 dataset
 importRange = 1:21; % Full Dataset
%importRange = 20; % Testing set. 

calcSelective = 1;
calcMetrics = 1;


%% Initializing and pathing

%paths.baseData = 'N:\LabUsers\kyzarm\data\NWB_SB\data_NWB'; % Dataset directory
% paths.nwb_sb = [paths.baseData fs 'STERNBERG']; % Native Directory
% paths.nwb_sc = [paths.baseData fs 'SCREENING']; % Native Directory

%paths.nwb_sb = 'N:\LabUsers\kyzarm\data\NWB_SB\data_Dandiset\000469'; % Dandiset testing
%paths.nwb_sc = 'N:\LabUsers\kyzarm\data\NWB_SB\data_Dandiset\000469'; % Dandiset testing

%to download: dandi download https://dandiarchive.org/dandiset/000469/draft
paths.nwb_sb = 'C:\data\000469'; % Dandiset testing
paths.nwb_sc = 'C:\data\000469'; % Dandiset testing

% Paths to code
%paths.code = 'C:\svnwork\neuro1\code\events\Sternberg\NWB-SB\NWB_SB_import_release';
paths.matnwb = 'C:\svnwork\matnwb-2.6.0.2';

paths.figOut = ['C:\temp\NWB_SB' fs 'figures'];

% Helpers
if(~isdeployed) 
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
  addpath(genpath([pwd fs 'helpers'])) % Should be in same folder as active script. 
else
    error('Unexpected error.')
end

pathCell = struct2cell(paths);
for i = 1:length(pathCell)
    addpath(genpath(pathCell{i}))
end

% Initialize NWB Package
%generateCore() for first instantiation of matnwb API
fprintf('Checking generateCore() ... ')
if isfile([paths.matnwb fs '+types' fs '+core' fs 'NWBFile.m'])
     fprintf('generateCore() already initialized.\n') %only need to do once
else 
    cd(paths.matnwb)
    generateCore();
    fprintf('generateCore() initialized.\n')
end 

%% Importing Datasets From Folder
% Note: Parallel loading is not yet supported due to the pools referencing
% matnwb files that are initialized outside the main matnwb folder. Could
% be a suggestion for future releases of matnwb. 
switch taskFlag
    case 1 % Screening
        internalFlag = taskFlag;
        [nwbAll_sc, importLog_sc] = NWB_importFromFolder(paths.nwb_sc, importRange,internalFlag);
    case 2 % Sternberg
        internalFlag = taskFlag;
        [nwbAll_sb, importLog_sb] = NWB_importFromFolder(paths.nwb_sb, importRange, internalFlag);    
    case 3 % Both
        internalFlag = 1;
        [nwbAll_sc, importLog_sc] = NWB_importFromFolder(paths.nwb_sc, importRange, internalFlag);
        internalFlag = 2;
        [nwbAll_sb, importLog_sb] = NWB_importFromFolder(paths.nwb_sb, importRange, internalFlag);
    otherwise
        error('Task flag not properly specified')
end

%% Extracting Single Units
switch taskFlag
    case 1 % Screening
        fprintf('Loading Screening\n')
        all_units_sc = NWB_SB_extractUnits(nwbAll_sc);
    case 2 % Sternberg
        fprintf('Loading Sternberg\n')
        all_units_sb = NWB_SB_extractUnits(nwbAll_sb);    
    case 3 % Both
        fprintf('Loading Screening\n')
        all_units_sc = NWB_SB_extractUnits(nwbAll_sc);
        fprintf('Loading Sternberg\n')
        all_units_sb = NWB_SB_extractUnits(nwbAll_sb);
    otherwise
        error('Task flag not properly specified')
end

%% SCREENING Stats. Loops over all cells.
paramsSC.doPlot = 0; % if =1, plot significant cells. 
paramsSC.plotAlways = 0;
paramsSC.exportFig = 0;
paramsSC.rateFilter = []; % Rate filter in Hz
paramsSC.runParallel = 0;
paramsSC.figOut = [paths.figOut fs 'stats_screening_test_out'];
if any(taskFlag == [1,3]) && calcSelective
    [sig_cells_sc, areas_sc] = NWB_calcSelective_SC(nwbAll_sc,all_units_sc,paramsSC);
end

%% STERNBERG Stats. Loops over all cells.
paramsSB.doPlot = 1;  % if =1, plot significant cells. 
paramsSB.plotAlways = 0;
paramsSB.exportFig = 0;
paramsSB.rateFilter = []; % Rate filter in Hz
paramsSB.runParallel = 0;
paramsSB.figOut = [paths.figOut fs 'stats_sternberg_test_out'];
if any(taskFlag == [2,3]) && calcSelective
    [sig_cells_sb, areas_sb] = NWB_calcSelective_SB(nwbAll_sb,all_units_sb,paramsSB);
end

%% Plot specific example cells shown in paper (Sternberg task)
% This ection plots example cells shown in the papers about this dataset (Sternberg task part)

% Kaminski 2017, Fig. 3A
%NWB_plotCell_SternbergMaintask(...)
    
% Kaminski 2017, Fig. 3B
% Kaminski 2020, Fig. 1E-G
% Kyzar et al 2023, Fig. XX

%% Plot specific example cells shown in paper (Screening task)
% This ection plots example cells shown in the papers about this dataset (Screening task part)

% Kaminski 2017, Fig. S1A,  top
%NWB_plotCell_SternbergScreenin(...)

% Kyzar et al 2023, Fig. XX

%% Metric Notes
%{
2017 Concept Cell Metrics:
SB: [Kaminsky2017: 93/651, 14.28%]
SC: [Kaminsky2017: 99/670, 14.77%]
2017 Replication Concept Cell Metics
SB2017: 93/651, 14.28%
SC2017: 133/681, 19.53%
Full Concept Cell Metrics:
SB: 119/902, 13.19%
SC: 164/907, 18.08%
%}
%% Plot Significant Areas
plotAreas = 1;
if taskFlag == 3 && plotAreas
    
    [sb_areaCount_all, sb_labels_sig, sb_areaCount_sig, sb_labels_tot] = NWB_countAreas(areas_sb,sig_cells_sb.concept_cells);
    [sc_areaCount_all, sc_labels_sig, sc_areaCount_sig, sc_labels_tot] = NWB_countAreas(areas_sc,sig_cells_sc);

    sb_prop_selective = (sb_areaCount_sig./sb_areaCount_all)*100;
    sc_prop_selective = (sc_areaCount_sig./sc_areaCount_all)*100;
    
    % Plotting selective proportion
    figure()
    bar([sb_prop_selective',sc_prop_selective'])
    ylim([0 100])
    yticks([0 50 100])
    
    legend({'sb','sc'})
    label_sb = ((string(sb_areaCount_sig') + '/') + num2str(sb_areaCount_all'))'; label_sb = strjoin(label_sb,' '); % Creates string of proportions for all areas
    label_sc = ((string(sc_areaCount_sig') + '/') + num2str(sc_areaCount_all'))'; label_sc = strjoin(label_sc,' ');
    set(gca, 'XTickLabel',sb_labels_tot)
    title('Proportion of selective units by area')
    xlabel(['SB: ' + label_sb 'SC: ' + label_sc])
    

    % Plotting total proportions (pie chart)
    colormap = [...
        0, 1, 0; ... % OFC
        0, 0, 1; ... % dACC
        1, 0, 0; ... % pre-SMA
        0, 1, 1; ... % Amg
        1, 1, 0 ...  % Hipp
    ];
    
    fig = figure();
    %tiledlayout('vertical')
    subplot(2,2,1);
    %nexttile
    pie(sb_areaCount_all, sb_labels_tot)
    title('SB: ' + label_sb)
    ax = gca();
    ax.Colormap = colormap;

    subplot(2,2,2);
    %nexttile
    pie(sc_areaCount_all,sc_labels_tot)
    title('SC: ' + label_sc)
    ax = gca();
    ax.Colormap = colormap;

end

%% Calculate Spike Sorting Metrics (Sternberg only)
if any(taskFlag == [1,3]) && calcMetrics
    QAfig = NWB_QA_graphs(nwbAll_sb, all_units_sb);
    QAfig.set("Visible","on")
    % Moving to left screen
    QAfig.WindowState = 'maximized';
    movegui(QAfig,'west') 
    QAfig.WindowState = 'maximized';
end










