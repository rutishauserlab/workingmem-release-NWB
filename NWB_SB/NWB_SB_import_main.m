%% NWB_SB_import_main
% A sample import script that computes the following:
% - Behavioral metrics
% - Spike sorting metrics
% - SU selectivity metrics (SB)
% - SU selectivity metrics (SC)
%
% SB refers to sternberg main task, and SC to sternberg screen task
% throughout.
%
%
% Michael Kyzar 6/30/23
 
clear; clc; close all
fs = filesep;
%% Parameters

% Operation Flags: Should either be  '1' (SCREENING), '2' (STERNBERG), or '3' (BOTH)
taskFlag = 2;

% subject IDs for dataset.
% importRange = 1:21; % Full Dataset: Kyzar et al 2023
importRange = [7 14 16]; % Sternberg Examples
% importRange = [4 7 15 16 21]; % Screening Examples
% importRange = 6:19; % Dataset: Kaminski et al 2017 
% importRange = 1:20; % Dataset: Kaminski et al 2020

%% Initializing and pathing
paths.baseData = 'C:\data\000469'; % Dataset directory
paths.nwb_sb = paths.baseData; % Dandiset Directory
paths.nwb_sc = paths.baseData; % Dandiset Directory
paths.code = 'C:\repos\workingmem-release-NWB';
paths.matnwb = 'C:\repos\matnwb-2.6.0.2';
paths.figOut = ['C:\temp' fs 'figures'];

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
load_all_waveforms = 1; % Extracts all by default. Set to '0' to only extract the mean waveform. 
switch taskFlag
    case 1 % Screening
        fprintf('Loading Screening\n')
        all_units_sc = NWB_SB_extractUnits(nwbAll_sc,load_all_waveforms);
    case 2 % Sternberg
        fprintf('Loading Sternberg\n')
        all_units_sb = NWB_SB_extractUnits(nwbAll_sb,load_all_waveforms);    
    case 3 % Both
        fprintf('Loading Screening\n')
        all_units_sc = NWB_SB_extractUnits(nwbAll_sc,load_all_waveforms);
        fprintf('Loading Sternberg\n')
        all_units_sb = NWB_SB_extractUnits(nwbAll_sb,load_all_waveforms);
    otherwise
        error('Task flag not properly specified')
end

%% SCREENING Params
paramsSC.doPlot = 1; % if =1, plot significant cells. 
paramsSC.plotAlways = 0; % Plot regardless of selectivity (warning: generates a lot of figures unless exportFig=1)
paramsSC.exportFig = 0;
paramsSB.exportType = 'png'; % File type for export. 'png' is the default. 
paramsSC.rateFilter = []; % Rate filter in Hz
paramsSC.figOut = [paths.figOut fs 'stats_screening_concept_example_eps'];
%% SCREENING Stats. Loops over all cells. 
paramsSC.calcSelective = 0;
if any(taskFlag == [1,3]) && paramsSC.calcSelective
    [sig_cells_sc, areas_sc] = NWB_calcSelective_SC(nwbAll_sc,all_units_sc,paramsSC);
end
%% SCREENING Examples. Loops over cells in Kyzar et al 2023. 
% For speed, specify the import range as [4 7 15 16 21] beforehand.
paramsSC.processExamples = 1;
if any(taskFlag == [1,3]) && paramsSC.processExamples
    % Kyzar et. al. 2023, Figure 5
    [sig_cells_sc_examples, areas_sc_examples] = NWB_SB_plotCell_Screening(nwbAll_sc,all_units_sc,paramsSC);
end

%% STERNBERG Params
paramsSB.doPlot = 1;  % if =1, plot significant cells. 
paramsSB.plotAlways = 0; % Plot regardless of selectivity (warning: generates a lot of figures unless exportFig=1)
paramsSB.plotMode = 1; % Which cell type to plot (1: Concept, 2: Maint, 3: Probe, 4: All)
paramsSB.exportFig = 0; 
paramsSB.exportType = 'png'; % File type for export. 'png' is the default. 
paramsSB.rateFilter =  0; % Rate filter in Hz. Setting to zero disables the filter. 
paramsSB.figOut = [paths.figOut fs 'stats_sternberg'];

%% STERNBERG Stats Loops over all cells. 
paramsSB.calcSelective = 0;
paramsSB.plotMode = 1; % Specify which type of cell to plot if only params.doPlot is 1 (1: Concept, 2: Maint, 3: Probe, 4: All types)
if any(taskFlag == [2,3]) && paramsSB.calcSelective
    [sig_cells_sb, areas_sb] = NWB_calcSelective_SB(nwbAll_sb,all_units_sb,paramsSB);
end

%% STERNBERG Examples. Loops over cells in Kyzar et al 2023. 
% For speed, specify the import range as [7 14 16] beforehand.
paramsSB.processExamples = 1;
if any(taskFlag == [2,3]) && paramsSB.processExamples
    [sig_cells_sc, areas_sc] = NWB_SB_plotCell_Sternberg(nwbAll_sb,all_units_sb,paramsSB);
end

%% Restate Metric Summaries:

if any(taskFlag == [1,3]) && paramsSC.calcSelective
    fprintf('SC: Concept Cells %d/%d\n\n',sum(sig_cells_sc),length(sig_cells_sc))
end
if any(taskFlag == [2,3]) && paramsSB.calcSelective
    fprintf('SB: Concept Cells %d/%d\n',sum(sig_cells_sb.concept_cells),length(sig_cells_sb.concept_cells))
    fprintf('SB: Maint Cells %d/%d\n',sum(sig_cells_sb.maint_cells),length(sig_cells_sb.maint_cells))
    fprintf('SB: Probe Cells %d/%d\n',sum(sig_cells_sb.probe_cells),length(sig_cells_sb.probe_cells))
end

%% Plot Significant Areas
plotAreas = 0;
if taskFlag == 3 && plotAreas
    
    [sb_areaCount_all, sb_labels_sig, sb_areaCount_sig, sb_labels_tot] = NWB_countAreas(areas_sb,sig_cells_sb.concept_cells);
    [sc_areaCount_all, sc_labels_sig, sc_areaCount_sig, sc_labels_tot] = NWB_countAreas(areas_sc,sig_cells_sc);
    [sb_areaCount_all_maint, sb_labels_sig_maint, sb_areaCount_sig_maint, sb_labels_tot_maint] = NWB_countAreas(areas_sb,sig_cells_sb.maint_cells);
    [sb_areaCount_all_probe, sb_labels_sig_probe, sb_areaCount_sig_probe, sb_labels_tot_probe] = NWB_countAreas(areas_sb,sig_cells_sb.probe_cells);

    sb_prop_selective = (sb_areaCount_sig./sb_areaCount_all)*100;
    sb_prop_selective_maint = (sb_areaCount_sig_maint./sb_areaCount_all_maint)*100;
    sc_prop_selective = (sc_areaCount_sig./sc_areaCount_all)*100;
    
    % Plotting selective proportion
    figure()
    bar([sb_prop_selective',sc_prop_selective'])
    ylim([0 100])
    yticks([0 50 100])
    
    legend({'sb','sc'})
    label_sb = ((string(sb_areaCount_sig') + '/') + num2str(sb_areaCount_all'))'; label_sb = strjoin(label_sb,' '); % Creates string of proportions for all areas
    label_sb_maint = ((string(sb_areaCount_sig_maint') + '/') + num2str(sb_areaCount_all_maint'))'; label_sb_maint = strjoin(label_sb_maint,' '); 
    label_sc = ((string(sc_areaCount_sig') + '/') + num2str(sc_areaCount_all'))'; label_sc = strjoin(label_sc,' ');
    set(gca, 'XTickLabel',sb_labels_tot)
    title('Proportion of selective units by area')
    xlabel(['SB: ' + label_sb 'SC: ' + label_sc])
    
    % Plotting total proportions (pie chart)
    colormap = [...
        0, 1, 0; ... % vmPFC
        0, 0, 1; ... % dACC
        1, 0, 0; ... % pre-SMA
        0, 1, 1; ... % Amg
        1, 1, 0 ...  % Hipp
    ];
    
    % Kyzar et. al. 2023, Figure 3.i-j
    fig = figure();
    
    subplot(2,1,1)
    pie(sb_areaCount_all, sb_labels_tot)
    title('SB: ' + label_sb)
    ax = gca();
    ax.Colormap = colormap;

    subplot(2,1,2)
    pie(sc_areaCount_all,sc_labels_tot)
    title('SC: ' + label_sc)
    ax = gca();
    ax.Colormap = colormap;

    %% Proportion of Selective Cells (Kyzar et al , Technical Validation)

    % SB: MTL Concept
    findMTL = find(contains(sb_labels_sig,'Amg') + contains(sb_labels_sig,'Hipp'));
    concept_MTL = sum(sb_areaCount_sig(findMTL));
    all_MTL = sum(sb_areaCount_all(findMTL));
    fprintf('MTL Concept Cells: %d/%d, (%.2f%%)\n',concept_MTL,all_MTL,concept_MTL/all_MTL*100)
    % SB: Persistent Firing, all concept cells
    hzPref_conceptOnly = sig_cells_sb.hzPref(logical(sig_cells_sb.concept_cells));
    hzNonPref_conceptOnly = sig_cells_sb.hzNonPref(logical(sig_cells_sb.concept_cells));

    [~, p_prefNonPref] = ttest(hzPref_conceptOnly,hzNonPref_conceptOnly,'Tail','right','Alpha',.05);
    fprintf('Persistent Firing : Pref %.2f ± %.2f | Non-Pref %.2f ± %.2f | p = %.4f\n',mean(hzPref_conceptOnly),std(hzPref_conceptOnly),mean(hzNonPref_conceptOnly),std(hzNonPref_conceptOnly), p_prefNonPref )
    % SB: MFC Maint
    findMFC_maint = find(contains(sb_labels_sig_maint,'vmPFC') + contains(sb_labels_sig_maint,'dACC') + contains(sb_labels_sig_maint,'pre-SMA'));
    maint_MFC = sum(sb_areaCount_sig_maint(findMFC_maint));
    maint_MFC_all = sum(sb_areaCount_all_maint(findMFC_maint));
    fprintf('MFC Maint Cells: %d/%d, (%.2f%%)\n',maint_MFC,maint_MFC_all,maint_MFC/maint_MFC_all*100)
    % SB: MFC probe
    findMFC_probe = find(contains(sb_labels_sig_probe,'vmPFC') + contains(sb_labels_sig_probe,'dACC') + contains(sb_labels_sig_probe,'pre-SMA'));
    probe_MFC = sum(sb_areaCount_sig_probe(findMFC_probe));
    probe_MFC_all = sum(sb_areaCount_all_probe(findMFC_probe));
    fprintf('MFC probe Cells: %d/%d, (%.2f%%)\n',probe_MFC,probe_MFC_all,probe_MFC/probe_MFC_all*100)
    % SC: MTL Concept
    findMTL = find(contains(sc_labels_sig,'Amg') + contains(sc_labels_sig,'Hipp'));
    concept_MTL = sum(sc_areaCount_sig(findMTL));
    all_MTL = sum(sc_areaCount_all(findMTL));
    fprintf('MTL Concept Cells (Screening): %d/%d, (%.2f%%)\n',concept_MTL,all_MTL,concept_MTL/all_MTL*100)
    % SC: MFC Concept
    findMFC = find(contains(sc_labels_sig,'vmPFC') + contains(sc_labels_sig,'dACC') + contains(sc_labels_sig,'pre-SMA'));
    concept_MFC = sum(sc_areaCount_sig(findMFC));
    all_MFC = sum(sc_areaCount_all(findMFC));
    fprintf('MFC Concept Cells (Screening): %d/%d, (%.2f%%)\n',concept_MFC,all_MTL,concept_MFC/all_MFC*100)

end

%% Calculate Spike Sorting Metrics 
calcMetrics = 0;
if any(taskFlag == [2,3]) && calcMetrics
    % Kyzar et. al., Figure 3.a-h
    is_sternberg = true;
    QAfig_sb = NWB_QA_graphs(nwbAll_sb, all_units_sb, is_sternberg);
    QAfig_sb.set("Visible","on")
    % Moving to left screen
    QAfig_sb.WindowState = 'maximized';
    movegui(QAfig_sb,'west') 
    QAfig_sb.WindowState = 'maximized';
end
if any(taskFlag == [1,3]) && calcMetrics % Calc metrics for Screening, no behavior
    % Kyzar et. al., Figure 3.a-h
    is_sternberg = false;
    QAfig_sc = NWB_QA_graphs(nwbAll_sc, all_units_sc, is_sternberg);
    QAfig_sc.set("Visible","on")
    % Moving to left screen
    QAfig_sc.WindowState = 'maximized';
    movegui(QAfig_sc,'west') 
    QAfig_sc.WindowState = 'maximized';
end


