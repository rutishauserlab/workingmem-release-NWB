function [sig_cells,areasSternberg] = NWB_calcSelective_SB(nwbAll, all_units, params)
%NWB_CALCSELECTIVE Takes the output of NWB_SB_extractUnits and runs
%selectivity tests for the sternberg task. Optionally plots trial aligned
%spikes

if isfield(params,'rateFilter') && params.rateFilter > 0
    rateFilter = params.rateFilter;
else
    rateFilter = [];
end

% Filtering for Global Rate (rateFilter should be a nonzero float)
aboveRate = ones(length(all_units),1);
if ~isempty(rateFilter)
    for i = 1:length(all_units)
        globalRate = length(all_units(i).spike_times)/(max(all_units(i).spike_times)-min(all_units(i).spike_times));
        rateBool = globalRate < rateFilter;
        if rateBool % If the rate is below the filter threshold
            aboveRate(i) = 0;
        end
    end
end
all_units = all_units(logical(aboveRate));


areasSternberg = cell(length(all_units),1);
concept_cells_sb = zeros(length(all_units),1); alphaLim = 0.05;
hzPref = zeros(length(all_units),1);
hzNonPref = zeros(length(all_units),1);
maint_cells_sb = zeros(length(all_units),1);
probe_cells_sb = zeros(length(all_units),1);
% Looping over all cells
for i = 1:length(all_units) 
    SU = all_units(i);
    subject_id = SU.subject_id;
    cellID = SU.unit_id;
    brain_area = nwbAll{SU.session_count}.general_extracellular_ephys_electrodes.vectordata.get('location').data.load(SU.electrodes);
    clusterID = nwbAll{SU.session_count}.units.vectordata.get('clusterID_orig').data.load(SU.unit_id);
    areasSternberg{i} = brain_area{:};
    fprintf('Processing: (%d/%d) Session SBID %d, Unit %d, Cluster %d ',i,length(all_units),subject_id,cellID,clusterID)
    
    % Loading stim timestamps and loads
    tsFix = num2cell(nwbAll{SU.session_count}.intervals_trials.vectordata.get('timestamps_FixationCross').data.load());
    tsEnc1 = num2cell(nwbAll{SU.session_count}.intervals_trials.vectordata.get('timestamps_Encoding1').data.load());
    tsEnc2 = num2cell(nwbAll{SU.session_count}.intervals_trials.vectordata.get('timestamps_Encoding2').data.load());
    tsEnc3 = num2cell(nwbAll{SU.session_count}.intervals_trials.vectordata.get('timestamps_Encoding3').data.load());
    tsMaint = num2cell(nwbAll{SU.session_count}.intervals_trials.vectordata.get('timestamps_Maintenance').data.load());
    tsProbe = num2cell(nwbAll{SU.session_count}.intervals_trials.vectordata.get('timestamps_Probe').data.load());
    ID_Enc1 = num2cell(nwbAll{SU.session_count}.intervals_trials.vectordata.get('loadsEnc1_PicIDs').data.load());
    ID_Enc2 = num2cell(nwbAll{SU.session_count}.intervals_trials.vectordata.get('loadsEnc2_PicIDs').data.load());
    ID_Enc3 = num2cell(nwbAll{SU.session_count}.intervals_trials.vectordata.get('loadsEnc3_PicIDs').data.load());
    ID_Probe = num2cell(nwbAll{SU.session_count}.intervals_trials.vectordata.get('loadsProbe_PicIDs').data.load());

    stim = cell2struct(...
        horzcat(tsFix,tsEnc1,tsEnc2,tsEnc3,tsMaint, tsProbe,...
        ID_Enc1,ID_Enc2,ID_Enc3,ID_Probe),...
        {'tsFix','tsEnc1','tsEnc2','tsEnc3','tsMaint','tsProbe',...
        'idEnc1','idEnc2','idEnc3','idProbe'},2);
    
    [idUnique, ~, ic] = unique([stim.idEnc1]);    
    uniqueCounts = histcounts(ic,'BinMethod','integers');
    
    HzEnc1 = NaN(length(stim),1);
    signalDelay = 0.2; % Delay of stimulus onset to effect. 
    stimOffset = 1; % Time past stimulus onset. End of picture presentation.
    for k = 1:length(HzEnc1)
        periodFilter = (SU.spike_times>(stim(k).tsEnc1+signalDelay)) & (SU.spike_times<(stim(k).tsEnc1+stimOffset));
        singleTrialSpikes = SU.spike_times(periodFilter);
        trialRate = length(singleTrialSpikes)/(stimOffset-signalDelay); % Firing rate across testing period.
        HzEnc1(k) = trialRate;
    end
    
    % Finding image with maximum response. Used in stage 2 test. 
    nUniqueStim = length(unique([stim.idEnc1]));
    mResp=nan(nUniqueStim,1);
    for k=1:nUniqueStim
        mResp(k) = mean(HzEnc1([stim.idEnc1]==k));
    end
    idMaxHz = find(mResp==max(mResp),1);
    id_Trial_maxOnly = [stim.idEnc1]; id_Trial_maxOnly([stim.idEnc1]~=idMaxHz) = -1;
    
    %% Significance Tests: Concept Cells
    % Count spikes in 200-1000ms window following stimulus onset
    % following the first encoding period. Use a 1-way ANOVA followed by a
    % t-test of the maximal response versus the non-selective responses.
    % Note that using a 1-way anova with two groups simplifies to a t-test
    p_ANOVA = 1; p_bANOVA = 1;%#ok<NASGU> % Preset as '1' to allow for paramsSB.plotAlways.
    p_ANOVA = anovan(HzEnc1,{string([stim.idEnc1])}, 'display','off','model', 'linear','alpha',alphaLim,'varnames','picID');
    if p_ANOVA < alphaLim % First test: 1-way ANOVA
        p_bANOVA = anovan(HzEnc1,{string(id_Trial_maxOnly)}, 'display','off','model', 'linear','alpha',alphaLim,'varnames','picID');
        if p_bANOVA < alphaLim % Second Test: Binarized 1-way ANOVA (simplifies to t-test)
        fprintf('| Concept -> SID %d, Unit %d p1:%.2f p2:%.2f',SU.subject_id,SU.unit_id,p_ANOVA,p_bANOVA)
        concept_cells_sb(i) = 1;
        end
    end
    % Saving rate data for additional tests between pref/non-pref trials
    hzPref(i) = mean(HzEnc1(logical(id_Trial_maxOnly)));
    hzNonPref(i) = mean(HzEnc1(~logical(id_Trial_maxOnly)));

    %% Significance Tests: Maintenance Cells
    % Perform a t-test between the mean firing rate during the maintenance
    % period and baseline period (fixation cross). If the cell was
    % idenfified as a concept cell, the maintenance activity of
    % non-selective images must be higher than the fixation baseline as
    % well.
    HzMaint = NaN(length(stim),1);
    maintDelay = 0; % Delay of maint onset to effect. 
    maintOffset = 2.5; % Time past maint onset. End of picture presentation.
    for k = 1:length(HzMaint)
        periodFilter = (SU.spike_times>(stim(k).tsMaint+maintDelay)) & (SU.spike_times<(stim(k).tsMaint+maintOffset));
        singleTrialSpikes = SU.spike_times(periodFilter);
        trialRate = length(singleTrialSpikes)/(maintOffset-maintDelay); % Firing rate across testing period.
        HzMaint(k) = trialRate;
    end
    HzFix = NaN(length(stim),1);
    fixDelay = 0; % Delay of fix onset to effect. 
    fixOffset = 0.5; % Time past fixation onset. End of picture presentation.
    for k = 1:length(HzFix)
        periodFilter = (SU.spike_times>(stim(k).tsFix+fixDelay)) & (SU.spike_times<(stim(k).tsFix+fixOffset));
        singleTrialSpikes = SU.spike_times(periodFilter);
        trialRate = length(singleTrialSpikes)/(fixOffset-fixDelay); % Firing rate across testing period.
        HzFix(k) = trialRate;
    end
    % % Performing paired t-test
    sidedness = 'right';
    [detect_maint, p_maint] = ttest(HzMaint,HzFix,'Tail',sidedness,'Alpha',alphaLim);
    if detect_maint && concept_cells_sb(i) % If the maint cell was previously marked as a concept cell
        non_selective_trials = ~id_Trial_maxOnly;
        [maint_cells_sb(i),p2_maint] = ttest(HzMaint(non_selective_trials),HzFix(non_selective_trials),'Tail',sidedness,'Alpha',alphaLim);
        if maint_cells_sb(i)
            fprintf('| Maint -> M:%.2fHz F:%.2fHz p:%.2f, p2:%.2f ',mean(HzMaint),mean(HzFix),p_maint,p2_maint)
        end
    elseif detect_maint && ~concept_cells_sb(i) % for pure maint cells. 
        maint_cells_sb(i) = 1;
        fprintf('| Maint -> M:%.2fHz F:%.2fHz p:%.2f ',mean(HzMaint),mean(HzFix),p_maint)
    end
    
    % % ALT METRIC: Permutation Test (Yields similar results)
    % permutations = 2000;
    % sidedness = 'larger'; % 'smaller': s1<s2; 'larger': s1>s2
    % [p_maint, ~, ~] = permutationTest(HzMaint,HzFix, permutations, 'sidedness', sidedness);
    % if p_maint < alphaLim && concept_cells_sb(i) % If the maint cell was previously marked as a concept cell
    %     non_selective_trials = ~id_Trial_maxOnly;
    %     [p2_maint,~,~] = permutationTest(HzMaint(non_selective_trials),HzFix(non_selective_trials), permutations, 'sidedness',sidedness);
    %     if p2_maint < alphaLim
    %         fprintf('| Maint -> M:%.2fHz F:%.2fHz p:%.2f, p2:%.2f ',mean(HzMaint),mean(HzFix),p_maint,p2_maint)
    %         maint_cells_sb(i) = 1;
    %     end
    % elseif p_maint < alphaLim && ~concept_cells_sb(i) % for pure maint cells. 
    %     maint_cells_sb(i) = 1;
    %     fprintf('| Maint -> M:%.2fHz F:%.2fHz p:%.2f ',mean(HzMaint),mean(HzFix),p_maint)
    % end


    %% Significance Tests: Probe Cells
    % Perform a t-test between the firing rates in the probe period and the
    % maint/enc periods. 
    HzProbe = NaN(length(stim),1);
    probeDelay = .2; % Delay of probe onset to effect. 
    probeOffset = 1; % Time past probe onset. End of picture presentation.
    for k = 1:length(HzProbe)
        periodFilter = (SU.spike_times>(stim(k).tsProbe+probeDelay)) & (SU.spike_times<(stim(k).tsProbe+probeOffset));
        singleTrialSpikes = SU.spike_times(periodFilter);
        trialRate = length(singleTrialSpikes)/(probeOffset-probeDelay); % Firing rate across testing period.
        HzProbe(k) = trialRate;
    end
    
    % Appending Enc/Maint periods for each trial.
    HzEnc = NaN(length(stim),1);
    encDelay = 0.2;
    encOffset = 1;
    % maintDelay = 0; % Delay of maint onset to effect. 
    % maintOffset = 2.5; % Time past maint onset. End of picture presentation.
    for k = 1:length(HzEnc)
        trialSpikeCounts = 0; totalTime = 0;
        % Enc 1
        periodFilter = (SU.spike_times>(stim(k).tsEnc1+encDelay)) & (SU.spike_times<(stim(k).tsEnc1+encOffset));
        singleTrialSpikes = SU.spike_times(periodFilter);
        trialSpikeCounts = trialSpikeCounts + length(singleTrialSpikes);
        totalTime = totalTime + (encOffset-encDelay);
        if stim(k).tsEnc2 > 0 % Enc 2
            periodFilter = (SU.spike_times>(stim(k).tsEnc2+encDelay)) & (SU.spike_times<(stim(k).tsEnc2+encOffset));
            singleTrialSpikes = SU.spike_times(periodFilter);
            trialSpikeCounts = trialSpikeCounts + length(singleTrialSpikes);
            totalTime = totalTime + (encOffset-encDelay);
        end
        if stim(k).tsEnc3 > 0 % Enc 3
            periodFilter = (SU.spike_times>(stim(k).tsEnc2+encDelay)) & (SU.spike_times<(stim(k).tsEnc2+encOffset));
            singleTrialSpikes = SU.spike_times(periodFilter);
            trialSpikeCounts = trialSpikeCounts + length(singleTrialSpikes);
            totalTime = totalTime + (encOffset-encDelay);
        end
        % % Maint (Deprecated, t-test calculated separately)
        % periodFilter = (SU.spike_times>(stim(k).tsMaint+maintDelay)) & (SU.spike_times<(stim(k).tsEnc2+maintOffset));
        % singleTrialSpikes = SU.spike_times(periodFilter);
        % trialSpikeCounts = trialSpikeCounts + length(singleTrialSpikes); 
        % totalTime = totalTime + (maintOffset-maintDelay);
       
        % Calculating Total Rate per trial
        trialRate = trialSpikeCounts/totalTime; % Firing rate across testing period.
        HzEnc(k) = trialRate;
    end

        
    % Performing paired one-tailed t-test
    [detect_probe1, p_probe1] = ttest(HzProbe,HzEnc,'Tail','right','Alpha',alphaLim);
    [detect_probe2, p_probe2] = ttest(HzProbe,HzMaint,'Tail','right','Alpha',alphaLim);
    if detect_probe1 && detect_probe2
        fprintf('| Probe -> P:%.2fHz E:%.2fHz M:%.2fHz p:%.2f p2:%.2f ',mean(HzProbe),mean(HzEnc),mean(HzMaint),p_probe1,p_probe2)
        probe_cells_sb(i) = 1;
    end

    

    %% Rasters, PSTH, & Selective Image
    % Flagging significant cells for plotting
    switch params.plotMode
        case 1 % Concept
            plotFlag = params.doPlot && concept_cells_sb(i);
        case 2 % Maint
            plotFlag = params.doPlot && maint_cells_sb(i);
        case 3 % Probe
            plotFlag = params.doPlot && probe_cells_sb(i);
        case 4 % All
            plotFlag = params.doPlot && (concept_cells_sb(i) || maint_cells_sb(i) || probe_cells_sb(i));
        otherwise
            warning('Plot mode not specified. Defaulting to concept cells.\n')
            in = input('Continue? (y/n)\n',"s");
            if any(strcmp(in,["Y","y"]))
                plotFlag = params.doPlot && concept_cells_sb(i);
            elseif any(strcmp(in,["N","n"]))
                fprintf('Aborting.\n')
                return
            else
                fprintf('Input not specified. Aborting.\n')
                return
            end
    end

    if plotFlag || params.plotAlways
        % Metrics
        subject_id = SU.subject_id;
        cellID = SU.unit_id;
        brain_area = nwbAll{SU.session_count}.general_extracellular_ephys_electrodes.vectordata.get('location').data.load(SU.electrodes);
        clusterID = nwbAll{SU.session_count}.units.vectordata.get('clusterID_orig').data.load(SU.unit_id);
        plabelStr = ['sub-' num2str(subject_id) '-ses-2'...
            ' | BA: ' strrep(brain_area{:},'_',' ') ...
            ' | Cell: ' num2str(cellID) ...
            ' | Elec: ' num2str(SU.electrodes) ...
            ' | p1=' num2str(p_ANOVA) ...
            ' | p2=' num2str(p_bANOVA)];
        
        % Set Params
        spikePlotParams=[];
        spikePlotParams.colors = { ...
            [1 0 0], ... % Red
            [.7 .7 .7], ... % Grey
            [0 0 1], ... % Blue
        	[1 0 1], ... % Magenta
            [0.4940 0.1840 0.5560] ... % Purple
            };  
        spikePlotParams.spikewidth  = 1.25;
        spikePlotParams.spikeheight = 1;
        spikePlotParams.binsizePSTH = .050;
        spikePlotParams.smoothKernelWidth = 1.5 * spikePlotParams.binsizePSTH;
        spikePlotParams.plotMarkerPos = 0;
        spikePlotParams.trimPadding = 1; % Used to overcome plotting artifacts encountered when using gaussian smoothing. 
        spikePlotParams.padding = 2 * spikePlotParams.smoothKernelWidth;


        % t before stim onset (In Seconds)
        StimBaselineEnc= .5; % for all encodings
        StimBaselineMaint = .5;  % before maintenance onset
        StimBaselineProbe = .5; % before probe onset
        % t after stim onset (In Seconds)
        StimOnTimeEnc  = 1; % after for all encodings
        StimOnTimeMaint = 2.5; % after maintenance onset
        StimOnTimeProbe = 2; % after probe onset

        
        % Fixes gaussian smoothing cutoff. 
        % % Offsets to include more data in gaussian smoothing kernel
        padding = .25;
        period_StimBaselineEnc= StimBaselineEnc + padding; % for all encodings
        period_StimBaselineMaint = StimBaselineMaint + padding;  % before maintenance onset
        period_StimBaselineProbe = StimBaselineProbe + padding; % before probe onset
        period_StimOnTimeEnc  = StimOnTimeEnc + padding; % after for all encodings
        period_StimOnTimeMaint = StimOnTimeMaint + padding; % after maintenance onset
        period_StimOnTimeProbe = StimOnTimeProbe + padding; % after probe onset

        % Defining Periods
        periods = [];
        periods.Enc1 = NWB_determineSBPeriods( nonzeros([stim.tsEnc1]), period_StimBaselineEnc, period_StimOnTimeEnc );
        periods.Enc2 = NWB_determineSBPeriods( nonzeros([stim.tsEnc2]), period_StimBaselineEnc, period_StimOnTimeEnc );
        periods.Enc3 = NWB_determineSBPeriods( nonzeros([stim.tsEnc3]), period_StimBaselineEnc, period_StimOnTimeEnc );
        periods.Maint = NWB_determineSBPeriods( nonzeros([stim.tsMaint]), period_StimBaselineMaint, period_StimOnTimeMaint );
        periods.Probe = NWB_determineSBPeriods( nonzeros([stim.tsProbe]), period_StimBaselineProbe, period_StimOnTimeProbe );
        
        

        figNum = i;
        currentFigure = figure(figNum);
        currentFigure.set('Name', plabelStr)
        %% === Encoding 1 == Plot Raster and PSTH
        offsets = 0; % Spike times offset
        indsPref = [stim.idEnc1]==idMaxHz; indsNonPref = ~indsPref;
        indsPref = find(indsPref); indsNonPref = find(indsNonPref);

        limitRange_forRaster = [ 0 StimBaselineEnc+StimOnTimeEnc + 2*padding ]; 
        limitRange_forPSTH = [ 0 StimBaselineEnc+StimOnTimeEnc + 2*padding ]; 
        markerPos=[StimBaselineEnc + padding  StimBaselineEnc + StimOnTimeEnc + padding ];
        params_Enc1 = spikePlotParams; params_Enc1.colors = params_Enc1.colors(1:2);
        [axRaster_Enc1,axPSTH_Enc1] = plotSpikeRaster_multiple_simplified( ...
            3, offsets, limitRange_forRaster, limitRange_forPSTH, ...
            markerPos, SU.spike_times, ...
            {periods.Enc1}, {indsPref, indsNonPref}, ...
            {'Preferred','Non-Preferred'},   [2 6 7], [2 6 1], params_Enc1);
        title('Enc1')
        xlabel(plabelStr)

        %% === Encoding 2 == Plot Raster and PSTH
        offsets = [0];
        limitRange_forRaster = [0 StimBaselineEnc+StimOnTimeEnc + 2*padding ]; 
        limitRange_forPSTH = [0 StimBaselineEnc+StimOnTimeEnc + 2*padding ];
        
        
        markerPos = [StimBaselineEnc + padding StimBaselineEnc+StimOnTimeEnc + padding ];% [StimBaselineEnc 1.5];

        indsEnc2_load2shown = find([ID_Enc2{:}]>0);    

        if length(indsEnc2_load2shown)>size(periods.Enc2,1)
            warning('more images then trial periods, missmatch?');
            indsEnc2_load2shown = indsEnc2_load2shown(1:size(periods.Enc2,1));
        end
        PicOrder_Enc2_shownOnly_inEnc2 = cell2mat(ID_Enc2(indsEnc2_load2shown));    % order of images shown during Enc2, onlly those shown in Enc2    
        PicOrder_Enc1_shownOnly_inEnc2 = cell2mat(ID_Enc1(indsEnc2_load2shown));

        indsPref_Enc2 = find(PicOrder_Enc2_shownOnly_inEnc2==idMaxHz);

        indsNonPref_Enc2_NonPrefEnc1 = find(PicOrder_Enc2_shownOnly_inEnc2 ~=idMaxHz &  PicOrder_Enc1_shownOnly_inEnc2 ~= idMaxHz );   % pref not shown in Enc2, and not shown in Enc1
        indsNonPref_Enc2_PrefEnc1 = find(PicOrder_Enc2_shownOnly_inEnc2    ~=idMaxHz &   PicOrder_Enc1_shownOnly_inEnc2 == idMaxHz);   % pref not shown in Enc2, but was shown in Enc1
        params_Enc2 = spikePlotParams; params_Enc2.colors = params_Enc2.colors([1 3 2]);
        [axRaster_Enc2,axPSTH_Enc2] = plotSpikeRaster_multiple_simplified( 3, offsets, limitRange_forRaster, limitRange_forPSTH, markerPos, ...
            SU.spike_times, {periods.Enc2}, {indsPref_Enc2, indsNonPref_Enc2_PrefEnc1, indsNonPref_Enc2_NonPrefEnc1 }, {'P_{Enc2 Only}', 'NP_{Enc2} P_{Enc1}','NP_{Enc2} NP_{Enc1}'}, [2 6 8], [2 6 2], params_Enc2 );
        title(['Enc2'])
        
  
        %% === Encoding 3 == Plot Raster and PSTH
        offsets = [0];
        limitRange_forRaster = [ 0 StimBaselineEnc+StimOnTimeEnc+ 2*padding  ];
        limitRange_forPSTH = [ 0 StimBaselineEnc+StimOnTimeEnc+ 2*padding  ]; 
        markerPos=[StimBaselineEnc+ padding  StimBaselineEnc+StimOnTimeEnc+ padding ]; %[StimBaselineEnc 1.5];
        
        indsEnc3_load3shown = find([ID_Enc3{:}] > 0);
        
        PicOrder_Enc3_shownOnly_inEnc3 = cell2mat(ID_Enc3(indsEnc3_load3shown));
        PicOrder_Enc2_shownOnly_inEnc3 = cell2mat(ID_Enc2(indsEnc3_load3shown));
        PicOrder_Enc1_shownOnly_inEnc3 = cell2mat(ID_Enc1(indsEnc3_load3shown));
        
        indsPref_Enc3 = find(PicOrder_Enc3_shownOnly_inEnc3==idMaxHz);
        
        indsNonPref_Enc3_NonPrefEnc12 = find(PicOrder_Enc3_shownOnly_inEnc3~=idMaxHz &  PicOrder_Enc2_shownOnly_inEnc3 ~= idMaxHz &  PicOrder_Enc1_shownOnly_inEnc3 ~= idMaxHz ); % pref not shown in Enc3, not in Enc2, not in Enc1
        
        indsNonPref_Enc3_PrefEnc2 = find(PicOrder_Enc3_shownOnly_inEnc3~=idMaxHz &  PicOrder_Enc2_shownOnly_inEnc3 == idMaxHz ); % pref not shown in Enc3, pref shown in Enc2
        indsNonPref_Enc3_PrefEnc1 = find(PicOrder_Enc3_shownOnly_inEnc3~=idMaxHz &  PicOrder_Enc1_shownOnly_inEnc3 == idMaxHz ); % pref not shown in Enc3, pref shown in Enc1
        params_Enc3 = spikePlotParams; params_Enc3.colors = params_Enc3.colors([1 3 2]);
        [axRaster_Enc3,axPSTH_Enc3] = plotSpikeRaster_multiple_simplified( 3, offsets, limitRange_forRaster, limitRange_forPSTH, markerPos, ...
            SU.spike_times, {periods.Enc3}, {indsPref_Enc3,[indsNonPref_Enc3_PrefEnc2; indsNonPref_Enc3_PrefEnc1],indsNonPref_Enc3_NonPrefEnc12 }, {'P','NP P_{E2} & NP P_{E1}','NP NP'}, [2 6 9], [2 6 3], params_Enc3 );
        title(['Enc3 ' ])


        %% === Maintenance == Plot Raster and PSTH
        offsets = [0];
        limitRange_forRaster = [0 StimBaselineMaint + StimOnTimeMaint + 2*padding ];
        limitRange_forPSTH = [0 StimBaselineMaint + StimOnTimeMaint + 2*padding ]; %
        markerPos=[StimBaselineMaint + padding, ...
            StimBaselineMaint + 1 + padding, ...
            StimBaselineMaint + 2 + padding];
            % StimBaselineMaint + StimOnTimeMaint + padding];% 
        
        
        % Temp In-memory all
        encFilter = (idMaxHz == [ID_Enc1{:}]) + (idMaxHz == [ID_Enc2{:}]) + (idMaxHz == [ID_Enc3{:}]) > 0;
        indsEnc = find(encFilter);
        indsNonEnc = find(~encFilter);

        params_Maint = spikePlotParams; params_Maint.colors = params_Maint.colors([3 2]);
        [axRaster_Maint,axPSTH_Maint] = plotSpikeRaster_multiple_simplified( 3, offsets, limitRange_forRaster, limitRange_forPSTH, markerPos, ...
            SU.spike_times, {periods.Maint}, {indsEnc,indsNonEnc}, {'P_{Enc1-3}','NP'},   [2 6 10], [2 6 4], params_Maint );

        % [axRaster_Maint,axPSTH_Maint] = plotSpikeRaster_multiple_simplified( 3, offsets, limitRange_forRaster, limitRange_forPSTH, markerPos, ...
        %     SU.spike_times, {periods.Maint}, {indsNonPref,indsPref}, {'NP_{Enc1}','P_{Enc1}'}, [2 6 4],  [2 6 10],spikePlotParams );
        title(['Maint'])
        % ylim([0 20])

        %% === Probe == Plot Raster and PSTH
        ProbeInOut = nwbAll{SU.session_count}.intervals_trials.vectordata.get('probe_in_out').data.load();
        probeID_mat = cell2mat(ID_Probe);
        indsProbe_pref_inmem = find( probeID_mat == idMaxHz & ProbeInOut==1 );
        indsProbe_pref_outmem = find( probeID_mat == idMaxHz & ProbeInOut==0 );
        indsProbe_nonpref_inmem = find( probeID_mat ~= idMaxHz & ProbeInOut==1 );
        indsProbe_nonpref_outmem = find( probeID_mat ~= idMaxHz & ProbeInOut==0 );
        
        offsets = 0;
        limitRange_forRaster = [ 0 StimBaselineProbe + StimOnTimeProbe + padding];
        limitRange_forPSTH = [ 0 StimBaselineProbe  + StimOnTimeProbe + padding]; %
        markerPos=[StimBaselineProbe + padding, ...
            StimBaselineProbe + 1 + padding, ...
            StimBaselineProbe + 2 + padding];

        params_probe = spikePlotParams; params_probe.colors = params_probe.colors([4 1 3 2]);

        [axRaster_Probe,axPSTH_Probe] = plotSpikeRaster_multiple_simplified( 3, offsets, limitRange_forRaster, limitRange_forPSTH, markerPos, ...
            SU.spike_times, {periods.Probe}, {indsProbe_pref_inmem, indsProbe_pref_outmem, indsProbe_nonpref_inmem,indsProbe_nonpref_outmem}, ...
            { 'P in', 'P out','NP in','NP out'},   [2 6 11], [2 6 5], params_probe ); 
        title(['Probe ']);
        

        %% === Selective Image ===    
        indPrefImg = idUnique(find(idUnique==idMaxHz,1)); % Index in StimulusTemplates at which the preferred image occurs
        current_session = all_units(i).session_count;
        StimulusTemplates = nwbAll{current_session}.stimulus_templates.get('StimulusTemplates') ;
       
        prefPathUnstripped = StimulusTemplates.order_of_images.data(indPrefImg).path; % Loading raw path
        stripPrepPath = split(prefPathUnstripped,'/'); % Stripping fileseps
        prefPath = stripPrepPath{end}; % Grabbing final entry
        
        prefImgload = StimulusTemplates.image.get(prefPath).data.load; % Loading using the key from order_of_images. 
        prefImg = permute(prefImgload,[2,3,1]); % Permuting for matlab compatibility
        
        subplot(2,6,6)
        image(prefImg);
        pbaspect([4 3 1])
        set(gca,'YTick',[])
        set(gca,'XTick',[])
        xlabel(strrep(prefPath,'_',' '))
        
        %% Plotting Spike pdf
        subplot(2,6,12)
        waveforms = SU.waveforms;
        if size(waveforms,1) > 1 % Plot spike pdf
            [~,~] = spikePDFestimate(waveforms);
        else % Only plot mean waveform if that's the only waveform available. 
            plot(waveforms)
            xlim([0 length(waveforms)])
            ylim([min(waveforms)*1.2 max(waveforms)*1.2])
        end
        ylabel('\muV')
        set(gca,'XTick',[])
        pbaspect([4 3 1])
        
        set(gca, 'XTick',0:100:256)
        set(gca, 'XTickLabel',0:1:256/100)
        xlabel('Time (ms)')
        ylabel('\muV')
        
        %% === Formatting Plot ===
        setCommonAxisRange( [axPSTH_Enc1 axPSTH_Enc2 axPSTH_Enc3 ,axPSTH_Maint,axPSTH_Probe], 2 );   % set common Y axis
        setCommonAxisRange( [axRaster_Enc1 axRaster_Enc2 axRaster_Enc3 ,axRaster_Maint,axRaster_Probe], 2 );   % set common Y axis
        fig = gcf;
        fig.OuterPosition = [0.25 0.25 1920 1080];
        fig.WindowState = 'maximized';
        movegui(gcf,'west')
        fig.WindowState = 'maximized';   

        axPSTH_Enc1.Legend.Position(2) = axPSTH_Enc1.Legend.Position(2) - 0.10;
        axPSTH_Enc2.Legend.Position(2) = axPSTH_Enc2.Legend.Position(2) - 0.10;
        axPSTH_Enc3.Legend.Position(2) = axPSTH_Enc3.Legend.Position(2) - 0.10;
        axPSTH_Maint.Legend.Position(2) = axPSTH_Maint.Legend.Position(2) - 0.10;
        axPSTH_Probe.Legend.Position(2) = axPSTH_Probe.Legend.Position(2) - 0.10;


        %% === Export ===
        if params.exportFig
            % set(gcf,'position',[-1920 817 1920 175])
            if ~isfield(params,'figOut') || isempty(params.figOut)
            figPath = 'C:\temp\figsSternberg\test_output'; 
            else
                figPath = params.figOut;
            end
            
            % Append selectivity to folder. Will store in separate folder
            % for each type. Allows for folders for overlapping cell types.
            if concept_cells_sb(i)
                figPath = [figPath '_concept']; %#ok<AGROW>
            end
            if maint_cells_sb(i)
                figPath = [figPath '_maint']; %#ok<AGROW>
            end
            if probe_cells_sb(i)
                figPath = [figPath '_probe']; %#ok<AGROW>
            end

            if ~isfolder(figPath)
                mkdir(figPath)
            end
            fName = ['SBID_' num2str(subject_id) '_Cell_' num2str(cellID) '_Cl_' num2str(clusterID) '_BA_' brain_area{:}];
            saveas(gcf, [figPath filesep fName '.' params.exportType ],params.exportType)
            fprintf('| Figure Saved')
            close(gcf);
        end
    end
    fprintf('\n') % New line for each cell
end
% Comparing concept cells for preferred and non-preferred trials. 
hzPref_conceptOnly = hzPref(logical(concept_cells_sb));
hzNonPref_conceptOnly = hzNonPref(logical(concept_cells_sb));
[~, p_prefNonPref] = ttest(hzPref_conceptOnly,hzNonPref_conceptOnly,'Tail','right','Alpha',alphaLim);
% p_prefNonPref = permutationTest(hzPref_conceptOnly,hzNonPref_conceptOnly,2000,'sidedness','larger');

fprintf('Total Concept: %d/%d (%.2f%%)\n',sum(concept_cells_sb),length(all_units),sum(concept_cells_sb)/length(all_units)*100)
fprintf('Persistent Firing : Pref %.2f ± %.2f | Non-Pref %.2f ± %.2f | p = %.4f\n',mean(hzPref_conceptOnly),std(hzPref_conceptOnly),mean(hzNonPref_conceptOnly),std(hzNonPref_conceptOnly), p_prefNonPref )
fprintf('Total Maint: %d/%d (%.2f%%)\n',sum(maint_cells_sb),length(all_units),sum(maint_cells_sb)/length(all_units)*100)
fprintf('Total Probe: %d/%d (%.2f%%)\n',sum(probe_cells_sb),length(all_units),sum(probe_cells_sb)/length(all_units)*100)
sig_cells.concept_cells = concept_cells_sb;
sig_cells.hzPref = hzPref;
sig_cells.hzNonPref = hzNonPref;
sig_cells.maint_cells = maint_cells_sb;
sig_cells.probe_cells = probe_cells_sb;
end

