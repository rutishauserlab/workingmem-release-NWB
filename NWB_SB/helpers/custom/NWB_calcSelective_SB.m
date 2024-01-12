function [sig_cells,areasSternberg] = NWB_calcSelective_SB(nwbAll, all_units, params)
%NWB_CALCSELECTIVE Takes the output of NWB_SB_extractUnits and runs
%selectivity tests for the sternberg task. Optionally plots trial aligned
%spikes

if isfield(params,'rateFilter') && ~isempty(params.rateFilter) && params.rateFilter > 0
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
concept_cells_sb = zeros(length(all_units),1); 
hzPref = zeros(length(all_units),1);
hzNonPref = zeros(length(all_units),1);
maint_cells_sb = zeros(length(all_units),1);
probe_cells_sb = zeros(length(all_units),1);
% Looping over all cells
for i = 1:length(all_units) 
    SU = all_units(i);
    subject_id = SU.subject_id;
    session_id = SU.session_id;
    identifier = SU.identifier;
    cellID = SU.unit_id;
    brain_area = nwbAll{SU.session_count}.general_extracellular_ephys_electrodes.vectordata.get('location').data.load(SU.electrodes);
    clusterID = nwbAll{SU.session_count}.units.vectordata.get('clusterID_orig').data.load(SU.unit_id);
    areasSternberg{i} = brain_area{:};
    fprintf('Processing: (%d/%d) sub-%s-ses-%s, Unit %d Cluster %d ',i,length(all_units),subject_id,session_id,cellID,clusterID)
    
    % Loading stim timestamps and loads
    tsFix = num2cell(nwbAll{SU.session_count}.intervals_trials.vectordata.get('timestamps_FixationCross').data.load());
    tsEnc1 = num2cell(nwbAll{SU.session_count}.intervals_trials.vectordata.get('timestamps_Encoding1').data.load());
    tsEnc2 = num2cell(nwbAll{SU.session_count}.intervals_trials.vectordata.get('timestamps_Encoding2').data.load());
    tsEnc3 = num2cell(nwbAll{SU.session_count}.intervals_trials.vectordata.get('timestamps_Encoding3').data.load());
    tsMaint = num2cell(nwbAll{SU.session_count}.intervals_trials.vectordata.get('timestamps_Maintenance').data.load());
    tsProbe = num2cell(nwbAll{SU.session_count}.intervals_trials.vectordata.get('timestamps_Probe').data.load());
    ID_Enc1 = num2cell(nwbAll{SU.session_count}.intervals_trials.vectordata.get('PicIDs_Encoding1').data.load());
    ID_Enc2 = num2cell(nwbAll{SU.session_count}.intervals_trials.vectordata.get('PicIDs_Encoding2').data.load());
    ID_Enc3 = num2cell(nwbAll{SU.session_count}.intervals_trials.vectordata.get('PicIDs_Encoding3').data.load());
    ID_Probe = num2cell(nwbAll{SU.session_count}.intervals_trials.vectordata.get('PicIDs_Probe').data.load());
    
    % Deriving CAT classifers (should be in range of 1-5 and 9
    CAT_Enc1 = cellfun(@(x) num2str(x),ID_Enc1,'UniformOutput',false); CAT_Enc1 = cellfun(@(x) str2double(x(1)),CAT_Enc1,'UniformOutput',false);
    CAT_Enc2 = cellfun(@(x) num2str(x),ID_Enc2,'UniformOutput',false); CAT_Enc2 = cellfun(@(x) str2double(x(1)),CAT_Enc2,'UniformOutput',false);
    CAT_Enc3 = cellfun(@(x) num2str(x),ID_Enc3,'UniformOutput',false); CAT_Enc3 = cellfun(@(x) str2double(x(1)),CAT_Enc3,'UniformOutput',false);
    CAT_Probe = cellfun(@(x) num2str(x),ID_Probe,'UniformOutput',false); CAT_Probe = cellfun(@(x) str2double(x(1)),CAT_Probe,'UniformOutput',false);

    stim = cell2struct(...
        horzcat(tsFix,tsEnc1,tsEnc2,tsEnc3,tsMaint, tsProbe,...
        ID_Enc1,ID_Enc2,ID_Enc3,ID_Probe,CAT_Enc1,CAT_Enc2,CAT_Enc3,CAT_Probe),...
        {'tsFix','tsEnc1','tsEnc2','tsEnc3','tsMaint','tsProbe',...
        'idEnc1','idEnc2','idEnc3','idProbe','CAT_Enc1','CAT_Enc2','CAT_Enc3','CAT_Probe'},2);
    
    [idUnique, ~, ic] = unique([stim.CAT_Enc1]);    
    uniqueCounts = histcounts(ic,'BinMethod','integers');
    
    %% Get all stimulus rates
    signalDelay = 0.2; % Delay of stimulus onset to effect. 
    stimOffset = 1; % Time past stimulus onset. End of picture presentation.
    HzEnc1 = NaN(length(stim),1);
    for k = 1:length(HzEnc1)
        periodFilter = (SU.spike_times>(stim(k).tsEnc1+signalDelay)) & (SU.spike_times<(stim(k).tsEnc1+stimOffset));
        singleTrialSpikes = SU.spike_times(periodFilter);
        trialRate = length(singleTrialSpikes)/(stimOffset-signalDelay); % Firing rate across testing period.
        HzEnc1(k) = trialRate;
    end
    nonZ_tsEnc2 = nonzeros([stim.tsEnc2]); nonZ_CAT_Enc2 = nonzeros([stim.CAT_Enc2]);
    HzEnc2 = NaN(length(nonZ_tsEnc2),1);
    for k = 1:length(HzEnc2)
        periodFilter = (SU.spike_times>(nonZ_tsEnc2(k)+signalDelay)) & (SU.spike_times<(nonZ_tsEnc2(k)+stimOffset));
        singleTrialSpikes = SU.spike_times(periodFilter);
        trialRate = length(singleTrialSpikes)/(stimOffset-signalDelay); % Firing rate across testing period.
        HzEnc2(k) = trialRate;
    end
    nonZ_tsEnc3 = nonzeros([stim.tsEnc3]); nonZ_CAT_Enc3 = nonzeros([stim.CAT_Enc3]);
    HzEnc3 = NaN(length(nonZ_tsEnc3),1);
    for k = 1:length(HzEnc3)
        periodFilter = (SU.spike_times>(nonZ_tsEnc3(k)+signalDelay)) & (SU.spike_times<(nonZ_tsEnc3(k)+stimOffset));
        singleTrialSpikes = SU.spike_times(periodFilter);
        trialRate = length(singleTrialSpikes)/(stimOffset-signalDelay); % Firing rate across testing period.
        HzEnc3(k) = trialRate;
    end
    nonZ_tsProbe = nonzeros([stim.tsProbe]); nonZ_CAT_Probe = nonzeros([stim.CAT_Probe]);
    HzProbe = NaN(length(nonZ_tsProbe),1);
    for k = 1:length(HzProbe)
        periodFilter = (SU.spike_times>(nonZ_tsProbe(k)+signalDelay)) & (SU.spike_times<(nonZ_tsProbe(k)+stimOffset));
        singleTrialSpikes = SU.spike_times(periodFilter);
        trialRate = length(singleTrialSpikes)/(stimOffset-signalDelay); % Firing rate across testing period.
        HzProbe(k) = trialRate;
    end

    %% Compile Stimulus Rates
    Hz_allTrials = [HzEnc1;HzEnc2;HzEnc3;HzProbe];
    CAT_allTrials = [[stim.CAT_Enc1]';nonZ_CAT_Enc2;nonZ_CAT_Enc3;nonZ_CAT_Probe];

    % Finding image with maximum response. Used in stage 2 test. 
    nUniqueStim = length(unique(CAT_allTrials));
    mResp=nan(nUniqueStim,1);
    for k=1:nUniqueStim
        mResp(k) = mean(Hz_allTrials(CAT_allTrials==k));
    end
    idMaxHz = find(mResp==max(mResp),1);
    id_Trial_maxOnly = CAT_allTrials; id_Trial_maxOnly(CAT_allTrials~=idMaxHz) = 0;
    
    %% Significance Tests: Category Cells
    % Count spikes in 200-1000ms window following all stimuli onsets (Enc1,
    % Enc2, Enc3, & Probe)
    % following the first encoding period. Use a 1-way ANOVA followed by a
    % t-test of the maximal response versus the non-selective responses.
    % Note that using a 1-way anova with two groups simplifies to a t-test
    alphaLim = 0.05;
    % sig_method = 'parametric';
    sig_method = 'non-parametric';
    if strcmp(sig_method,'parametric')
        p_ANOVA = 1; p_bANOVA = 1;%#ok<NASGU> % Preset as '1' to allow for paramsSB.plotAlways.
        p_ANOVA = anovan(Hz_allTrials,{string(CAT_allTrials)}, 'display','off','model', 'linear','alpha',alphaLim,'varnames','picID');
        if p_ANOVA < alphaLim % First test: 1-way ANOVA
            p_bANOVA = anovan(Hz_allTrials,{string(id_Trial_maxOnly)}, 'display','off','model', 'linear','alpha',alphaLim,'varnames','picID');
            if p_bANOVA < alphaLim % Second Test: Binarized 1-way ANOVA (simplifies to t-test)
            fprintf('| Category -> sub-%s-ses-%s, Unit %d p1:%.2f p2:%.2f',SU.subject_id,SU.session_id,SU.unit_id,p_ANOVA,p_bANOVA)
            concept_cells_sb(i) = 1;
            end
        end
    elseif strcmp(sig_method,'non-parametric')
        p_ANOVA = 1; p_bANOVA = 1;%#ok<NASGU> % Preset as '1' to allow for paramsSB.plotAlways.
        groups = cellstr(string(CAT_allTrials));
        p_ANOVA = randanova1(Hz_allTrials,groups, 2000);
        if p_ANOVA < alphaLim % First test: 1-way ANOVA            
            a = Hz_allTrials(CAT_allTrials==idMaxHz); b = Hz_allTrials(CAT_allTrials~=idMaxHz); 
            p_perm = permutationTest(a,b,2000,'sidedness','larger');
            if p_perm < alphaLim % Second Test: Binarized 1-way ANOVA (simplifies to t-test)
            fprintf('| Category -> sub-%s-ses-%s, Unit %d p1:%.2f p2:%.2f',SU.subject_id,SU.session_id,SU.unit_id,p_ANOVA,p_perm)
            concept_cells_sb(i) = 1;
            end
            p_bANOVA = p_perm; % For var references later in the script.
        end      
    else
        error('Significance method not specified')
    end


    % Saving rate data for additional tests between pref/non-pref trials
    hzPref(i) = mean(Hz_allTrials(logical(id_Trial_maxOnly)));
    hzNonPref(i) = mean(Hz_allTrials(~logical(id_Trial_maxOnly)));


    %% Rasters, PSTH, & Selective Image
    % Flagging significant cells for plotting
    switch params.plotMode
        case 1 % Category
            plotFlag = params.doPlot && concept_cells_sb(i);
        case 2 % Maint
            plotFlag = params.doPlot && maint_cells_sb(i);
        case 3 % Probe
            plotFlag = params.doPlot && probe_cells_sb(i);
        case 4 % All
            plotFlag = params.doPlot && (concept_cells_sb(i) || maint_cells_sb(i) || probe_cells_sb(i));
        otherwise
            warning('Plot mode not specified. Defaulting to category cells.\n')
            in = input('Continue? (y/n)\n',"s");
            if any(strcmp(in,["Y","y"]))
                plotFlag = params.doPlot && concept_cells_sb(i);
            elseif any(strcmp(in,["N","n"]))
                fprintf('Aborting.\n')
                return
            else
                fprintf('Answer not specified. Aborting.\n')
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
        StimOnTimeEnc  = 2; % after for all encodings
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
        indsPref = [stim.CAT_Enc1]==idMaxHz; indsNonPref = ~indsPref;
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
        PicOrder_Enc2_shownOnly_inEnc2 = cell2mat(CAT_Enc2(indsEnc2_load2shown));    % order of images shown during Enc2, onlly those shown in Enc2    
        PicOrder_Enc1_shownOnly_inEnc2 = cell2mat(CAT_Enc1(indsEnc2_load2shown));

        indsPref_Enc2 = find(PicOrder_Enc2_shownOnly_inEnc2==idMaxHz);

        indsNonPref_Enc2_NonPrefEnc1 = find(PicOrder_Enc2_shownOnly_inEnc2 ~=idMaxHz &  PicOrder_Enc1_shownOnly_inEnc2 ~= idMaxHz );   % pref not shown in Enc2, and not shown in Enc1
        indsNonPref_Enc2_PrefEnc1 = find(PicOrder_Enc2_shownOnly_inEnc2    ~=idMaxHz &   PicOrder_Enc1_shownOnly_inEnc2 == idMaxHz);   % pref not shown in Enc2, but was shown in Enc1
        params_Enc2 = spikePlotParams; params_Enc2.colors = params_Enc2.colors([1 2]);
        [axRaster_Enc2,axPSTH_Enc2] = plotSpikeRaster_multiple_simplified( 3, offsets, limitRange_forRaster, limitRange_forPSTH, markerPos, ...
            SU.spike_times, {periods.Enc2}, {indsPref_Enc2, sort([indsNonPref_Enc2_PrefEnc1;indsNonPref_Enc2_NonPrefEnc1]) }, {'Pref', 'NonPref'}, [2 6 8], [2 6 2], params_Enc2 );
        title(['Enc2'])
        
  
        %% === Encoding 3 == Plot Raster and PSTH
        offsets = [0];
        limitRange_forRaster = [ 0 StimBaselineEnc+StimOnTimeEnc+ 2*padding  ];
        limitRange_forPSTH = [ 0 StimBaselineEnc+StimOnTimeEnc+ 2*padding  ]; 
        markerPos=[StimBaselineEnc+ padding  StimBaselineEnc+StimOnTimeEnc+ padding ]; %[StimBaselineEnc 1.5];
        
        indsEnc3_load3shown = find([ID_Enc3{:}] > 0);
        
        PicOrder_Enc3_shownOnly_inEnc3 = cell2mat(CAT_Enc3(indsEnc3_load3shown));
        PicOrder_Enc2_shownOnly_inEnc3 = cell2mat(CAT_Enc2(indsEnc3_load3shown));
        PicOrder_Enc1_shownOnly_inEnc3 = cell2mat(CAT_Enc1(indsEnc3_load3shown));
        
        indsPref_Enc3 = find(PicOrder_Enc3_shownOnly_inEnc3==idMaxHz);
        
        indsNonPref_Enc3_NonPrefEnc12 = find(PicOrder_Enc3_shownOnly_inEnc3~=idMaxHz &  PicOrder_Enc2_shownOnly_inEnc3 ~= idMaxHz &  PicOrder_Enc1_shownOnly_inEnc3 ~= idMaxHz ); % pref not shown in Enc3, not in Enc2, not in Enc1
        
        indsNonPref_Enc3_PrefEnc2 = find(PicOrder_Enc3_shownOnly_inEnc3~=idMaxHz &  PicOrder_Enc2_shownOnly_inEnc3 == idMaxHz ); % pref not shown in Enc3, pref shown in Enc2
        indsNonPref_Enc3_PrefEnc1 = find(PicOrder_Enc3_shownOnly_inEnc3~=idMaxHz &  PicOrder_Enc1_shownOnly_inEnc3 == idMaxHz ); % pref not shown in Enc3, pref shown in Enc1
        params_Enc3 = spikePlotParams; params_Enc3.colors = params_Enc3.colors([1 2]);
        [axRaster_Enc3,axPSTH_Enc3] = plotSpikeRaster_multiple_simplified( 3, offsets, limitRange_forRaster, limitRange_forPSTH, markerPos, ...
            SU.spike_times, {periods.Enc3}, {indsPref_Enc3,[indsNonPref_Enc3_PrefEnc2; indsNonPref_Enc3_PrefEnc1;indsNonPref_Enc3_NonPrefEnc12] }, {'Pref','NonPref'}, [2 6 9], [2 6 3], params_Enc3 );
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
        encFilter = (idMaxHz == [CAT_Enc1{:}]) + (idMaxHz == [CAT_Enc2{:}]) + (idMaxHz == [CAT_Enc3{:}]) > 0;
        indsEnc = find(encFilter);
        indsNonEnc = find(~encFilter);

        params_Maint = spikePlotParams; params_Maint.colors = params_Maint.colors([1 2]);
        [axRaster_Maint,axPSTH_Maint] = plotSpikeRaster_multiple_simplified( 3, offsets, limitRange_forRaster, limitRange_forPSTH, markerPos, ...
            SU.spike_times, {periods.Maint}, {indsEnc,indsNonEnc}, {'Pref','NonPref'},   [2 6 10], [2 6 4], params_Maint );

        title(['Maint'])

        %% === Probe == Plot Raster and PSTH
        ProbeInOut = nwbAll{SU.session_count}.intervals_trials.vectordata.get('probe_in_out').data.load();
        probeID_mat = cell2mat(CAT_Probe);
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

        params_probe = spikePlotParams; params_probe.colors = params_probe.colors([1 2]);

        [axRaster_Probe,axPSTH_Probe] = plotSpikeRaster_multiple_simplified( 3, offsets, limitRange_forRaster, limitRange_forPSTH, markerPos, ...
            SU.spike_times, {periods.Probe}, {[indsProbe_pref_inmem;indsProbe_pref_outmem], [indsProbe_nonpref_inmem;indsProbe_nonpref_outmem]}, ...
            { 'Pref', 'NonPref'},   [2 6 11], [2 6 5], params_probe ); 
        title(['Probe ']);
        

        %% === Selective Image ===    
        indPrefImg = idUnique(find(idUnique==idMaxHz,1)); % Index in StimulusTemplates at which the preferred image occurs
        % NOTE: Category not explicitly stated. Dependent on task version. 
        catLabel = sprintf('Selective Category:\n%d',indPrefImg);

        subplot(2,6,6)
        text(0.5,0.5,catLabel); axis off
        pbaspect([4 3 1])
        set(gca,'YTick',[])
        set(gca,'XTick',[])
        % xlabel(strrep(prefPath,'_',' '))
        
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
        % setCommonAxisRange( [axRaster_Enc1 axRaster_Enc2 axRaster_Enc3 ,axRaster_Maint,axRaster_Probe], 2 );   % set common Y axis
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
                figPath = [figPath '_category']; %#ok<AGROW>
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
            fName = [identifier '_Cell_' num2str(cellID) '_Cl_' num2str(clusterID) '_BA_' brain_area{:}];
            saveas(gcf, [figPath filesep fName '.' params.exportType ],params.exportType)
            fprintf('| Figure Saved')
            close(gcf);
        end
    end
    fprintf('\n') % New line for each cell
end
% Comparing category cells for preferred and non-preferred trials. 
hzPref_conceptOnly = hzPref(logical(concept_cells_sb));
hzNonPref_conceptOnly = hzNonPref(logical(concept_cells_sb));
[~, p_prefNonPref] = ttest(hzPref_conceptOnly,hzNonPref_conceptOnly,'Tail','right','Alpha',alphaLim);

fprintf('Total Category Cells: %d/%d (%.2f%%)\n',sum(concept_cells_sb),length(all_units),sum(concept_cells_sb)/length(all_units)*100)
sig_cells.concept_cells = concept_cells_sb;
sig_cells.hzPref = hzPref;
sig_cells.hzNonPref = hzNonPref;
sig_cells.maint_cells = maint_cells_sb;
sig_cells.probe_cells = probe_cells_sb;
end

