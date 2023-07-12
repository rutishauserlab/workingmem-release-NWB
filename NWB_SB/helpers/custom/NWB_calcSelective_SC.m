function [sig_cells_sc,areasScreening] = NWB_calcSelective_SC(nwbAll,all_units, params)
%NWB_CALCSELECTIVE_SC Takes the output of NWB_SB_extractUnits and runs
%selectivity tests for the screening task. Optionally does plotting of
%response histograms. 

if isfield(params,'rateFilter')
    rateFilter = params.rateFilter;
else
    rateFilter = [];
end

% Filtering for Global Rate (rateFilter should be a nonzero float. Otherwise, all cells are used.)
aboveRate = ones(length(all_units),1);
if ~isempty(rateFilter)
    for i = 1:length(all_units)
        globalRate = length(all_units(i).spike_times)/(max(all_units(i).spike_times)-min(all_units(i).spike_times));
        rateBool = globalRate < rateFilter ;
        if rateBool % If the rate is below the filter threshold
            aboveRate(i) = 0;
        end
    end
end
all_units = all_units(logical(aboveRate));


areasScreening = cell(length(all_units),1);
sig_cells_sc = zeros(length(all_units),1); alphaLim = 0.05;
if isfield(params,'runParallel') && params.runParallel == 1 % Setting to serial or parallel processing. 
    parforArg = Inf; % Note: status messages will not occur in order
else
    parforArg = 0;
end
for i = 1:length(all_units)%, parforArg) % Parallel currently disabled. 
    SU = all_units(i);
    subject_id = SU.subject_id;
    clusterID = nwbAll{SU.session_count}.units.vectordata.get('clusterID_orig').data.load(SU.unit_id);
    brain_area = nwbAll{SU.session_count}.general_extracellular_ephys_electrodes.vectordata.get('location').data.load(SU.electrodes);
    areasScreening{i} = brain_area{:};
    fprintf('Processing: (%d/%d) Session SCID %d, Unit %d, Cluster %d ',i,length(all_units),SU.subject_id,SU.unit_id,clusterID)
    
    % Load Trial timestamps & Stimulus IDs
    tsIn = num2cell(nwbAll{SU.session_count}.intervals_trials.start_time.data.load());
    idIn = num2cell(nwbAll{SU.session_count}.stimulus_presentation.get('StimulusPresentation').data.load() + 1); % Offset to 1-indexing
    stim = cell2struct(horzcat(tsIn,idIn), {'ts','id'}, 2);
    
    [idUnique, ~, ic] = unique([stim.id]);
    % uniqueCounts = histcounts(ic,'BinMethod','integers');

    HzTrial = NaN(length(stim),1);
    signalDelay = 0.2; % Delay of stimulus onset to effect. 
    stimOffset = 1; % Time past stimulus onset. End of picture presentation.
    for k = 1:length(HzTrial)
        periodFilter = (SU.spike_times>(stim(k).ts+signalDelay)) & (SU.spike_times<(stim(k).ts+stimOffset));
        singleTrialSpikes = SU.spike_times(periodFilter);
        trialRate = length(singleTrialSpikes)/(stimOffset-signalDelay); % Firing rate across testing period.
        HzTrial(k) = trialRate;
    end

    % Finding image with maximum response. Used in stage 2 test. 
    nUniqueStim = length(idUnique);
    mResp = nan(nUniqueStim,1);
    stdResp = nan(nUniqueStim,1);
    semResp = nan(nUniqueStim,1);
    for k=1:nUniqueStim
        stimRates =HzTrial([stim.id]==idUnique(k));
        mResp(k) = mean(stimRates);
        stdResp(k) = std(stimRates); % Might not be needed
        semResp(k) = stdResp(k)/sqrt(length(stimRates));
    end
    
    idMaxHz = idUnique(find(mResp==max(mResp),1)); % Find the first entry where the max value is reached. 
    id_Trial_maxOnly = [stim.id]; id_Trial_maxOnly([stim.id]~=idMaxHz) = -1;

    % Significance Tests
    p_ANOVA = 1; p_bANOVA = 1; %#ok<NASGU> % Preset as '1' to allow for paramsSB.plotAlways.
    p_ANOVA = anovan(HzTrial,{string([stim.id])}, 'display','off','model', 'linear','alpha',alphaLim,'varnames','picID');
    if p_ANOVA < alphaLim % First test: 1-way ANOVA
        p_bANOVA = anovan(HzTrial,{string(id_Trial_maxOnly)}, 'display','off','model', 'linear','alpha',alphaLim,'varnames','picID');
        if p_bANOVA < alphaLim % Second Test: Binarized 1-way ANOVA (simplifies to t-test)
        fprintf('| SIGNIFICANT: Session SBID %d, Unit %d ',SU.subject_id,SU.unit_id)
        sig_cells_sc(i) = 1 ;
        end
    end

    % % If we wanted the permuted versions mentioned in K2017:
    % p_permANOVA = randanova1(HzTrial,string([stim.id]),2000);
    % a = HzTrial([stim.id]~=idMaxHz); b = HzTrial([stim.id]==idMaxHz);
    % p_permTest = permutationTest(a,b,2000);

    % Plotting Histogram & Selective Image
    isTuned = sig_cells_sc(i);
    if (isTuned && params.doPlot) || params.plotAlways
        pTitleStr = ['sub-' num2str(subject_id) '-ses-1'...
            ' | BA: ' strrep(brain_area{:},'_',' ') ...
            ' | Cell: ' num2str(SU.unit_id) ...
            ' | Elec: ' num2str(SU.electrodes) ...
            ' | p1=' sprintf('%.2f',p_ANOVA) ...
            ' | p2=' sprintf('%.2f',p_bANOVA)];
        % Bar Plot w/ Error Bars
        figure('Name',pTitleStr)
        subplot(1,7,1:5)
        collapseID = 1; % Whether to collapse the histogram to not consider pic ID in the x axis, rather the picture itself. Used for debugging
        if collapseID
            idPlot = 1:length(idUnique); %#ok<UNRCH>
        else
            idPlot = idUnique; %#ok<UNRCH>
        end
        bar_h = bar(idPlot,mResp);
        
        hold on
        er = errorbar(idPlot,mResp,semResp,semResp);
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';  
        
        hold off
        xlabel('Picture ID')
        ylabel('Firing Rate (Hz)')
        title(pTitleStr)
        if collapseID
            tickLabels = num2cell(idUnique); tickLabels = cellfun(@(x) num2str(x),tickLabels,'UniformOutput',false);
            xticks(idPlot)
            xticklabels(tickLabels)
            h=gca; h.XAxis.TickLength = [0 0];
        end

        % Getting Selective Image
        idPrefImg = idMaxHz; % Index in StimulusTemplates at which the preferred image occurs. 
        
        % Adding markers to max firing rate entry
        
        ast_height_factor = 1.05;
        find_max_idPlot = idPlot(find(mResp==max(mResp),1));
        ast_height = (mResp(find_max_idPlot) + semResp(find_max_idPlot))*ast_height_factor;
        if isTuned
            % Adding asterisk
            text(idPlot(find_max_idPlot),ast_height,'\bf*','FontSize',13);
        end
        % Color selective bar
        bar_h.FaceColor = 'flat';
        bar_h.CData(idPlot(find_max_idPlot),:) = [1 0 0];
        ylim([0 ast_height*1.1])

        current_session = all_units(i).session_count;
        StimulusTemplates = nwbAll{current_session}.stimulus_templates.get('StimulusTemplates') ;
       

        prefPath = ['image_' num2str(idPrefImg)];
        prefImgload = StimulusTemplates.image.get(prefPath).data.load; % Loading using the key from order_of_images. 
        prefImg = permute(prefImgload,[2,3,1]); % Permuting for matlab compatibility
        
        % Plotting preferred image
        subplot(1,7,6)
        image(prefImg);
        xlabel(strrep(prefPath,'_',' '))
        
        % pbaspect([4 3 1])
        set(gca,'YTick',[])
        set(gca,'XTick',[])
        
        % % Plotting waveforms PDF. 
        subplot(1,7,7)
        waveforms = SU.waveforms;
        if size(waveforms,1) > 1 % Plot spike pdf
            [~,~] = spikePDFestimate(waveforms);
        else % Only plot mean waveform if that's the only waveform available. 
            plot(waveforms)
            xlim([0 length(waveforms)])
            ylim([min(waveforms)*1.2 max(waveforms)*1.2])
        end

        % Setting x axis to ms (Hardcoded to 100kHz sample rate, 256 samples)
        set(gca, 'XTick',0:100:256)
        set(gca, 'XTickLabel',0:1:256/100)
        xlabel('Time (ms)')
        ylabel('\muV')

        % set(gcf,'position',[-1920 800 1920 180]) 
         
        shg
        if params.exportFig
            set(gcf,'position',[-1920 800 1920 180])
            if ~isfield(params,'figOut') || isempty(params.figOut)
            figPath = 'C:\temp\figsSternberg\'; 
            else
                figPath = params.figOut;
            end
            if ~isfolder(figPath)
                mkdir(figPath)
            end
            subject_id = SU.subject_id;
            cellID = SU.unit_id;
            
            clusterID = nwbAll{SU.session_count}.units.vectordata.get('clusterID_orig').data.load(SU.unit_id);
            fName = ['SCID_' num2str(subject_id) '_Cell_' num2str(cellID) '_Cl_' num2str(clusterID) '_BA_' brain_area{:}];
            saveas(gcf, [params.figOut filesep fName '.' params.exportType ], params.exportType)
            fprintf('| Figure Saved')
            close(gcf);
        end
        
    end
    fprintf('\n') % New line for each cell
end
    fprintf('Total Significant:%d/%d (%.2f%%)\n',sum(sig_cells_sc),length(all_units),sum(sig_cells_sc)/length(all_units)*100)
end



