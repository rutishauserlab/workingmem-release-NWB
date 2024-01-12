

%% Access data loaded from import script

nwb = nwbAll_sb{:};
tsMaint = nwb.intervals_trials.vectordata.get('timestamps_Maintenance').data.load();

electrode_labels = string(nwb.general_extracellular_ephys_electrodes.vectordata.get('location').data.load());
electrode_origChan = string(nwb.general_extracellular_ephys_electrodes.vectordata.get('origChannel').data.load());


%% Load trial accuracy data (need to randomly subsample load 1/3 trials later)
% Just load the provided data for now:
loadPath = "C:\svnwork\neuro1\code\events\Sternberg\NWB-SBCAT\NWB_SBCAT_import_demo\helpers\internal\PAC\P60cs_RH1_stratTrials.mat";
loadIn = load(loadPath);
stratTrials = loadIn.stratTrials;


%%
range = 49;
lfp_raw = nwb.acquisition.get('LFPs').data.load();
lfp_electrodes = int32(nwb.acquisition.get('LFPs').electrodes.data.load())+1; % Offset to 1-indexing
lfp_labels = electrode_labels(lfp_electrodes);
lfp_origChan = electrode_origChan(lfp_electrodes);

for i = range
    lfp_ID = i; % Just a filler ID for now. 
    lfp_chan = lfp_raw(lfp_ID,:);
    lfp_start_ts = nwb.acquisition.get('LFPs').starting_time;
    lfp_srate = nwb.acquisition.get('LFPs').starting_time_rate; % Should be 400
    lfp_ts = lfp_start_ts + (0:length(lfp_chan))*(1/lfp_srate); % tspan with 1/hz increments, used in ts filter
    
    
    %% [sample X trial] array
    trial_pre = 0.5;
    trial_post = 3;
    
    lfp_trials = []; % Acquire LFP segments per trial
    for k=1:length(tsMaint)
        in_period = lfp_ts >= (tsMaint(k)-trial_pre) & lfp_ts <= (tsMaint(k)+trial_post);
        period_lfp = lfp_chan(in_period)';
        lfp_trials = [lfp_trials,period_lfp]; %#ok<AGROW> % It's just 140 trials. 
    end
    
    % Separate into loads 1 and 3
    lfp_trials_L1 = lfp_trials(:,stratTrials(:,1));
    lfp_trials_L3 = lfp_trials(:,stratTrials(:,2));
    
    
    %% Call section
    whichSession = 'P60cs';
    whichChannel = 'RH1'; 
    
    % load data from session and channel noted above
    
    % select time window of -.5 - 3s around delay period onset from all trials
    
    % CFC should be computed per load condition using the settings listed below
    
    % load stratTrials to know which trials to compute from each condition
    % (column 1 = correct Load 1 trials; column 2 = correct Load 3 trials)
    % these trials have been randomly subsampled so that the same number of
    % trials is used in both conditions (after excluding artifact trials)
    
    % generate a datMat variable that contains the delay periods from all
    % considered trials (see help of function below)
    
    % frequencies of interest for later phase2power plot LF = [4 6]; HF = [70 140]
    
    %%
    % datMat: samples x trials
    srate = lfp_srate; % Should be 400
    n_surrogates = 200; % 200 by default
    n_bins = 18;
    LF_steps = 2:2:14;
    LF_bw = 2;
    HF_steps = 30:5:150;
    tcutEdge = 0.5;
    
    [comdlgrm_1, comdlgrm_z_1, phase2power_1] = cfc_tort_comodulogram(lfp_trials_L1,srate,n_surrogates,n_bins,LF_steps,LF_bw,HF_steps,tcutEdge);
    [comdlgrm_3, comdlgrm_z_3, phase2power_3] = cfc_tort_comodulogram(lfp_trials_L3,srate,n_surrogates,n_bins,LF_steps,LF_bw,HF_steps,tcutEdge);
    
    %% Prep Histogram
    % Mean across freq of interest
    LF_interest = [4 6]; % LF Band of interest
    HF_interest = [70 140]; % HF Band of interest
    % Load 1
    tempMean_LF = mean(phase2power_1(find(LF_steps==LF_interest(1)):find(LF_steps==LF_interest(2)),:,:),1);
    tempMean_HF = mean(tempMean_LF(:,find(HF_steps==HF_interest(1)):find(HF_steps==HF_interest(2)),:),2);
    phase2power_1_HFLF = squeeze(tempMean_HF);
    % Load 3
    tempMean_LF = mean(phase2power_3(find(LF_steps==LF_interest(1)):find(LF_steps==LF_interest(2)),:,:),1);
    tempMean_HF = mean(tempMean_LF(:,find(HF_steps==HF_interest(1)):find(HF_steps==HF_interest(2)),:),2);
    phase2power_3_HFLF = squeeze(tempMean_HF);
    


    %% Plot
    figure(i);
    
    % Load 1
    subplot(2,2,1)
    bar( [phase2power_1_HFLF;phase2power_1_HFLF],'BarWidth',1) 
    xlabel('Theta phase (deg)')
    ylabel('Gamma amplitude (norm.)')
    xlim([0.5 36])
    ylim([0 0.07])
    axis square
    
    subplot(2,2,2)
    contourf(LF_steps, HF_steps,comdlgrm_1',500,'linecolor','none')
    title('Load 1')
    xlabel('Frequency for phase (Hz)')
    ylabel('Frequency for amplitude (Hz)')
    clim([0 7e-4])
    colorbar
    set(gca,'fontsize',10)
    axis square
    
    % Load 3
    subplot(2,2,3)
    bar([phase2power_3_HFLF;phase2power_3_HFLF],'BarWidth',1)
    xlabel('Theta phase (deg)')
    ylabel('Gamma amplitude (norm.)')
    xlim([0.5 36])
    ylim([0 0.07])
    axis square
    
    subplot(2,2,4)
    contourf(LF_steps, HF_steps,comdlgrm_3',500,'linecolor','none')
    title('Load 3')
    xlabel('Frequency for phase (Hz)')
    ylabel('Frequency for amplitude (Hz)')
    clim([0 7e-4])
    colorbar
    set(gca,'fontsize',10)
    axis square
end
