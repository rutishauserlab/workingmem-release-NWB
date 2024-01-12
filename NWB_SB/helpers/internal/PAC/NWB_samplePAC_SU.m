function [figOut] = NWB_samplePAC_SU(nwbAll, all_units, paths)
% - nwbAll: a cell array of all nwb files imported in NWB_SBCAT_import_main
% operation: 




% whichSession = ‘P101T’;
% whichNeuron = ‘RPH7_2’;
% whichChannel = ‘RPH5’;
% whichLoad = 3;


fs = filesep;

sample_sub = 26;
sample_ses = 1;
% lfp_ID = 34;
% su_id = 16;


% Access data loaded from import script
% Check if session of interest was imported
sub_ses_arr = NaN(length(nwbAll),2);
for k = 1:length(nwbAll)
    sub_ses_arr(k,1) = str2double(nwbAll{k}.general_subject.subject_id);
    sub_ses_arr(k,2) = str2double(nwbAll{k}.general_session_id);
end

find_example = find(sample_sub == sub_ses_arr(:,1) & sample_ses == sub_ses_arr(:,2),1);
if isempty(find_example)
    warning('Example not found. Please import session (5) sub-4_ses-1_P60CS.')
    figOut = [];
    return
else
    nwb = nwbAll{find_example};
end

for i = 33
    for j = 20
        lfp_ID = i;
        su_ID = j;     
        % Load Maint ts
        tsMaint = nwb.intervals_trials.vectordata.get('timestamps_Maintenance').data.load();
        % Electrode labels for troubleshooting
        electrode_labels = string(nwb.general_extracellular_ephys_electrodes.vectordata.get('location').data.load());
        electrode_origChan = string(nwb.general_extracellular_ephys_electrodes.vectordata.get('origChannel').data.load());
        lfp_electrodes = int32(nwb.acquisition.get('LFPs').electrodes.data.load())+1; % Offset to 1-indexing; 
        lfp_labels = electrode_labels(lfp_electrodes);
        lfp_origChan = electrode_origChan(lfp_electrodes);
        su_electrodes = nwb.units.electrodes.data.load()+1; % Offset to 1-indexing
        su_labels = electrode_labels(su_electrodes);
        su_origChan = electrode_origChan(su_electrodes);
        %% Load & Format LFPs
        fprintf('Loading LFP: Orig chan %s (%s)\n',lfp_origChan(i),lfp_labels(i))
        fprintf('Loading SU: Orig chan %s (%s)\n',su_origChan(j),su_labels(j))
        lfp_raw = nwb.acquisition.get('LFPs').data.load();
        lfp_chan = lfp_raw(lfp_ID,:);
        lfp_start_ts = nwb.acquisition.get('LFPs').starting_time;
        lfp_srate = nwb.acquisition.get('LFPs').starting_time_rate; % Should be 400
        lfp_ts = lfp_start_ts + (0:length(lfp_chan))*(1/lfp_srate); % tspan with 1/hz increments, used in ts filter
        
        trial_pre = 1.5;
        trial_post = 4;
        
        lfp_trials = []; % Acquire LFP segments per trial
        for k=1:length(tsMaint)
            in_period = lfp_ts >= (tsMaint(k)-trial_pre) & lfp_ts <= (tsMaint(k)+trial_post);
            period_lfp = lfp_chan(in_period)';
            lfp_trials = [lfp_trials,period_lfp]; %#ok<AGROW> % It's just 140 trials. 
        end
        
        %% Load & Format SUs
        spike_ts_all = nwb.units.spike_times().data.load();
        spike_ts_index = nwb.units.spike_times_index.data.load();
        
        if su_ID == 1
            su_ts = spike_ts_all(1:spike_ts_index(su_ID));
        else
            su_ts = spike_ts_all(spike_ts_index(su_ID-1)+1:spike_ts_index(su_ID));
        end
        
        
        binSize = 1/lfp_srate;
        trials_pointProcess = []; 
        for k=1:length(tsMaint)
            % Extract trial spikes
            in_trial = su_ts > (tsMaint(k)-trial_pre) & su_ts <= (tsMaint(k)+trial_post);
            su_ts_trial = su_ts(in_trial) - (tsMaint(k)-trial_pre); % Spike times in trial relative to period onset. All values should be positive
            % trial_tspan = tsMaint(i) + (-trial_pre:(1/lfp_srate):trial_post);
        
            % Converting to point process
            binaryProcess = zeros(1,(trial_post+trial_pre)*lfp_srate);
            for j=1:length(su_ts_trial)
                binIndex = ceil(su_ts_trial(j)/binSize);
                binaryProcess(binIndex) = 1;
            end
            % Append point process to [point process X trials] matrix
            trials_pointProcess = [trials_pointProcess,binaryProcess'];
        end
        
        %% Filtering to correct Load 3 Trials & Removing Artifacts. 
        is_load_3 = nwb.intervals_trials.vectordata.get('loads').data.load() == 3;
        is_correct = nwb.intervals_trials.vectordata.get('response_accuracy').data.load();
        artifact_path = [paths.code fs 'helpers' fs 'internal' fs 'PAC' fs 'artif_trials.mat'];
        loadIn = load(artifact_path); artifact_trials = loadIn.artif_trials;
        is_not_artifact = ~ismember(1:length(tsMaint),artifact_trials)';
        
        keep_trials = is_load_3 & is_correct & is_not_artifact;
        lfp_trials_L3filtered = lfp_trials(:,keep_trials);
        trials_pointProcess_L3filtered = trials_pointProcess(:,keep_trials);
        
        %% Call plotting function
        srate = lfp_srate;
        low_fb = [3 7];
        HF_steps = 70:5:140;
        tcutEdge = 1.5;
        
        [figOut] = plot_PAC_neuron(trials_pointProcess_L3filtered,lfp_trials_L3filtered,srate,low_fb,HF_steps,tcutEdge);
    end
end
fprintf('debugging')
end

