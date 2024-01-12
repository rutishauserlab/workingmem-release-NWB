function [comdlgrm, comdlgrm_z, phase2power] = cfc_tort_comodulogram(datMat,srate,n_surrogates,n_bins,LF_steps,LF_bw,HF_steps,tcutEdge)

%% [comdlgrm, comdlgrm_z] = cfc_tort_comodulogram(LFP,fs,n_surrogates,n_bins,LF_steps,LF_bw,HF_steps,tcutEdge)
%
% Compute PAC comodulogram with Tort's method in a single channel and condition
%
% Input
% datMat: 2d matrix containing time samples of interest per trial (samples x trials)
% srate: sampling rate
% n_surrogates: number of surrogates for obtaining z-scored comodulogram (0 = no surrogates computed; default: 0)
% n_bins: number of bins to compute modulation index (default: 18 bins)
% LF_steps: center frequencies for phase signal in Hz (default: 2:2:14)
% LF_bw: Bandwidth for phase signal in Hz (default: 2)
% HF_steps: center frequencies for amplitude signal in Hz (default: 30:5:150; bandwidth is determined by phase signal frequency)
% tcutEdge: time to cut off at the edges of each trial to prevent filter artifacts in s (full time will be cutoff at the beginning and end of trial; default: 0 (no cutoff)) 
%
% Output
% comdlgrm: raw MI comodulugram; LF_steps x HF_steps
% comdlgrm_z: z-scored MI comodulugram; if n_surrogates > 0; LF_steps x HF_steps
% JD 2023

rng('shuffle')

%% defaults
if nargin < 3 || isempty(n_surrogates)
    n_surrogates = 0; % number of surrogates; if 0, no surrogates are computed
end

if nargin < 4 || isempty(n_bins)
    n_bins = 18;
end

if nargin < 5 || isempty(LF_steps)
    LF_steps = 2:2:14; % center frequencies for phase signal in Hz
end

if nargin < 6 || isempty(LF_bw)
    LF_bw = 2;% 2; bandwidth in Hz
end

if nargin < 7 || isempty(HF_steps)
    HF_steps = 30:5:150; % center frequencies for amplitude signal in Hz (bandwidth is dependent on corresponding phase signal)
end

if nargin < 8 || isempty(tcutEdge)
    tcutEdge = 0; % time removed at the beginning and end of each trial after filtering and Hilbert to prevent edge artifacts; total time removed = tcutEdge*2
end

cutTrial = logical(tcutEdge);


n_trials       = size(datMat,2);
n_samples_long = size(datMat,1);

%% initalize output matrix
clear comdlgrm
n_HF = length(HF_steps);
n_LF = length(LF_steps);
comdlgrm = nan(n_LF,n_HF);
comdlgrm_surr_mean = nan(n_LF,n_HF);
comdlgrm_surr_std = nan(n_LF,n_HF);
phase2power = nan(n_LF,n_HF,n_bins);
%% Make filter input

EEG = [];
EEG.srate  = srate;
EEG.pnts   = n_samples_long;
EEG.trials = n_trials;
EEG.nbchan = 1;
EEG.data   = datMat;

%% Loop through all frequency pairs
for i_phase = 1:n_LF
    disp('Filtering and computing hilbert transform for phase signal...');

    % filter phase signal
    LF_bp = [LF_steps(i_phase)-LF_bw/2 LF_steps(i_phase)+LF_bw/2];
    EEG_p = pop_eegfiltnew(EEG,LF_bp(1),LF_bp(2));
    
    
    %% Hilbert
    phase_long = angle(hilbert(squeeze(EEG_p.data)));
    clear EEG_p
    
    %% Cutting out time window of interest
    if cutTrial
        disp('Cutting out time window of interest in each trial for further analysis...')
        t2cut_samples = tcutEdge*srate;
        phase_toi     = phase_long(t2cut_samples+1:end-t2cut_samples,:);
    else
        phase_toi = phase_long;
    end
    
    % Transfer data to one long vector
    n_samples = size(phase_toi,1);
    numpoints = n_samples * n_trials;
    phase     = reshape(phase_toi,numpoints,1);
    
    clear phase_toi phase_long
    
    %%
    for i_amplitude = 1:n_HF
        fprintf('Computing PAC between %dHz amplitude and %dHz phase frequency...\n',HF_steps(i_amplitude),LF_steps(i_phase));

        HF_bp = [HF_steps(i_amplitude)-LF_steps(i_phase) HF_steps(i_amplitude)+LF_steps(i_phase)];
        EEG_A = pop_eegfiltnew(EEG,HF_bp(1),HF_bp(2));
        
        %% Hilbert
        amplitude_long = abs(hilbert(squeeze(EEG_A.data)));
        clear EEG_A
        
        %% Cutting out time window of interest
        if cutTrial
            t2cut_samples = tcutEdge*srate;
            amplitude_toi = amplitude_long(t2cut_samples+1:end-t2cut_samples,:);
        else
            amplitude_toi = amplitude_long;
        end
        % Transfer data to one long vector
        amplitude = reshape(amplitude_toi,numpoints,1);
        
        clear amplitude_toi amplitude_long
        
        %% Code for calculating MI like Tort
        
        phase_degrees = rad2deg(phase); % Phases in degrees
        % Bining the phases
        step_length = 360/n_bins;
        phase_bins = -180:step_length:180;
        [~,phase_bins_ind] = histc(phase_degrees,phase_bins);
        clear phase_degrees
        
        % Averaging amplitude time series over phase bins
        amplitude_bins = nan(n_bins,1);
        
        for bin = 1:n_bins
            amplitude_bins(bin,1) = mean(amplitude(phase_bins_ind==bin));
        end
        
        % Normalize amplitudes
        P = amplitude_bins./repmat(sum(amplitude_bins),n_bins,1);
        
        % Compute modulation index and store in comodulogram
        mi = 1+sum(P.*log(P))./log(n_bins);
        comdlgrm(i_phase,i_amplitude) = mi;
        phase2power(i_phase,i_amplitude,:) = P;
        
        clear amplitude_bins mi P
        
        %% Compute surrogates
        if n_surrogates
            
             mi_surr = nan(n_surrogates,1);

            % reshape back to trials for shuffling
            amplitude_trials = reshape(amplitude,n_samples,n_trials);
            clear amplitude
            
            %% compute surrogate values
            disp('Computing surrogate data...');
            for s=1:n_surrogates
                randind = randperm(n_trials);
                surrogate_amplitude = reshape(amplitude_trials(:,randind),numpoints,1);
                amplitude_bins_surr = zeros(n_bins,1);
                for bin = 1:n_bins
                    amplitude_bins_surr(bin,1) = mean(surrogate_amplitude(phase_bins_ind==bin));
                end
                P_surr = amplitude_bins_surr./repmat(sum(amplitude_bins_surr),n_bins,1);
                mi_surr(s) = 1+sum(P_surr.*log(P_surr))./log(n_bins);
            end
            
            %% fit gaussian to surrogate data, uses normfit.m from MATLAB Statistics toolbox
            [surrogate_mean,surrogate_std]=normfit(mi_surr);
            
            comdlgrm_surr_mean(i_phase,i_amplitude) = surrogate_mean;
            comdlgrm_surr_std(i_phase,i_amplitude) = surrogate_std;
            
            clear mi_surr
        end
        
    end %loop amplitude
end %loop phase

% z-transform raw comod
comdlgrm_z = (comdlgrm - comdlgrm_surr_mean) ./ comdlgrm_surr_std;
disp('Done')


