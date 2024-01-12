function [figOut] = plot_PAC_neuron(spike_times,datMat,srate,low_fb,HF_steps,tcutEdge)

%% plot_PAC_neuron(spike_times,datMat,srate,low_fb,HF_steps,tcutEdge)
% Plots the data table used to compute the GLMs for PAC neuron selection
%
% Input:
% spike_times:  matrix containing the spike times from a given neuron (samples x trials). The
%               matrix should contain zeros and ones. Ones indicate a spike
%               at a certain sample in a trial.
% datMat:       LFP data for a given channel (samples x trials)
% srate:        sampling rate
% low_fb:       frequency band for the low frequency signal in Hz [lower edge, higher edge];
% HF_steps:     center frequencies for the higher amplitude in Hz [lower edge : stepsize : higher end]
% tcutEdge:     time to cut off at the edges of each trial to prevent filter artifacts in s (full time will be cutoff at the beginning and end of trial;  0 = no cutoff) 
%
% JD 2023

n_phase_bins = 10; 
rng('shuffle')
spike_times = logical(spike_times);

%% init
n_trials       = size(datMat,2);
n_samples_long = size(datMat,1);

EEG = [];
EEG.srate  = srate;
EEG.pnts   = n_samples_long;
EEG.trials = n_trials;
EEG.nbchan = 1;
EEG.data   = datMat;

EEG_LF = pop_eegfiltnew(EEG,low_fb(1),low_fb(2));
clear EEG

%% Hilbert
disp('Computing Hilbert transform...')
lf_phase = angle(hilbert(squeeze(EEG_LF.data)));

clear EEG_LF
%% Wavelet convolution parameters

f            = HF_steps;
wl_time      = -1:1/srate:1;
wl_time_half = (length(wl_time)-1)/2;
n_freq       = length(f);

% FFT parameters (use next-power-of-2)
n_samples_wl          = length(wl_time);
n_samples_data        = n_samples_long;
n_samples_convolution = n_samples_wl+n_samples_data-1;
n_samples_conv_pow2   = pow2(nextpow2(n_samples_convolution));
wavelet_cycles        = 7;

% Computation and FFT of wavelets
fprintf('Computation and fft of wavelets for all %d specified frequencies...\n',n_freq)
% compute wavelets for each frequency
wavelets = zeros(n_samples_wl,n_freq);
for fi = 1:n_freq
    wavelets(:,fi) = exp(2*1i*pi*f(fi).*wl_time) .* exp(-wl_time.^2./(2*(wavelet_cycles/(2*pi*f(fi)))^2));
end

% Fourier transform of the wavelets
wavelets_fft = fft(wavelets, n_samples_conv_pow2);

%% wavelet transform
hf_amp = nan(n_samples_data, n_trials, n_freq);
disp('Computing wavelet transform..')

% FFT of data (note: this doesn't change on frequency iteration)
fft_datMat = fft(datMat,n_samples_conv_pow2);

% compute convolution for each frequency
for i_freq=1:n_freq
    
    % duplicate wavelets to match number of channel
    wl = repmat(wavelets_fft(:,i_freq), [1 n_trials]);
    
    % run convolution
    convResult = ifft(wl.*fft_datMat,n_samples_conv_pow2);
    convResult = convResult(1:n_samples_convolution,:); % here the extra points from the power-of-2 FFT are removed
    convResult = convResult(wl_time_half+1:end-wl_time_half,:);
    
    % Put averaged data to tf-matrix
    hf_amp(:,:,i_freq) = abs(convResult);
    
end %freq

%% Cut data to final toi
t2cut_samples = tcutEdge*srate;
lf_phase = lf_phase(t2cut_samples+1:end-t2cut_samples,:);
hf_amp = hf_amp(t2cut_samples+1:end-t2cut_samples,:,:);
spike_times = spike_times(t2cut_samples+1:end-t2cut_samples,:);

n_samples = size(lf_phase,1);

%% Make long vectors
lf_phase = reshape(lf_phase, [n_samples*n_trials 1]);
hf_amp = reshape(hf_amp, [n_samples*n_trials n_freq]);
spike_times = reshape(spike_times, [n_samples*n_trials 1]);

%% normalize gamma amp and average across all frequencies
hf_amp = zscore(hf_amp); % normalize across all trials
hf_amp = mean(hf_amp,2); % average across freqs

%% determine spike count in phase and amplitude categories
% gamma
med_amp = median(hf_amp);
gamma_cat = zeros(size(hf_amp));
gamma_cat(hf_amp < med_amp) = 1;
gamma_cat(hf_amp >= med_amp) = 2;
n_gamma_bins = max(gamma_cat);

%theta
step_length = 2*pi/n_phase_bins;
phase_bins = -pi:step_length:pi;
[~,phase_bins_ind] = histc(lf_phase,phase_bins);

Spike_count = zeros(n_phase_bins*n_gamma_bins,1);
Gamma_amp = zeros(n_phase_bins*n_gamma_bins,1);
Lf_phase = zeros(n_phase_bins*n_gamma_bins,1);
counter = 0;
for i_phasebin = 1:n_phase_bins
    for i_ampbin = 1:n_gamma_bins
        counter = counter + 1;
        Spike_count(counter) = sum(spike_times(gamma_cat == i_ampbin & phase_bins_ind == i_phasebin));
        Gamma_amp(counter) = i_ampbin;
        Lf_phase(counter) = phase_bins(i_phasebin);
    end
end
Gamma_amp = categorical(Gamma_amp);
Cos_theta = cos(Lf_phase);
Sin_theta = sin(Lf_phase);

T = table(Spike_count,Lf_phase,Cos_theta,Sin_theta,Gamma_amp);


%% GLM
disp('Computing GLM comparisons...')
model1_full = fitglme(T,'Spike_count ~ (Cos_theta + Sin_theta) *  Gamma_amp', 'Distribution','Poisson','FitMethod','Laplace');
model2_noInt = fitglme(T,'Spike_count ~ Cos_theta + Sin_theta +  Gamma_amp', 'Distribution','Poisson','FitMethod','Laplace');
compare_stats_int = compare(model2_noInt,model1_full);

model3_noGamma = fitglme(T,'Spike_count ~ (Cos_theta + Sin_theta) *  Gamma_amp - Gamma_amp', 'Distribution','Poisson','FitMethod','Laplace');
compare_stats_meg = compare(model3_noGamma,model1_full);

%% plot
disp('Plotting..')
figOut = figure('position',[100 100 600 600]);

edges = rad2deg(T.Lf_phase(T.Gamma_amp == '1'));
edges(end+1) = -edges(1);
sc_highgamma = T.Spike_count(T.Gamma_amp == '2');
sc_lowgamma = T.Spike_count(T.Gamma_amp == '1');
maxsc = max([sc_highgamma;sc_lowgamma])+5;
subplot(221)
histogram('BinEdges', edges', 'BinCounts', sc_highgamma)
set(gca,'ylim', [0 maxsc],'xtick',[-180 0 180],'fontsize',18)
title('High gamma power')
ylabel('Count')
xlabel('Theta Phase')
subplot(223)
histogram('BinEdges', edges', 'BinCounts', sc_lowgamma)
set(gca,'ylim', [0 maxsc],'xtick',[-180 0 180],'fontsize',18)
title('Low gamma power')
xlabel('Theta Phase')

txt_int = sprintf('LRstat = %.2f, p = %.4f',compare_stats_int.LRStat(2),compare_stats_int.pValue(2));
txt_meg = sprintf('LRstat = %.2f, p = %.4f',compare_stats_meg.LRStat(2),compare_stats_meg.pValue(2));
x = 300;
y = max(sc_lowgamma)/2;
text(x,y+30,'Model comparisons:','fontsize',18,'FontWeight','bold')
text(x,y+23,'Model 1 vs Model 2 (Interaction):','fontsize',18)
text(x, y+18, txt_int,'fontsize',18)
text(x,y+10,'Model 1 vs Model 3 (Gamma):','fontsize',18)
text(x, y+5, txt_meg,'fontsize',18)



