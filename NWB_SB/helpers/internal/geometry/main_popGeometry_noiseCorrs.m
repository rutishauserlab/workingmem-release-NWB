% MAIN_POPGEOMETRY_NOISECORRS MATLAB script for studying population coding.
% This script analyzes the geometry of population coding with emphasis on
% noise correlations and signal/noise axis angles in a two-class problem.
% This script is intended to serve as an example of the analysis used in
% Daume et. al. 

% Functions:
%   noiseCorr_genData - Generates data with specified noise correlations.
%       Inputs: numDims, numDataPoints, corrModes
%       Outputs: generatedData, labels
%   popGeometry_defineAxes - Defines signal and noise axes using classifiers.
%       Inputs: data, labels
%       Outputs: signalAxis, noiseAxis, errorRate, angleMetrics
%   popGeometry_plotSignalNoiseAxes2D - Plots the signal and noise axes in 2D.
%       Inputs: signalAxis, noiseAxis
%       Outputs: 2D plot figure
%   popGeometry_plotSignalNoiseAxes3D - Plots the signal and noise axes in 3D.
%       Inputs: signalAxis, noiseAxis
%       Outputs: 3D plot figure
%   getNoiseVar_alongSignalAxis - Computes the variance of noise along the signal axis.
%       Inputs: data, signalAxis
%       Outputs: noiseVariance
%   LinearClassifier_getSignalAxis - Trains a linear classifier and identifies the signal axis.
%       Inputs: data, labels
%       Outputs: signalAxis, classifierParameters

% Examples:
%   To use this script, set the parameters and call the functions in order.
%   Example of generating data:
%       [data, labels] = noiseCorr_genData(2, 1, 200, 2, 1);
%   Example of defining axes:
%       [signalAxis, noiseAxis, ~, ~] = popGeometry_defineAxes(data, labels);

% Author: Ueli Rutishauser, Jonathan Daume
% Date: August 2023

%% Parameters
Pdims = 2;  % nr of dimensions
Pdims_withTuning = 1;  % nr of dims that have tuning to stimulus
NperGroup = 200;  % nr datapoints per class
addNoiseCommon = 1; % common noise correlations. 0 off, 1 on 
addNoiseSignaldependent = 1; %signal dependent noise correlations
Ngroups = 2; %nr of classes (should not be changed for this script)
modeNoiseAxis = 2;  %1 PCA, 2 covariance matrix based (all data), 3 covariance based (each cluster separately,averaged)

%%
modeClassifier = 2;  % 1 LDA, 2 SVM

%% Simulate the data
[data, labels, dataScrambled] = noiseCorr_genData(Pdims, Pdims_withTuning, NperGroup, addNoiseCommon, Ngroups, addNoiseSignaldependent );

indsGrp1 = find(labels==0);
indsGrp2 = find(labels==1);

CorrMatrix = corr(data);
CorrMatrix_scramble = corr(dataScrambled);

%corr matrix separately for each class
CorrMatrix1 = corr( data(indsGrp1,:) );
CorrMatrix2 = corr( data(indsGrp2,:) );
CorrMatrix_Mean_Noise = (CorrMatrix1+CorrMatrix2)./2;


figure(5);
subplot(2,2,1)
imagesc(CorrMatrix,[0 1]);
colorbar;
title('corr orig data');

subplot(2,2,2)
imagesc(CorrMatrix_scramble,[0 1]);
colorbar;
title('corr scr (signal corr only)');

subplot(2,2,3)
imagesc(CorrMatrix_Mean_Noise,[0 1]);
colorbar;
title('corr per class');


subplot(2,2,4)
imagesc(CorrMatrix-CorrMatrix_scramble,[0 1]);
colorbar;
title('corr full - corr scramble');

%% Train classifier, define signal axis based on classifier trained, define noise axes
%Actual data
[errRate, signalVector_toUse, noiseVector1, noiseVector2, meanAngle, projValues_Cl1, projValues_Cl2, fHyperplane, fBoundLine, PP] ...
    = popGeometry_defineAxes(modeClassifier, data, labels, Pdims, NperGroup, modeNoiseAxis);

%After scramble
[errRate_scr, signalVector_toUse_scr, noiseVector1_scr, noiseVector2_scr, meanAngle_scr, projValues_Cl1_scr, projValues_Cl2_scr, fHyperplane_scr, fBoundLine_scr, PP_scr] ...
    = popGeometry_defineAxes(modeClassifier, dataScrambled, labels, Pdims, NperGroup, modeNoiseAxis);


%% plot data, hyperplane, signalvector
if Pdims == 2
    h1 = figure('position',[300 300 800 400]);
    subplot(121)
    popGeometry_plotSignalNoiseAxes2D(data, NperGroup, fBoundLine, signalVector_toUse, noiseVector1)
    title('Intact');
    axis square
    set(gca,'fontsize',16,'xtick',[],'ytick',[])
    legend off
    box on
    
    
    subplot(122)
    popGeometry_plotSignalNoiseAxes2D(dataScrambled, NperGroup, fBoundLine_scr, signalVector_toUse_scr, noiseVector1_scr)
    title('Removed');
    axis square
    set(gca,'fontsize',16,'xtick',[],'ytick',[])
    box on
    % only plots the first 3 dimensions of the data. So this plot isnt very meaningful for more than 3 dims
elseif Pdims==3
    figure(2);
    popGeometry_plotSignalNoiseAxes3D(data, NperGroup, fHyperplane, signalVector_toUse, noiseVector1, noiseVector2)
    title('data intact');
    
    figure(3);
    popGeometry_plotSignalNoiseAxes3D(dataScrambled, NperGroup, fHyperplane_scr, signalVector_toUse_scr, noiseVector1_scr, noiseVector2_scr)
    title('data scr');
    
end



%% plot signal and noise axes projections and properties
figure(4);
subplot(3,3,1)
histogram(projValues_Cl1);
hold on
histogram(projValues_Cl2);
hold off
mstd = [std(projValues_Cl1) std(projValues_Cl2)];
title( ['data std=' num2str(mstd)] );

subplot(3,3,2)
histogram(projValues_Cl1_scr);
hold on
histogram(projValues_Cl2_scr);
hold off
mstd_scr = [std(projValues_Cl1_scr) std(projValues_Cl2_scr)];
title( ['scramble std=' num2str(mstd_scr)] );

subplot(3,3,3)
bar([1 2], [mean(mstd) mean(mstd_scr)] );
set(gca,'XTickLabel', {'data','scr'});
ylabel('std on signal axis');

subplot(3,3,4)
bar([1 2], rad2deg([meanAngle meanAngle_scr]) );
set(gca,'XTickLabel', {'data','scr'});
ylabel('signal-noise angle');
ylim([0 90]);

subplot(3,3,5)
bar([1 2], [errRate errRate_scr]);
set(gca,'XTickLabel', {'data','scr'});
ylabel('error rate');
title(['classifier=' num2str(modeClassifier)]);

subplot(3,3,6)
bar([1 2], [PP PP_scr]);
set(gca,'XTickLabel', {'data','scr'});
ylabel('projected precision (PP)');

