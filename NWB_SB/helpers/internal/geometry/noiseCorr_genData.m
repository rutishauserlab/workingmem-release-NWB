% Simulate data to examine effect of noise correlations in populatin of neurons
%
% addNoiseCommon: add noise correlations
% NperGroup: number datapoints per group
% Ngroups: number of groups/classes. Default is =2, can also accomodate =3.
% addNoiseSignaldependent: 0 no, 1 add signal dependent noise
%
function [data, labels, dataScrambled] = noiseCorr_genData(Pdims, Pdims_withTuning, NperGroup, addNoiseCommon, Ngroups, addNoiseSignaldependent)

if nargin<5
    Ngroups=2;
end
if nargin<6
    addNoiseSignaldependent=0;
end

%% parameters
%variance. Randomize so not all the same
sigmaDefault = rand(Ngroups,Pdims)*2+1;
sigmaData_ofGrp = ones(Ngroups,Pdims) .* sigmaDefault;

%mean.  Each entry is a different group
muTuned   = [10  5 20 20]; %mean for each class for tuned dimensions
muUntuned = [0  0 0 0]; %mean for each class for untuned dimensions

muData_ofGrp=[];
for k=1:Ngroups
    muData_ofGrp(k,1:Pdims_withTuning) = ones(1,Pdims_withTuning) * muTuned(k);
    if Pdims-Pdims_withTuning>0
        muData_ofGrp(k,Pdims_withTuning+1:Pdims) = ones(1,Pdims-Pdims_withTuning) * muUntuned(k);
    end
end

sigmaNoiseCommon = 2;  % so that std is variable
muNoiseCommon = 0;

sigmaNoiseSignal=[2 0.2 0.2 0.2];  % each class
muNoiseSignal = 0;

%% simulate the data
labels=[];
data=[];
for k=1:Ngroups
    labelsOfGroup = ones(NperGroup,1)*(k-1);

    dataOfGroup = randn(NperGroup, Pdims) .* sigmaData_ofGrp(k,:) + muData_ofGrp(k,:);
    
    %add signal correlations
    if addNoiseSignaldependent
        %add signal dependent corrs to the tuned dimensions
            for jj=1:size(dataOfGroup,1)    

                noiseCorrSignal = randn*sigmaNoiseSignal(k) + muNoiseSignal;    
                dataOfGroup(jj,:) = dataOfGroup(jj,:) + noiseCorrSignal;
            end

    end

    labels = [labels; labelsOfGroup];
    data = [ data; dataOfGroup ];
end


%%
if addNoiseCommon
    for k=1:size(data,1)    
        noiseCorrCommon = randn*sigmaNoiseCommon + muNoiseCommon;    
        data(k,:) = data(k,:) + noiseCorrCommon;
    end
end


%% scramble, within group
dataScrambled = data;
for j=1:Pdims
    for k=1:Ngroups
        indsOfGroup = find( labels==(k-1) );
        dataScrambled(indsOfGroup,j) = dataScrambled(indsOfGroup( randperm(length(indsOfGroup))), j);
    end
end


