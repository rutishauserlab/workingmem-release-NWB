function [projValues] = getNoiseVar_alongSignalAxis(data, signalVector)

projValues=[];
for k=1:size(data,1)
    projValues(k) = dot( data(k,:), signalVector)/norm(signalVector);
end