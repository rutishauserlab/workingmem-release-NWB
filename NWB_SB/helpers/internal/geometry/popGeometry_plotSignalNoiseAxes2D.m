function popGeometry_plotSignalNoiseAxes2D(data, NperGroup, fBoundLine, signalVector_toUse, noiseVector1)
scatter( data(1:NperGroup,1), data(1:NperGroup,2),[],[244 225 129]/255,'filled'); 
hold on
scatter( data(NperGroup+1:end,1), data(NperGroup+1:end,2), [],[223 159 207]/255,'filled');
% plot decision boundary
minVal = min(data(:,1))+7;
maxVal = max(data(:,1))-8;

line([minVal maxVal],[fBoundLine(minVal) fBoundLine(maxVal)],'color',[0.3 0.3 0.3],'linewidth', 2)



%==
meanX = mean(data(:,1));
meanY = mean(data(:,2));

%plot signalVector
scaleFactor=10;
hSignal=line( [0 signalVector_toUse(1)]*scaleFactor + meanX, [ 0 signalVector_toUse(2)]*scaleFactor +  meanY, 'color', [121 224 223]/255, 'linewidth', 3);

%plot noiseVectors
scaleFactor=10;
hNoise1=line( [0 noiseVector1(1)]*scaleFactor + meanX, [ 0 noiseVector1(2)]*scaleFactor + meanY, 'color', [152 103 255]/255, 'linewidth', 3); 

legend([hSignal hNoise1], {'Signal axis', 'Noise axis'});

% %==same range so things look orthogonal
mVal = max(abs(data(:)));
rangeCommon=[-mVal mVal];
xlim(rangeCommon);
ylim(rangeCommon);

xlabel('Neuron 1');
ylabel('Neuron 2');