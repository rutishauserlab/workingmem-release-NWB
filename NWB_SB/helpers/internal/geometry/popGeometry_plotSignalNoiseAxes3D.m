function popGeometry_plotSignalNoiseAxes3D(data, NperGroup, fHyperplane, signalVector_toUse, noiseVector1, noiseVector2)

plot3( data(1:NperGroup,1), data(1:NperGroup,2), data(1:NperGroup,3), 'ro');
hold on
plot3( data(NperGroup+1:end,1), data(NperGroup+1:end,2),data(NperGroup+1:end,3), 'bo');
hold off

%== plot the hyperplane. Define 4 points
maxVals = max( abs(data) );
corners = [0 maxVals(2) maxVals(3);
           0 -maxVals(2) maxVals(3);
           0 -maxVals(2) -maxVals(3);
           0 maxVals(2) -maxVals(3)];

planePoints=[];
for k=1:4
    hypPoint = zeros(1, size(data,2));
    hypPoint(1:3) = corners(k,:);
    hypPoint(1) = fHyperplane(hypPoint);

    planePoints = [ planePoints; hypPoint];
end

patch( planePoints(:,1)', planePoints(:,2)', planePoints(:,3)', 'red', 'FaceAlpha', 0.2 );
%==

%plot signalVector
scaleFactor=10;
hSignal=line( [0 signalVector_toUse(1)]*scaleFactor, [ 0 signalVector_toUse(2)]*scaleFactor,[ 0 signalVector_toUse(3)]*scaleFactor, 'color', 'm', 'linewidth', 3);
scaleFactor=-10;
line( [0 signalVector_toUse(1)]*scaleFactor, [ 0 signalVector_toUse(2)]*scaleFactor,[ 0 signalVector_toUse(3)]*scaleFactor, 'color', 'm', 'linewidth', 3);

%plot noiseVectors
scaleFactor=10;
hNoise1=line( [0 noiseVector1(1)]*scaleFactor, [ 0 noiseVector1(2)]*scaleFactor,[ 0 noiseVector1(3)]*scaleFactor, 'color', 'g', 'linewidth', 3);
hNoise2=line( [0 noiseVector2(1)]*scaleFactor, [ 0 noiseVector2(2)]*scaleFactor,[ 0 noiseVector2(3)]*scaleFactor, 'color', 'g', 'linewidth', 3);

legend([hSignal hNoise1 hNoise2], {'Signal axis', 'Noise axis 1', 'Noise axis 2'});

%==same range so things look orthogonal
mVal = max(abs(data(:)));
rangeCommon=[-mVal mVal];
xlim(rangeCommon);
ylim(rangeCommon);
zlim(rangeCommon);

xlabel('x');
ylabel('y');
zlabel('z');
