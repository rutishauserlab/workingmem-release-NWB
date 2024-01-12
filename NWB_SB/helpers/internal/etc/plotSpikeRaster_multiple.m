%
% plot raster and PSTH with multiple alignment points 
%
% prepares the data appropriately and then calls plotSpikeRasterMain
%
%urut/MPI/may11
%
function spikesToPlotAll = plotSpikeRaster_multiple(  spikesOfGrp, offsets, limitRange  )
%offsets: offset to previous group, in s
%limitRange: one row for each group, from/to in s of spikes to keep

spikesToPlotAll=[];

%collapse across k, keep j
for k=1:length(spikesOfGrp ) %cell array,each has spikes aligned differently
    spikesToPlot = spikesOfGrp{k};
    
    %only pick from this range
    indsToUse = find( spikesToPlot(:,2)>=limitRange(k,1) & spikesToPlot(:,2)<limitRange(k,2));
    
    if k==1
        spikesToPlotAll= spikesToPlot(indsToUse,:); 
        spikesToPlotAll(:,2) =  spikesToPlotAll(:,2) + offsets(1);
    else
        % append
        spikesToPlotShifted=spikesToPlot;
        spikesToPlotShifted(:,2) = spikesToPlot(:,2) + offsets(k-1);   %shift spike times
        spikesToPlotAll = [ spikesToPlotAll; spikesToPlotShifted(indsToUse,:)]; %#ok<AGROW>
    end
end

%hs=plotSpikeRasterMain( spikesToPlot) ; %, 'colortill', colortill, 'colors', colors, 'range', [1:nrTrialsTot], 'spikeheight', 1, 'spikewidth', 1 );
