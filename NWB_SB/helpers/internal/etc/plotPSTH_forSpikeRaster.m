%
%plot PSTH based on spikesToPlot structure that is used to plot rasters
%Also has options to deal with multiple-align rasters
%
%timeRangesToPlot: if set, only plot these ranges (from/to, one per row)
%binsize to use for PSTH
%triallength in ms
%grpRanges till which trial nr is each grp (inclusive), 
%ie grpRanges=[10 20] means 1...10 and 11..20
%
%plotMode: 1 raw, 2 smoothened
%
%timeRangesToSmooth: use these ranges to smooth data, but plot only the range in timeRangesToPlot; this is to avoid smoothing artifacts at borders
%
%smoothMethod: 1 gaussian kernel (non-causal), 2 exponential (causal)
%
%urut/MPI/may11
function [rangeAll,lineHandles] = plotPSTH_forSpikeRaster( spikesToPlot, timeRangesToPlot, binsize, triallength, grpRanges, plotMode,colors, smoothKernelWidth, timeRangesToSmooth, smoothMethod  )
rangeAll=[0 0];
lineHandles=[];

if nargin<8
    smoothKernelWidth=.200;
end
if nargin<9
    timeRangesToSmooth= timeRangesToPlot;% [timeRangesToPlot(1)+smoothKernelWidth timeRangesToPlot(2) + smoothKernelWidth];
end
if nargin<10
    smoothMethod=1;
end
% timeRangesToPlot = [timeRangesToPlot_in(1)+smoothKernelWidth timeRangesToPlot_in(2) + smoothKernelWidth]; % Testing resize


%grpRanges is the number of subgroups, equivalent to colortill in
%plotSpikeRasterMain.m
kSize=10;
switch(smoothMethod)
    case 1
        smoothKernel = getGaussianKernel(smoothKernelWidth/binsize, kSize);
    case 2
        smoothKernel = getExponentialKernel(smoothKernelWidth/binsize, kSize, 1);
end



grpRanges=[0 grpRanges];
nrGrps = length(grpRanges)-1;

edges=0:binsize:triallength; %+binsize; %1 bin more for trash
    

for k=1:nrGrps
    indsToUse = find( spikesToPlot(:,1)>grpRanges(k) & spikesToPlot(:,1)<=grpRanges(k+1) );
    
    allTimestamps = spikesToPlot(indsToUse,2);
    nrTrials=length(unique(spikesToPlot(indsToUse,1)));

    binnedTotal = histc(allTimestamps, edges);
    binnedTotal = binnedTotal / nrTrials;
    binnedTotal=1/binsize*binnedTotal;
    
    binned = binnedTotal(1:end-1); %remove last element produced by histc because it contains 0 which falsifies our significance calculation.
    %edges = edges(1:end-1);

    if plotMode == 2 & isempty(timeRangesToPlot)
        binned = conv(binned, smoothKernel );
        binned = binned(2:end);
    end
    
    if k>1
        hold on
    end
    
    %plot separatly for each time range part
    if ~isempty( timeRangesToPlot )        
        indEdgesToPlot = [];
        for j=1:size( timeRangesToPlot )
            if j>1
                hold on;
            end
            
            indEdgesToPlot = find( (edges+binsize/2)>=timeRangesToPlot(j,1) & (edges+binsize/2)<timeRangesToPlot(j,2) & edges<triallength);

            if plotMode==2
                
                indEdgesToSmooth = find( (edges+binsize/2)>=timeRangesToSmooth(j,1) & (edges+binsize/2)<timeRangesToSmooth(j,2) & edges<triallength);
                
                
                smoothened = conv(binned(indEdgesToSmooth), smoothKernel );
                if smoothMethod==1
                    toPlot = smoothened(kSize+1:end-kSize);
                else
                    toPlot = smoothened(1:end-kSize);
                end
                
                % discard whatever is in timeRangesToSmooth but not in timeRangesToPlot
                if length(toPlot) ~= length(indEdgesToPlot)
                    indsToRemove = setdiff(indEdgesToSmooth, indEdgesToPlot);
                    
                    indToRemoveRemapped=[];
                    for jj=1:length(indsToRemove)
                        indToRemoveRemapped(jj)=find( indsToRemove(jj) == indEdgesToSmooth);
                    end
                    indKeep = setdiff( 1:length(toPlot), indToRemoveRemapped);
                    
                    toPlot = toPlot(indKeep);
                end
            else
                toPlot = binned(indEdgesToPlot);
            end            
            lineHandles(k) = plot(edges(indEdgesToPlot)+binsize/2, toPlot,'color', [colors{k}], 'linewidth',2);     
            rangeAll = updateYRange(rangeAll, toPlot);
        end
    else
        lineHandles(k) = plot(edges+binsize/2, binned,'color', [colors{k}], 'linewidth',2);
        rangeAll = updateYRange(rangeAll, binned);
    end
end
hold off;


%===
function rangeAll = updateYRange(rangeAll, toPlot)
if rangeAll(1)>min(toPlot)
    rangeAll(1)=min(toPlot);
end
if rangeAll(2)<max(toPlot)
    rangeAll(2)=max(toPlot);
end
