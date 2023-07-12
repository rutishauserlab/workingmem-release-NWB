%
%prepare groups of trials for plotting of raster (with color info)
%
%indsOfGrps: cell array of list of periods (each a list of trials)
%
%urut/april13; initially derived from getTimestampsBubbles.m
function [spikesToPlot, colortill, nrTrialsTot] = getTimestamps_forSortedRasterPlot( timestampsOfCell, indsOfGrps)

spikesToPlot=[];
trialNr=0;
colortill=[];


for k=1:length(indsOfGrps)
    
    spikesOfCat = getRelativeTimestamps(timestampsOfCell, indsOfGrps{k} );
       
    if k==1
        colortill(k) = length(spikesOfCat);
    else
        if k==2
            offset=1;
        else
            offset=0;
        end
        colortill(k) = colortill(k-1)+length(spikesOfCat)+offset;
    end
    
    for kk=1:length(spikesOfCat) %for each trial
        trialNr=trialNr+1;
        spikesToPlot = [ spikesToPlot; [repmat(trialNr,length(spikesOfCat{kk}), 1) spikesOfCat{kk}] ];    
    end
end

nrTrialsTot=trialNr;
