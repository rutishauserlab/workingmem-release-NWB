%
% plot raster and PSTH, grouped by multiple conditions indicated by
% different colors. Pulls together various code needed to do this with one
% simple call.
%
% 
% urut/Oct'22
function [axRaster,axPSTH]=plotSpikeRaster_multiple_simplified( plotMode, offsets, limitRange_forRaster, limitRange_forPSTH, markerPos, timestampsOfCell, periods_all, periods_indx_toPlot, grpStrs, subplot_Raster, subplot_PSTH, spikePlotParams )

% plotMode   1 Raster only, 2 PSTH only, 3 both
% periods_all  list of all available periods. Cell Array, one entry per alignment point
% periods_indx_toPlot   cell array of indices of to be plotted periods. Each entry is one grouping. Same indices are used for all align points
% subplot_Raster = [ subplotX subplotY subplotstouse ]
% subplot_PSTH = [ subplotX subplotY subplotstouse ]

% Parameters
%spikewidth  = 1.5;
%spikeheight = 1;
%smoothKernelWidth = 200;
%colors1 = {[1.0 0 0],[1.0 0.8 0.6], [0.2 0.8 1.0],[0 0 1]};   % splitup into 4
%binsizePSTH = 250; % bin size for PSTH
axRaster=[];
axPSTH=[];

spikesToPlot_all=[];
for alignPointsId=1:length(periods_all)
    periods_ofAlignPoint = periods_all{alignPointsId};
    
    % prepare timestamps of groups
    periods_ofGrp=[];
    for grpId=1:length(periods_indx_toPlot)
        periods_ofGrp{grpId} = periods_ofAlignPoint( periods_indx_toPlot{grpId}, :);
    end
    
    [spikesToPlot_align1, colortill1, nrTrialsTot1] = getTimestamps_forSortedRasterPlot( timestampsOfCell, periods_ofGrp);
    
    spikesToPlot_all{alignPointsId} = spikesToPlot_align1;
end

spikesToPlot_multiple = plotSpikeRaster_multiple(  spikesToPlot_all, offsets, limitRange_forRaster  );
        
if plotMode==1 || plotMode==3
    %Raster

    axRaster = subplot( subplot_Raster(1), subplot_Raster(2), subplot_Raster(3:end) );
    
    plotSpikeRasterMain( spikesToPlot_multiple, 'colortill', colortill1, 'colors', spikePlotParams.colors1, 'range', [1:nrTrialsTot1], 'spikeheight', spikePlotParams.spikeheight, 'spikewidth', spikePlotParams.spikewidth );
    h=plotTrialMarkers( markerPos, [0.7 0.7 0.7],[], 1 );
    ylabel('trial nr (re-sorted)');
end

if plotMode==2 || plotMode==3
    %PSTH
    axPSTH = subplot( subplot_PSTH(1), subplot_PSTH(2), subplot_PSTH(3:end) );

    triallengthCombined=limitRange_forPSTH(end);
%  spikesToPlot, timeRangesToPlot, binsize, triallength, grpRanges, plotMode,colors, smoothKernelWidth, timeRangesToSmooth, smoothMethod 
    plotModePSTH=2; %1 raw, 2 smooth
    [~,hs_raster]  = plotPSTH_forSpikeRaster( spikesToPlot_multiple, limitRange_forPSTH, spikePlotParams.binsizePSTH, triallengthCombined, colortill1, plotModePSTH, spikePlotParams.colors1, spikePlotParams.smoothKernelWidth );
    h=plotTrialMarkers( markerPos, [0.7 0.7 0.7],[],1 );
    legend(hs_raster, grpStrs, 'Location', 'SouthEast' );
end

