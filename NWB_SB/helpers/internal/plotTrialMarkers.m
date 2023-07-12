%
%plot stop/start markers on a raster.
%
%plotModeBottom: if enabled (1), plot the lines at the back, not over what is in the plot
%
%
function h = plotTrialMarkers( markerPos, col,lw, plotModeBottom )
if nargin<2
    col='r';
end
if nargin<3
    lw=2;
end
if isempty(lw)
    lw=2;
end
if nargin<4
    plotModeBottom=0;
end

h=[];
for i=1:length(markerPos)
    maxY=ylim;
    h(i) = line([markerPos(i) markerPos(i)],[maxY(1) maxY(2)]);
    set(h,'color',col);
    set(h','linewidth',lw);
    set(h,'LineStyle','--')
end

if plotModeBottom
    uistack(h,'bottom');
end
         