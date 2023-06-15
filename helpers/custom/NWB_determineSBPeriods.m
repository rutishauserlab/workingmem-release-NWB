%determines the timestamp for all trials specified with inds from the event file
%
%inds is onset of stimulus
%inds+1 is offset of stimulus
%
%
%before/after Offset is in ms. it determines the length of the trial, relative to the timepoint specified by ind
%
%urut/jan07
% SSullivan july/17 slightly modified for SB analysis
function periods = NWB_determineSBPeriods( tsEvents, beforeOffset, afterOffset )

periods = zeros(length(tsEvents),3);
for i=1:length(tsEvents)
   periods(i,1:3) = [ i tsEvents(i)-beforeOffset tsEvents(i)+afterOffset];
end

