%
%reference timestamps to beginning of the trial
%
%periods: each row is one trial. 3 columns: trial nr, from, to.
%returns a cell array; each item contains the timestamps of one trial. the
%number of trials is equal to the number of rows in periods.
%
%urut/dec07
function trialsTimestamps = getRelativeTimestamps(timestampsOfCell, periods)

trialsTimestamps = getTimestampsOfTrials( timestampsOfCell, periods(:,2:3) );

for i=1:length(trialsTimestamps)
    trialsTimestamps{i} = (trialsTimestamps{i} - periods(i,2)); %remove offset and convert to ms
end