%
% calculates the coefficient of variation metric CV2, as described in Holt et al 1996, J Neurophysiol
%
% ISIs are in sec and expected to be in order
%
%
%returns:
%CV2: mean CV2 over all ranges
%CV2_all: CV2 for each pair of consecutive ISIs
%CV2_all_diff: the difference between each consecutive ISIs, useful to calculate the CV2 only for sub-ranges of intervals
%
%
%ignoreMode: 0 take all, 1 ignore very large ISIs; use this if only sub-periods of the spike-train are used, so these have "wholes" that create large artificial ISIs
%
%
%urut/jan13
function [CV2, CV2_all, CV2_all_sum, CV2_all_diff] = calcCV2(ISIs, ignoreMode)
if nargin<2
    ignoreMode=0;
end

upperLim=1; % max ISI duration in sec

CV2_all=[];
validISI=[];
for j=2:length(ISIs)
    
    CV2_all_sum(j-1) = ISIs(j)+ISIs(j-1);
    CV2_all_diff(j-1) = ISIs(j)-ISIs(j-1);
    
    CV2_all(j-1) = 2*abs(CV2_all_diff(j-1))/CV2_all_sum(j-1) ;
    
    if ISIs(j)>upperLim | ISIs(j-1)>upperLim
        validISI(j-1)=0;
    else
        validISI(j-1)=1;
    end
end

if ignoreMode==0
    CV2 = mean( CV2_all );
else
    CV2 = mean( CV2_all(find(validISI==1)) );
end