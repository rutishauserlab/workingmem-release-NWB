%
% For a set of axis (subplots), determine and set common axis range
%
% axList: list of axis, i.e. [ax1 ax2 ax3]   with ax1=subplot(...)
% axisToSet: 1 X, 2 Y
%
% Implemented following https://www.mathworks.com/matlabcentral/answers/164039-readjusting-y-limits-on-subplots
%
%urut/Oct22
function setCommonAxisRange( axList, axisToSet )

if axisToSet==1
    str='Xlim';
else
    str='Ylim';
end

Lims_all = cell2mat(get(axList,str));
Lims_new = [min(Lims_all(:,1)) max(Lims_all(:,2))];
set(axList, str, Lims_new)
