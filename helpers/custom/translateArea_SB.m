%
% translate area string to areacode, or vice-versa
% 
function [areaCode, areaLabel] = translateArea_SB(areaLabel, areaCode, condenseFlag)
if nargin<3
    condenseFlag = 0;
end
if nargin<2
    areaCode=[];
end

if ~isempty(areaCode) && isempty(areaLabel) % Convert area codes to labels
    switch (areaCode)
        case 1
            areaLabel='hippocampus_right';
        case 2
            areaLabel= 'hippocampus_left';
        case 3
            areaLabel='amygdala_right';
        case 4
            areaLabel='amygdala_left';
        case 5
            areaLabel= 'dorsal_anterior_cingulate_cortex_right';
        case 6
            areaLabel= 'dorsal_anterior_cingulate_cortex_left';                
        case 7
            areaLabel= 'pre_supplementary_motor_area_right';
        case 8
            areaLabel= 'pre_supplementary_motor_area_left';
        case 12
            areaLabel = 'orbitofrontal_cortex_left';
        case 13
            areaLabel = 'orbitofrontal_cortex_right';
        case 0
            areaLabel= 'TBD';
        otherwise
            warning(['Area code cannot be mapped, assigning NA: ' areaCode]);
            areaLabel='NA';
    end
end
if isempty(areaCode) && ~isempty(areaLabel) && condenseFlag == 0 % Convert labels to area codes
    switch (areaLabel)
        case 'hippocampus_right'
            areaCode=1;
        case 'hippocampus_left'
            areaCode=2;
        case 'amygdala_right'
            areaCode=3;
        case 'amygdala_left'
            areaCode=4;
        case 'dorsal_anterior_cingulate_cortex_right'
            areaCode=5;
        case 'dorsal_anterior_cingulate_cortex_left'
            areaCode=6;
        case 'pre_supplementary_motor_area_right'
            areaCode=7;
        case 'pre_supplementary_motor_area_left'
            areaCode=8;
        case 'orbitofrontal_cortex_left'
            areaCode=12;
        case 'orbitofrontal_cortex_right'
            areaCode=13;
        otherwise
            warning(['Label cannot be mapped, assigning NaN: ' areaLabel]);
            areaCode=NaN;
    end
end
if isempty(areaCode) && ~isempty(areaLabel) && condenseFlag == 1 % Condense labels into one abbreviated region. 
    switch (areaLabel)
        case {'hippocampus_right','hippocampus_left'}
            areaLabel='Hipp';
        case {'amygdala_right','amygdala_left'}
            areaLabel= 'Amg';
        case {'dorsal_anterior_cingulate_cortex_right','dorsal_anterior_cingulate_cortex_left'}
            areaLabel='dACC';
        case {'pre_supplementary_motor_area_right','pre_supplementary_motor_area_left'}
            areaLabel='pre-SMA';
        case {'orbitofrontal_cortex_left','orbitofrontal_cortex_right'}
            areaLabel= 'OFC';
        otherwise
            warning(['Label cannot be mapped, assigning NaN: ' areaLabel]);
            areaLabel='NA';
    end
end

if isempty(areaCode) && isempty(areaLabel) % Null generator
    warning('Null generated')
   areaLabel='NA';
   areaCode=nan;
end
end