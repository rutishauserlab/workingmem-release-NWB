function str_out = condenseAreas(str_in)
% Remaps strings of brain areas stored in the SB NWB format to a
% de-lateralized naming convension (e.g. 'amygdala_left' -> 'amygdala')
% Input: A string of the lateralized brain area.
% Output: A string of the de-lateralized brain area. 
str_mid = str_in;
switch (str_mid)
    case {'amygdala_left','amygdala_right'}
        str_mid = 'Amy';
    case {'dorsal_anterior_cingulate_cortex_left','dorsal_anterior_cingulate_cortex_right'}
        str_mid = 'dACC';
    case {'hippocampus_left','hippocampus_right'}
        str_mid = 'Hippo';
    case {'pre_supplementary_motor_area_left','pre_supplementary_motor_area_right'}
        str_mid = 'preSMA';
    case {'ventral_medial_prefrontal_cortex_left','ventral_medial_prefrontal_cortex_right'}
        str_mid = 'vmPFC';
    otherwise
        error('Area cannot be remapped. Are you using the right input?')
end
str_out = str_mid;
end