function [sig_cells_sb, areas_sb] = NWB_SB_plotCell_Sternberg(nwbAll,all_units,params)
%NWB_SB_PLOTCELL_STERNBERG Plots example cells from Kyzar et al 2023
% Outputs results of calcSelective_SB

%% Speficy unit filter
unit_examples = [... % Specify by sub-id, ses-id, unit_id
    % See Daume et al. Fig 3a
    5, 1, 29; ... % (6) sub-5-ses-1 cell 29 % Category 2: Animals
    ];

%% Apply unit filter

% Getting unit ids
unit_sub_ids ={all_units.subject_id}'; unit_sub_ids = cellfun(@(x) str2double(x), unit_sub_ids,'UniformOutput',false); unit_sub_ids = cell2mat(unit_sub_ids);
unit_sub_ses_ids ={all_units.session_id}'; unit_sub_ses_ids = cellfun(@(x) str2double(x), unit_sub_ses_ids,'UniformOutput',false); unit_sub_ses_ids = cell2mat(unit_sub_ses_ids);
unit_sub_ses_cell_ids = double([all_units.unit_id]');

% Creating filter for unit examples
in_sub = ismember(unit_sub_ids,unit_examples(:,1));
in_ses = ismember(unit_sub_ses_ids,unit_examples(:,2));
in_unit = ismember(unit_sub_ses_cell_ids,unit_examples(:,3));
is_example = in_sub.*in_ses.*in_unit;

if sum(is_example) < size(unit_examples,1)
    warning('All examples not loaded. See import range. [6]')
elseif sum(is_example) > size(unit_examples,1)
    error('Too many neurons filtered. Manually diagnose.')
end

% Applying filter
all_units_filtered = all_units(logical(is_example));

% Updating session count in units object
temp_units = all_units_filtered;
unique_subjects = unique([temp_units.session_count]);
for i = 1:length(unique_subjects)
    find_uniques = find([temp_units.session_count] == unique_subjects(i));
    for j = find_uniques
        all_units_filtered(j).session_count = i ;
    end
end
clear temp_units 

% Filtering files
nwbAll_filtered = nwbAll(unique_subjects);

[sig_cells_sb, areas_sb] = NWB_calcSelective_SB(nwbAll_filtered,all_units_filtered,params);
end