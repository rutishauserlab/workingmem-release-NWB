function [sig_cells_sb, areas_sb] = NWB_SB_plotCell_Screening(nwbAll,all_units,params)
%NWB_SB_PLOTCELL_SCREENING Plots example cells from Kyzar et al 2023
% Outputs results of calcSelective_SB

%% Speficy unit filter
unit_examples = [... % Specify by SBID, cellID
    % Kaminski et al 2017, Figure 3.a-b
    14, 14; ... % Halle Berry
    14, 35; ... % Tiger Woods
    % Kyzar et al 2023, Figure 5.a-b
    7, 9; ... % Cat
    16, 16 ... % Delicate Arch
    ];

%% Apply unit filter
if any(~ismember(unit_examples(:,1),[all_units.subject_id]))
    warning('All examples not loaded. See import range. [4 7 16 21]')
end

% Filtering units
unit_mat = [[all_units.subject_id]',[all_units.unit_id]'];
inc_units = ismember(unit_mat,unit_examples, 'rows');
all_units_filtered = all_units(inc_units);

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