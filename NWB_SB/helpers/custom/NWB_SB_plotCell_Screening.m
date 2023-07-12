function [sig_cells_sc, areas_sc] = NWB_SB_plotCell_Screening(nwbAll,all_units,params)
%NWB_SB_PLOTCELL_SCREENING Plots example cells from Kyzar et al 2023
% Outputs results of calcSelective_SC

%% Speficy unit filter
unit_examples = [... % Specify by SBID, cellID
    % Kyzar et al 2023, Fig 5.a-e
    4, 16; ... % Shoe House
    7, 24; ... % Zebras
    15, 3; ... % Moon man
    16, 22; ... % Delicate Arch
    21, 15; ... % The Pike Ferris Wheel 
    % Kaminski et al 2017, SF 1
    7, 24; ... % Zebras
    12, 15; ... % Leonardo DiCaprio
    14, 69; ... % Kanye West
    16, 14; ... % iPhone
    16, 8 ... % Bill Gates
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


[sig_cells_sc, areas_sc] = NWB_calcSelective_SC(nwbAll_filtered,all_units_filtered,params);
end