function [all_units] =NWB_SB_extractUnits(nwbAll)
%NWB_SB_extractUnits Takes a cell array of loaded nwb files and returns all
%   single unit information across all the files. 
%   nwbAll: cell array of loaded nwb files
%   importRange:  
%
%   mkyzar 4/27/2023
all_units = {};
for i=1:length(nwbAll)
    unit_ids = num2cell(nwbAll{i}.units.id.data.load() + 1); % Convert to 1-based indexing
    session_count = num2cell(i.*ones(length(unit_ids),1));
    subject_id = num2cell(str2double(nwbAll{i}.general_subject.subject_id).*ones(length(unit_ids),1));
    electrodes_ind = num2cell(nwbAll{i}.units.electrodes.data.load() + 1); % Convert to 1-based indexing
    fprintf('Loading SUs: Subject ID %d (%d/%d)...',subject_id{1},i,length(nwbAll))
    

    % Organize spike times
    spike_times_session = nwbAll{i}.units.spike_times.data.load();
    spike_times_ind = nwbAll{i}.units.spike_times_index.data.load();
    
    spike_times_cells = cell(length(unit_ids),1);
    spike_times_cells{1} = spike_times_session(1:spike_times_ind(1));
    for j = 2:length(spike_times_cells)
        spike_times_cells{j} = spike_times_session(spike_times_ind(j-1)+1:spike_times_ind(j));
    end
    
    % Import Waveforms
    % wf_mean_all = nwbAll{i}.units.waveform_mean.data.load();
    wf_all = nwbAll{i}.units.waveforms.data.load();
    wf_ind = nwbAll{i}.units.waveforms_index.data.load();
    wf_ind_ind = nwbAll{i}.units.waveforms_index_index.data.load();

    % Compare num_waveforms to num_spikes
    if size(wf_all,2) ~= length(spike_times_session)
        error('Number of spikes does not equal number of waveforms')
    elseif size(wf_all,2) ~= max(wf_ind)
        error('Waveform indices exceed the number of waveforms.')
    end

    % Performs double indexing
    wf_cells = cell(length(unit_ids),1);
    wf_cells{1} = wf_all(:,1:wf_ind(wf_ind_ind(1)))';
    for j = 2:length(wf_cells)
        wf_cells{j} = wf_all(:,wf_ind(wf_ind_ind(j)-1)+1:wf_ind(wf_ind_ind(j)))';
    end


    session_units = horzcat(session_count,subject_id,unit_ids,electrodes_ind,spike_times_cells, wf_cells);
    all_units = vertcat(all_units,session_units); %#ok<AGROW>
    fprintf(' Loaded \n')
end
all_units = cell2struct(all_units,{'session_count','subject_id','unit_id','electrodes','spike_times','waveforms'},2); % << Follows this column title format.
end