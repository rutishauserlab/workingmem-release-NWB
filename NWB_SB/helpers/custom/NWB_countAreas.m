function [areaCount_all, labels_sig, areaCount_sig, labels_tot] = NWB_countAreas(allAreas,sig_cells)
%NWB_COUNTAREAS Provides a count of all selective cells in the following areas of
%interest:
% OFC
% dACC
% pre-SMA
% Amyg
% Hippo
% Inputs: Area labels for all Sternberg cells; A logical array specifiying
% if a cell was significant. 
% Outputs: Condensed area counts per area and the related area labels. 
sigAreas = allAreas(sig_cells>0);
[~, allCondensed] = cellfun(@(x) translateArea_SB(x,[],1),allAreas,'UniformOutput',false); % Condensing left and right to one area.
[~, sigCondensed] = cellfun(@(x) translateArea_SB(x,[],1),sigAreas,'UniformOutput',false); % Condensing left and right to one area.
% areaLabel, areaCode, condenseFlag
[labels_tot_unordered, ~, uniqueInds_all] = unique(allCondensed);
[labels_sig_unordered, ~, uniqueInds_sig] = unique(sigCondensed);

% Organize Areas
ind_OFC_tot = find(strcmp(labels_tot_unordered,'vmPFC'));
ind_OFC_sig = find(strcmp(labels_sig_unordered,'vmPFC'));
ind_ACC_tot = find(strcmp(labels_tot_unordered,'dACC'));
ind_ACC_sig = find(strcmp(labels_sig_unordered,'dACC'));
ind_SMA_tot = find(strcmp(labels_tot_unordered,'pre-SMA'));
ind_SMA_sig = find(strcmp(labels_sig_unordered,'pre-SMA'));
ind_AMYG_tot = find(strcmp(labels_tot_unordered,'Amg'));
ind_AMYG_sig = find(strcmp(labels_sig_unordered,'Amg'));
ind_HP_tot = find(strcmp(labels_tot_unordered,'Hipp'));
ind_HP_sig = find(strcmp(labels_sig_unordered,'Hipp'));

ind_label_tot = [ind_OFC_tot ind_ACC_tot ind_SMA_tot ind_AMYG_tot ind_HP_tot];  % Changing order to desired format 
ind_label_sig = [ind_OFC_sig ind_ACC_sig ind_SMA_sig ind_AMYG_sig ind_HP_sig];  % Changing order to desired format 
labels_tot = labels_tot_unordered(ind_label_tot); labels_tot = cellfun(@(x) strrep(x,'_',' '), labels_tot, 'UniformOutput', false);
labels_sig = labels_sig_unordered(ind_label_sig); labels_sig = cellfun(@(x) strrep(x,'_',' '), labels_sig, 'UniformOutput', false);

areaCount_all = histcounts(uniqueInds_all,'BinMethod','integers'); areaCount_all = areaCount_all(ind_label_tot); % Changing order to desired format
areaCount_sig = histcounts(uniqueInds_sig,'BinMethod','integers'); areaCount_sig = areaCount_sig(ind_label_sig); % Changing order to desired format



end

