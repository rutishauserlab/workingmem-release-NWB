function [nwbAll, importLog] = NWB_importFromFolder(datasetPath, importRange, taskFlag)
%NWB_importFromFolder Returns a cell array containing all imported nwb
%files in the specified folder and import range. 
%   pathIn: Path to desired folder
%   importRange: Range of SB/SC IDs to import. 
%
%   mkyzar 4/15/2023
%% Compiling subject list
fs = filesep;
if nargin < 2
    session_type = 1; % Defaults to loading SB
else
    session_type = taskFlag;
end
if ~ismember(taskFlag,[1,2])
    error('Task not specified. Must be ''1'' (Sternberg) or ''2'' (Screening)')
end

if isfile([datasetPath fs 'dandiset.yaml']) % Check if the dandiset is being accessed. 
    fileStruct = dir([datasetPath fs 'sub-*']);
    folderNamesTotal = sort_nat({fileStruct.name}'); % Alphanumeric sort (3rd Party).
    if ~isempty(importRange)
        preStr1 = string(importRange);
        preStr2 = string(importRange); preStr2(:) = 'sub-';
        preStrFs = string(importRange); preStrFs(:) = fs;
        strSearch = strcat(preStr2,preStr1);
        folderNamesAppended = cellfun(@(x) strcat(x,preStrFs(1)),folderNamesTotal );
        folderList = folderNamesTotal(contains(folderNamesAppended,strcat(strSearch,preStrFs)))';
        
        preStrFileExt = string(importRange); preStrFileExt(:) = ['_ses-' num2str(session_type) '_ecephys+image.nwb'];
        
        fileList = strcat(folderList,preStrFs,folderList,preStrFileExt);
    else
        error('Please specify import range.')
    end
else % Unexported NWB data
    fileStruct = dir([datasetPath fs '*.nwb']);
    fileNamesTotal = sort_nat({fileStruct.name}'); % Alphanumeric sort (3rd Party).
    if ~isempty(importRange) % Search directory for subset
        % Create array of search strings
        preStr1 = string(importRange);
        preStr2 = string(importRange); preStr2(:) = 'ID_';
        preStr3 = string(importRange); preStr3(:) = '_';
        strSearch = strcat(preStr2,preStr1,preStr3);
        fileList = fileNamesTotal(contains(fileNamesTotal,strSearch)); 
    else % Defaults to all available subjects
        fileList = fileNamesTotal; 
    end
end

%% Import Loop 
% (Parallel Functionality is work in progress)

% Parallel is currently non-functional. Looks like data in AppData\Local\Temp\... is being co-accessed, leading to the conflict. 
% Might want to bring this up as a potential addition to matnwb. 

% % % Set the number of workers to use
% % num_workers = 4;
% % % Create a parallel pool with the specified number of workers
% % importPool = gcp('nocreate'); % If no pool, do not create new one.
% % if isempty(importPool)
% %     importPool = parpool(num_workers);
% % end
% % 
% % % Adding an instance of matnwb to each channel to eliminate file access conflicts. 
% % matnwb_dir = dir(fullfile(matnwbPath, '**\*.*'));  %get list of files and folders in any subfolder
% % matnwb_dir = matnwb_dir(~[matnwb_dir.isdir]);
% % matnwb_files = strcat({matnwb_dir.folder}', filesep ,{matnwb_dir.name}');
% % addAttachedFiles(importPool,matnwb_files)

importLog = cell(length(fileList),1); % For import error logging. 
nwbAll = cell(length(fileList),1);
for i = 1:length(fileList) % Imports in parallel to save time. 
    filePath = [datasetPath fs fileList{i}];
    try
        if ~isfile(filePath)
            warning('Specified file does not exist: %s',fileList{i})
            continue
        end
        fprintf('Reading %s ... ',fileList{i})
        nwbAll{i} = nwbRead(filePath); % Importing nwb objects to cell array
        logOut = sprintf('Read Successful: %s\n',fileList{i});
        fprintf('%s',logOut)
        importLog{i} = logOut;
    catch e
        warning('Error found for: %s\n',filePath)
        disp(e.message)
        logOut = sprintf('Error found for: %s\n',filePath);
        importLog{i} = logOut;
    end
end
% enddelete(importPool) % Closing the pool

nwbAll = nwbAll(~cellfun('isempty',nwbAll)); % Removing empty cells in the case of an error
if length(nwbAll) ~= length(importRange)
    warning('Number of imported files not equal to import range')
end
end
