function [nwbAll, importLog] = NWB_importFromFolder_SBCAT(datasetPath, importRange)
%NWB_importFromFolder Returns a cell array containing all imported nwb
%files in the specified folder and import range. 
%   pathIn: Path to desired folder
%   importRange: Range of SB/SC IDs to import. 
%
%   mkyzar 4/15/2023
%% Compiling subject list
fs = filesep;

if isfile([datasetPath fs 'dandiset.yaml']) % Check if the dandiset is being accessed. 
    fileStruct = dir([datasetPath fs '**' fs '*.nwb']);
    fullFilePaths = cellfun(@(x,y) strcat(x,fs,y),{fileStruct.folder},{fileStruct.name},'UniformOutput',false)';
    fullFilePaths_sorted = sort_nat(fullFilePaths); % Alphanumeric sort (3rd Party).
    if ~isempty(importRange) % Access subset. 
        fileList = fullFilePaths_sorted(importRange);
    else % Defaults to all available subjects
        fileList = fullFilePaths_sorted; 
    end
else % Unformatted NWB data. All-in-one folder
    fileStruct = dir([datasetPath fs '*.nwb']);
    fullFilePaths = cellfun(@(x,y) strcat(x,fs,y),{fileStruct.folder},{fileStruct.name},'UniformOutput',false)';
    fullFilePaths_sorted = sort_nat(fullFilePaths); % Alphanumeric sort (3rd Party).
    if ~isempty(importRange) % Search directory for subset
        fileList = fullFilePaths_sorted(importRange); % Defaults to all values. 
    else % Defaults to all available subjects
        fileList = fullFilePaths_sorted; 
    end
end

if isempty(importRange)
    importRange = 1:length(fileList);
end

%% Import Loop 
% (Parallel import currently non-functional)
importLog = cell(length(fileList),1); % For import error logging. 
nwbAll = cell(length(fileList),1);
for i = 1:length(fileList) 
    filePath = [fileList{i}];
    fileName = split(filePath,fs); fileName = fileName(end);
    try
        if ~isfile(filePath)
            warning('Specified file does not exist: %s',fileList{i})
            continue
        end
        fprintf('(%d) Reading %s ... ',importRange(i),fileName{:})
        nwbAll{i} = nwbRead(filePath); % Importing nwb objects to cell array
        logOut = sprintf('Read Successful: %s\n',fileName{:});
        fprintf('%s',logOut)
        importLog{i} = logOut;
    catch e
        warning('Error found for: %s\n',filePath)
        disp(e.message)
        logOut = sprintf('Error found for: %s\n',filePath);
        importLog{i} = logOut;
    end
end

nwbAll = nwbAll(~cellfun('isempty',nwbAll)); % Removing empty cells in the case of an error
if length(nwbAll) ~= length(importRange)
    warning('Number of imported files not equal to import range')
end
end
