%   Recursively loads NIRS data from a folder structure.
%
%   [taskData, restData, time] = loadData(currentFolder, taskData, restData)
%
%   This helper function is used by NIRSAnalysis to read raw text data
%   files, process numeric conversion (handling comma/dot decimal issues),
%   and separate task vs rest data into structured arrays. It can
%   recursively navigate subdirectories to aggregate data from multiple
%   subjects or conditions.
%
%   Inputs:
%       currentFolder - (string/char) Path to the folder containing .txt data files
%       taskData      - (struct) Existing task data struct (used for recursion)
%       restData      - (struct) Existing rest data struct (used for recursion)
%
%   Outputs:
%       taskData      - (struct) Aggregated task data structured as
%                       taskData.GroupName.ConditionName.ChannelName
%       restData      - (struct) Aggregated rest data structured as
%                       restData.GroupName.ConditionName.ChannelName
%       time          - (double array) Time vector extracted from first column
%
%   Details:
%       - Expects .txt files with tab-delimited values.
%       - First row indicates task/rest columns (1 = task, 0 = rest)
%       - Second row contains headers for the data columns.
%       - Handles comma vs dot decimal marks automatically.
%       - Converts all numeric data to doubles.
%       - Concatenates data across files for the same group/condition/channel.
%
%   Example usage:
%       [taskData, restData, time] = loadData("C:\NIRSdata", struct(), struct());
%
%   Author: Dino Soldic
%   Email: dino.soldic@urjc.es
%   Date: 2025-10-24
%
%   See also NIRSAnalysis

function [taskData, restData, time] = loadData(currentFolder, taskData, restData)

    % Get all subdirectories
    items = dir(currentFolder);
    subdirs = items([items.isdir]);
    subdirs = subdirs(~ismember({subdirs.name}, {'.', '..'}));

    if isempty(subdirs)
        files = dir(fullfile(currentFolder, '*.txt'));

        if ~isempty(files)

            for f = 1:numel(files)

                % Extract grp and cond
                parts = strsplit(currentFolder, filesep);
                groupName = replace(parts{end - 1}, "-", "");
                condName = replace(parts{end}, "-", "");

                % Load data
                rawData = readcell(fullfile(files(f).folder, files(f).name), 'Delimiter', '\t');

                % Convert everything from row 3 onward to string
                strData = string(rawData(3:end, :));

                % Replace commas with dots in all strings at once
                strData = strrep(strData, ",", ".");

                % Convert numeric columns to double
                numData = str2double(strData);

                % Extract time
                time = numData(:, 1);

                % Determine task/rest indices
                taskIdx = [rawData{1, :}] == 1;

                % Get headers for task columns
                headerRow = split(string(rawData(2, taskIdx)), " ");
                headers = replace(headerRow(:, :, 2), ",", "_");

                % Extract task and rest numeric data
                rawTask = numData(:, taskIdx);
                taskIdx(1) = true; % skip time column for rest
                rawRest = numData(:, ~taskIdx);

                % check fields
                if ~isfield(taskData, groupName)
                    taskData.(groupName) = struct();
                end

                if ~isfield(taskData.(groupName), condName)
                    taskData.(groupName).(condName) = struct();
                end

                if ~isfield(restData, groupName)
                    restData.(groupName) = struct();
                end

                if ~isfield(restData.(groupName), condName)
                    restData.(groupName).(condName) = struct();
                end

                % fill data
                for structIdx = 1:numel(headers)
                    % Initialize field if it doesn't exist
                    if ~isfield(taskData.(groupName).(condName), headers{structIdx})
                        taskData.(groupName).(condName).(headers{structIdx}) = [];
                    end

                    if ~isfield(restData.(groupName).(condName), headers{structIdx})
                        restData.(groupName).(condName).(headers{structIdx}) = [];
                    end

                    taskData.(groupName).(condName).(headers{structIdx}) = [taskData.(groupName).(condName).(headers{structIdx}), rawTask(:, structIdx)];
                    restData.(groupName).(condName).(headers{structIdx}) = [restData.(groupName).(condName).(headers{structIdx}), rawRest(:, structIdx)];
                end

                disp("Loaded: " + groupName + '-' + condName + '-' + files(f).name);
            end

        end

    else
        % Recurse into subdirectories
        for idx = 1:numel(subdirs)
            [taskData, restData, time] = loadData(fullfile(currentFolder, subdirs(idx).name), taskData, restData);
        end

    end

end
