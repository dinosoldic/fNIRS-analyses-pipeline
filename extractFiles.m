%   Recursively organizes HRF files into subject/group/condition folders.
%   Prepares data structure compatible with loadData for NIRS analysis.
%   Not general-purpose; requires user-specific adjustments.
%
%   Author: Dino Soldic
%   Email: dino.soldic@urjc.es
%   Date: 2025-10-24
%

function extractFiles(parent_root, output_root)

    if nargin < 2
        parent_root = "";
        output_root = "";
    end

    % Get all first-level directories inside parent_root
    topDirs = dir(parent_root);
    topDirs = topDirs([topDirs.isdir]);
    topDirs = topDirs(~ismember({topDirs.name}, {'.', '..'}));

    % Loop over each group
    for i = 1:numel(topDirs)
        currentGroup = fullfile(parent_root, topDirs(i).name);
        processFolder(currentGroup, output_root);

    end

    disp("✅ Done scanning all groups.")
end

function processFolder(currentFolder, output_root)

    % Get all subdirectories
    items = dir(currentFolder);
    subdirs = items([items.isdir]);
    subdirs = subdirs(~ismember({subdirs.name}, {'.', '..'}));

    if isempty(subdirs)
        % Leaf folder: check for HFR files
        files = dir(fullfile(currentFolder, '*HRF*'));
        files = files(~contains({files.name}, 'HRFMean')); % exclude

        if ~isempty(files)

            for f = 1:numel(files)

                % Split the path into parts
                parts = strsplit(currentFolder, filesep);

                % Find the folders for name and cond
                subIdx = find(startsWith(parts, 'Sujeto') | startsWith(parts, 'FM'), 1, 'last');
                condIdx = find(contains(parts, 'back'), 1, 'last');

                % Extract the subject name
                subName = parts{subIdx};

                if contains(subName, 'FM')
                    group = 'Fibromyalgia';
                else
                    group = 'Healthy-Controls';
                end

                condName = parts{condIdx};

                if contains(condName, '0-')
                    cond = 'back0';
                elseif contains (condName, '1-')
                    cond = 'back1';
                elseif contains (condName, '2-')
                    cond = 'back2';
                end

                % Extract numeric part from subject name
                numPart = regexp(subName, '\d+', 'match', 'once'); % e.g. '1', '10'

                % Zero-pad if single digit
                numPart = sprintf('%02d', str2double(numPart)); % '01', '02', '10', etc.

                filename = ['sub-', numPart, '_', group, '_', cond, '.txt'];
                destFolder = fullfile(output_root, group, cond);

                if ~exist(destFolder, 'dir')
                    mkdir(destFolder)
                end

                copyfile(fullfile(currentFolder, files(f).name), fullfile(destFolder, filename));

                disp("Extracted: " + fullfile(currentFolder, files(f).name))
            end

        end

    else
        % Recurse into subdirectories
        for i = 1:numel(subdirs)
            processFolder(fullfile(currentFolder, subdirs(i).name), output_root);
        end

    end

end
