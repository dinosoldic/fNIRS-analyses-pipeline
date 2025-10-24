%   Export NIRS statistical results to CSV.
%
%   exportNIRS() prompts the user to select a result file and a folder to
%   save the exported CSV files. Optionally, prompts for time vector and
%   whether conditions have been subtracted.
%
%   Inputs:
%       results    - Struct containing statistical analysis results
%       savepath   - Folder path to save CSV files
%       time       - Vector of time points
%       isSubCond  - Logical, true if conditions have been subtracted
%
%   Example usage:
%       exportNIRS(); % prompts user for file/folder
%       exportNIRS(results, 'C:\Exports', results.time, true);
%
%   Author: Dino Soldic
%   Email: dino.soldic@urjc.es
%   Date: 2025-10-24
%
%   See also NIRSAnalysis

function exportNIRS(results, savepath, time, isSubCond)

    if nargin < 1 || ~isstruct(results)
        [file, path] = uigetfile(".mat", "Select result file");
        if file == 0, error("Operation Canceled"); end
        load(fullfile(path, file), "results");
    end

    if nargin < 2 || ~ischar(savepath)
        savepath = uigetdir("pwd", "Select folder to save data");
        if savepath == 0, error("Operation Canceled"); end
    end

    if nargin < 3 || ~isa(time, "double") || ~isvector(time)
        time = results.time;
    end

    if nargin < 4 || ~islogical(isSubCond)
        isSubCond = questdlg("Have conditions been substracted", "Condition State", "Yes", "No", "No");
        if strcmp(isSubCond, 'Yes'), isSubCond = true; else, isSubCond = false; end
    end

    % extract vars
    nChan = length(results.channel.name);
    nTime = length(time);

    dataHeaders = ["Task Type", "Analysis", "Group", "Condition", "Time", strings(1, nChan * 2)];

    for c = 1:nChan
        dataHeaders(5 + c * 2 - 1) = results.channel.name(c) + "_t";
        dataHeaders(5 + c * 2) = results.channel.name(c) + "_p";
    end

    isAvg = contains(results.type.analysis, "avg");

    displacement = 1;

    switch results.type.analysis
        case {"muit-group", "muit-group-avg"} % groups

            nCond = results.condition.code(end);
            group = strjoin(results.group.name, "-");

            dataCell = cell(nCond * nTime, nChan * 2 + 5); % cols = chans + time, task, grp, cond * p+t stats

            for condIdx = 1:nCond
                condition = results.condition.name(condIdx);

                for timeIdx = 1:nTime

                    dataCell{displacement, 1} = results.type.data;
                    dataCell{displacement, 2} = results.type.analysis;
                    dataCell{displacement, 3} = group;
                    dataCell{displacement, 4} = condition;
                    dataCell{displacement, 5} = time(timeIdx);

                    for chanIdx = 1:nChan

                        if isAvg
                            dataCell{displacement, 5 + 2 * chanIdx - 1} = results.stats.t(timeIdx, condIdx);
                            dataCell{displacement, 5 + 2 * chanIdx} = results.stats.p(timeIdx, condIdx);
                        else
                            dataCell{displacement, 5 + 2 * chanIdx - 1} = results.stats.t(chanIdx, timeIdx, condIdx);
                            dataCell{displacement, 5 + 2 * chanIdx} = results.stats.p(chanIdx, timeIdx, condIdx);
                        end

                    end

                    displacement = displacement + 1;

                end

            end

            dataT = cell2table(dataCell, "VariableNames", dataHeaders);

            if isAvg

                if isSubCond
                    filename = 'results_task_avg_group_substracted.csv';
                else
                    filename = 'results_task_avg_group.csv';
                end

            else

                if isSubCond
                    filename = 'results_task_group_substracted.csv';
                else
                    filename = 'results_task_group.csv';
                end

            end

            writetable(dataT, fullfile(savepath, filename));
            disp("Exported:'" + filename(1:end - 4) + "' to csv");

        case {"muit-condition", "muit-condition-avg"} % conditions

            nCond = results.condition.code(end);
            group = "";

            dataCell = cell(nCond * nTime, nChan * 2 + 5); % cols = chans + time, task, grp, cond * p+t stats

            for condIdx = 1:nCond
                condition = results.condition.name(condIdx);

                for timeIdx = 1:nTime

                    dataCell{displacement, 1} = results.type.data;
                    dataCell{displacement, 2} = results.type.analysis;
                    dataCell{displacement, 3} = group;
                    dataCell{displacement, 4} = condition;
                    dataCell{displacement, 5} = time(timeIdx);

                    for chanIdx = 1:nChan

                        if isAvg
                            dataCell{displacement, 5 + 2 * chanIdx - 1} = results.stats.t(timeIdx, condIdx);
                            dataCell{displacement, 5 + 2 * chanIdx} = results.stats.p(timeIdx, condIdx);
                        else
                            dataCell{displacement, 5 + 2 * chanIdx - 1} = results.stats.t(chanIdx, timeIdx, condIdx);
                            dataCell{displacement, 5 + 2 * chanIdx} = results.stats.p(chanIdx, timeIdx, condIdx);
                        end

                    end

                    displacement = displacement + 1;

                end

            end

            dataT = cell2table(dataCell, "VariableNames", dataHeaders);

            if isAvg

                if isSubCond
                    filename = 'results_task_avg_condition_substracted.csv';
                else
                    filename = 'results_task_avg_condition.csv';
                end

            else

                if isSubCond
                    filename = 'results_task_condition_substracted.csv';
                else
                    filename = 'results_task_condition.csv';
                end

            end

            writetable(dataT, fullfile(savepath, filename));
            disp("Exported:'" + filename(1:end - 4) + "' to csv");

        case {"mudt", "mudt-avg"} % within

            nGrp = results.group.code(end);
            nCond = results.condition.code(end);

            dataCell = cell(nGrp * nCond * nTime, nChan * 2 + 5); % cols = chans + time, task, grp, cond * p+t stats

            for grpIdx = 1:nGrp
                group = results.group.name(grpIdx);

                for condIdx = 1:nCond
                    condition = results.condition.name(condIdx);

                    for timeIdx = 1:nTime

                        dataCell{displacement, 1} = results.type.data;
                        dataCell{displacement, 2} = results.type.analysis;
                        dataCell{displacement, 3} = group;
                        dataCell{displacement, 4} = condition;
                        dataCell{displacement, 5} = time(timeIdx);

                        for chanIdx = 1:nChan

                            if isAvg
                                dataCell{displacement, 5 + 2 * chanIdx - 1} = results.stats.t(timeIdx, condIdx, grpIdx);
                                dataCell{displacement, 5 + 2 * chanIdx} = results.stats.p(timeIdx, condIdx, grpIdx);
                            else
                                dataCell{displacement, 5 + 2 * chanIdx - 1} = results.stats.t(chanIdx, timeIdx, condIdx, grpIdx);
                                dataCell{displacement, 5 + 2 * chanIdx} = results.stats.p(chanIdx, timeIdx, condIdx, grpIdx);
                            end

                        end

                        displacement = displacement + 1;

                    end

                end

            end

            dataT = cell2table(dataCell, "VariableNames", dataHeaders);

            if isAvg

                if isSubCond
                    filename = 'results_task_avg_withinSub_substracted.csv';
                else
                    filename = 'results_task_avg_withinSub.csv';
                end

            else

                if isSubCond
                    filename = 'results_task_withinSub_substracted.csv';
                else
                    filename = 'results_task_withinSub.csv';
                end

            end

            writetable(dataT, fullfile(savepath, filename));
            disp("Exported:'" + filename(1:end - 4) + "' to csv");

    end

end
