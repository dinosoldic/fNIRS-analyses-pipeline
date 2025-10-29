%   Perform NIRS statistical analyses on task or rest (coming soon) data.
%
%   NIRSAnalysis() prompts the user to select a folder containing NIRS
%   data files. It loads all subject/condition data, allows the user to
%   select data type(s) and channels, choose analyses, and then computes
%   various Mass Univariate Permutation-based t-tests.
%
%   NIRSAnalysis(ALLDATA) uses a preloaded ALLDATA struct containing
%   fields:
%       ALLDATATASK - structured task data
%       ALLDATAREST - structured rest data
%       time        - time vector
%
%   The directory MUST be organized as follows:
%
%     Study
%     ├── Group1
%     │   ├── Condition1
%     │   │   ├── data1
%     │   │   └── dataN
%     │   └── Condition2
%     │       ├── data1
%     │       └── dataN
%     ├── Group2
%     │   ├── Condition1
%     │   │   ├── data1
%     │   │   └── dataN
%     │   └── Condition2
%     │       ├── data1
%     │       └── dataN
%     └── Group3
%         ├── Condition1
%         │   ├── data1
%         │   └── dataN
%         └── Condition2
%             ├── data1
%             └── dataN
%
%   When prompted the Study folder should be selected to load data from.
%
%   Inputs:
%       ALLDATA (optional) - Struct containing preloaded task/rest data
%
%   Outputs:
%       Results are saved to files in user-selected directories.
%       The following analysis types are supported:
%           1. Mass Univariate Independent T-test - condition
%           2. Mass Univariate Independent T-test - group
%           3. Mass Univariate Dependent T-test - within groups
%           4. Mass Univariate Dependent T-test - across conditions (no groups)
%           5. Average Mass Univariate Independent T-test - condition
%           6. Average Mass Univariate Independent T-test - group
%           7. Average Mass Univariate Dependent T-test - within groups
%           8. Average Mass Univariate Dependent T-test - across conditions (no groups)
%
%   Notes:
%       - Subtraction of a baseline condition is supported for all analysis types.
%       - Permutation testing (Monte Carlo, 5000 permutations by default) is used
%         to compute p-values for robust statistical inference.
%       - Parallel computing is leveraged if MATLAB Parallel Toolbox is available.
%
%   Example usage:
%       % Run analysis from scratch
%       NIRSAnalysis();
%
%       % Run analysis using preloaded data
%       load('ALLDATA.mat');
%       NIRSAnalysis(ALLDATA);
%
%   Author: Dino Soldic
%   Email: dino.soldic@urjc.es
%   Date: 2025-10-29
%
%   See also plotNIRS, exportNIRS

function NIRSAnalysis(ALLDATA)

    % load data
    if nargin < 1

        % get data path
        rootPath = uigetdir(pwd, "Select folder with NIRS data");
        if rootPath == 0, error("Operation Canceled"), end

        ALLDATATASK = struct();
        ALLDATAREST = struct();

        [ALLDATATASK, ALLDATAREST, time] = loadData(rootPath, ALLDATATASK, ALLDATAREST);

        disp("✅ Done loading data.");
    else
        ALLDATATASK = ALLDATA.ALLDATATASK;
        ALLDATAREST = ALLDATA.ALLDATAREST;
        time = ALLDATA.time;
    end

    % select task
    analDataOpt = questdlg('Do you wish to analyze rest or task data?', 'Data Selection', 'Task', 'Rest', 'Task');
    if strcmp(analDataOpt, 'Task'), AnalData = ALLDATATASK; isTask = true; elseif strcmp(analDataOpt, 'Rest'), AnalData = ALLDATAREST; isTask = false; else, disp('Operation canceled. Shutting down'); return, end

    % save results path
    saveDir = uigetdir(pwd, "Select folder where results will be saved");
    if saveDir == 0, saveDir = pwd; end

    saveDirCSV = [saveDir, '\csv'];
    if ~isfolder(saveDirCSV), mkdir(saveDirCSV); end

    % Extract labels
    groupNames = fieldnames(AnalData);
    condNames = fieldnames(AnalData.(groupNames{1}));
    dataHeaders = split(fieldnames(AnalData.(groupNames{1}).(condNames{1})), "_");

    dataTypes = unique(dataHeaders(:, 1));

    chanLabels = unique(strcat(dataHeaders(:, 2), '-', dataHeaders(:, 3)));
    splitLabels = split(chanLabels, '-');
    num1 = str2double(splitLabels(:, 1));
    num2 = str2double(splitLabels(:, 2));
    [~, sortIdx] = sortrows([num1, num2]);
    chanLabels = chanLabels(sortIdx); % sorted

    % Select data and chan
    [dataOpt, ~] = listdlg('ListString', dataTypes, 'PromptString', 'Select the data type to analyze:', 'SelectionMode', 'multiple');
    if isempty(dataOpt), disp('Operation canceled. Shutting down'); return, end

    [chanOpt, ~] = listdlg('ListString', chanLabels, 'PromptString', 'Select channels to analyze:', 'SelectionMode', 'multiple');
    if isempty(chanOpt), disp('Operation canceled. Shutting down'); return, end

    % Mix types and chans
    selectedDataTypes = dataTypes(dataOpt);
    selectedChannels = chanLabels(chanOpt);
    chanCombos = strings(length(dataOpt) * length(chanOpt), 1);
    idx = 1;

    for d = 1:length(dataOpt)

        for c = 1:length(chanOpt)
            chanCombos(idx) = selectedDataTypes(d) + "_" + replace(selectedChannels(c), "-", "_");
            idx = idx + 1;
        end

    end

    % Select analysis
    if isTask
        analTypes = {"Compare Conditions (between groups)", "Compare Groups (fuse conditions)", "Compare Conditions (within groups)", "Compare Conditions (fuse groups)", ...
                         "Avg. Compare Conditions (between groups)", "Avg. Compare Groups (fuse conditions)", "Avg. Compare Conditions (within groups)", "Avg. Compare Conditions (fuse groups)"};
        [analOpt, ~] = listdlg('ListString', analTypes, 'PromptString', 'Select analyses:', 'SelectionMode', 'multiple');
        if isempty(analOpt), disp('Operation canceled. Shutting down'); return, end

        [condOpt, ~] = listdlg('ListString', condNames, 'PromptString', {'Select condition for Mass', 'Univariate T test:'}, 'SelectionMode', 'multiple');
        if isempty(condOpt), disp('Operation canceled. Shutting down'); return, end
        condOpt = sort(condOpt); % force lowest to highest

        doSubstractCond = questdlg('Do you wish to substract any condition?', 'Condition Substraction', 'Yes', 'No', 'Yes');
        if strcmp(doSubstractCond, 'Yes'), doSubstractCond = true; else, doSubstractCond = false; end

        if doSubstractCond

            while true
                [subCondOpt, ~] = listdlg('ListString', condNames, 'PromptString', {'Select condition to substract', 'from:'}, 'SelectionMode', 'single');
                if ~isempty(subCondOpt), break, end
            end

        end

    end

    % amount of permutation for montecarlo p values (5000 good - Maris & Oostenveld 2007)
    nPerm = 5000;
    p = 0.05;

    try

        if isempty(gcp('nocreate'))
            myCluster = parcluster('local');
            maxWorkers = myCluster.NumWorkers;
            parpool('local', maxWorkers);
        end

        useParallel = true;
        disp("Parallel pool initialized with " + num2str(gcp().NumWorkers) + " workers.");
    catch ME
        warning(ME.identifier, 'Parallel pool could not be started. Running in serial mode.\nReason: %s', ME.message);
        useParallel = false;
    end

    %% run analyses
    % mass uni ind T (muit) - condition
    % Compares condition effects between groups
    if any(analOpt == 1)

        % preallocate
        nCond = length(condOpt);
        nChan = numel(chanCombos);
        nTime = length(time);

        tVals = zeros(nChan, nTime, nCond);
        pVals = zeros(nChan, nTime, nCond);
        hVals = zeros(nChan, nTime, nCond);

        for condIdx = 1:nCond
            cond = condNames{condOpt(condIdx)};

            for chanIdx = 1:nChan
                chan = chanCombos{chanIdx};

                if doSubstractCond
                    groupOneData = AnalData.(groupNames{1}).(cond).(chan) - AnalData.(groupNames{1}).(condNames{subCondOpt}).(chan);
                    groupTwoData = AnalData.(groupNames{2}).(cond).(chan) - AnalData.(groupNames{2}).(condNames{subCondOpt}).(chan);
                else
                    groupOneData = AnalData.(groupNames{1}).(cond).(chan);
                    groupTwoData = AnalData.(groupNames{2}).(cond).(chan);

                end

                for timeIdx = 1:nTime

                    % t-test across subjects at this time point
                    [~, ~, ~, stats] = ttest2(groupOneData(timeIdx, :), groupTwoData(timeIdx, :));

                    allData = [groupOneData(timeIdx, :), groupTwoData(timeIdx, :)];

                    nTmpts1 = length(groupOneData(timeIdx, :));
                    permT = zeros(nPerm, 1);

                    if useParallel % run in parallel if toolbox installed for performance

                        % compute permuted t tests
                        parfor permIdx = 1:nPerm

                            permDataIdx = randperm(length(allData));
                            permGroupOneData = allData(permDataIdx(1:nTmpts1));
                            permGroupTwoData = allData(permDataIdx(nTmpts1 + 1:end));

                            [~, ~, ~, permStats] = ttest2(permGroupOneData, permGroupTwoData);
                            permT(permIdx) = permStats.tstat;
                        end

                    else

                        for permIdx = 1:nPerm

                            permDataIdx = randperm(length(allData));
                            permGroupOneData = allData(permDataIdx(1:nTmpts1));
                            permGroupTwoData = allData(permDataIdx(nTmpts1 + 1:end));

                            [~, ~, ~, permStats] = ttest2(permGroupOneData, permGroupTwoData);
                            permT(permIdx) = permStats.tstat;
                        end

                    end

                    tVals(chanIdx, timeIdx, condIdx) = stats.tstat;
                    pVals(chanIdx, timeIdx, condIdx) = mean(abs(permT) >= abs(stats.tstat));
                    hVals(chanIdx, timeIdx, condIdx) = pVals(chanIdx, timeIdx, condIdx) < p;

                end

            end

            disp("Mass-Uni Independent T-test for condition " + cond + " done");

        end

        % Make results file and save
        results = struct();
        results.type.analysis = "muit-condition";
        results.type.data = string(analDataOpt);
        results.type.substracted = doSubstractCond;
        if doSubstractCond, results.type.substractedCond = condNames{subCondOpt}; end
        results.stats.t = tVals;
        results.stats.p = pVals;
        results.stats.h = hVals;
        results.group.code = [];
        results.group.name = [];
        results.condition.code = [];
        results.condition.name = [];
        results.channel.code = [];
        results.channel.name = [];
        results.time = time;

        for condIdx = 1:nCond
            results.condition.code = [results.condition.code; condIdx];
            results.condition.name = [results.condition.name; string(condNames{condOpt(condIdx)})];
        end

        for chanIdx = 1:nChan
            results.channel.code = [results.channel.code; chanIdx];
            results.channel.name = [results.channel.name; string(chanCombos{chanIdx})];
        end

        for grpIdx = 1:numel(groupNames)
            results.group.code = [results.group.code; grpIdx];
            results.group.name = [results.group.name; string(groupNames{grpIdx})];
        end

        % export & save
        exportNIRS(results, saveDirCSV, time, doSubstractCond);

        if doSubstractCond

            if isTask
                save(fullfile(saveDir, "results_task_condition_substracted.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_condition_substracted.mat"), "results");
            end

        else

            if isTask
                save(fullfile(saveDir, "results_task_condition.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_condition.mat"), "results");
            end

        end

    end

    % mass uni ind T (muit) - group
    % Ignores condition separation and checks whether there's an effect group-wise
    if any(analOpt == 2)

        % preallocate
        nCond = numel(condNames(condOpt));
        nChan = numel(chanCombos);
        nTime = length(time);

        tVals = zeros(nChan, nTime);
        pVals = zeros(nChan, nTime);
        hVals = zeros(nChan, nTime);

        for chanIdx = 1:nChan
            chan = chanCombos{chanIdx};
            groupOneData = [];
            groupTwoData = [];

            for condIdx = 1:nCond

                if doSubstractCond
                    groupOneData = [groupOneData, AnalData.(groupNames{1}).(condNames{condOpt(condIdx)}).(chan) - AnalData.(groupNames{1}).(condNames{subCondOpt}).(chan)]; %#ok<*AGROW>
                    groupTwoData = [groupTwoData, AnalData.(groupNames{2}).(condNames{condOpt(condIdx)}).(chan) - AnalData.(groupNames{2}).(condNames{subCondOpt}).(chan)];

                else
                    groupOneData = [groupOneData, AnalData.(groupNames{1}).(condNames{condOpt(condIdx)}).(chan)];
                    groupTwoData = [groupTwoData, AnalData.(groupNames{2}).(condNames{condOpt(condIdx)}).(chan)];
                end

            end

            for timeIdx = 1:nTime

                % t-test across subjects at this time point
                [~, ~, ~, stats] = ttest2(groupOneData(timeIdx, :), groupTwoData(timeIdx, :));

                allData = [groupOneData(timeIdx, :), groupTwoData(timeIdx, :)];

                nTmpts1 = length(groupOneData(timeIdx, :));
                permT = zeros(nPerm, 1);

                if useParallel % run in parallel if toolbox installed for performance

                    % compute permuted t tests
                    parfor permIdx = 1:nPerm

                        permDataIdx = randperm(length(allData));
                        permGroupOneData = allData(permDataIdx(1:nTmpts1));
                        permGroupTwoData = allData(permDataIdx(nTmpts1 + 1:end));

                        [~, ~, ~, permStats] = ttest2(permGroupOneData, permGroupTwoData);
                        permT(permIdx) = permStats.tstat;
                    end

                else

                    for permIdx = 1:nPerm

                        permDataIdx = randperm(length(allData));
                        permGroupOneData = allData(permDataIdx(1:nTmpts1));
                        permGroupTwoData = allData(permDataIdx(nTmpts1 + 1:end));

                        [~, ~, ~, permStats] = ttest2(permGroupOneData, permGroupTwoData);
                        permT(permIdx) = permStats.tstat;
                    end

                end

                tVals(chanIdx, timeIdx) = stats.tstat;
                pVals(chanIdx, timeIdx) = mean(abs(permT) >= abs(stats.tstat));
                hVals(chanIdx, timeIdx) = pVals(chanIdx, timeIdx) < p;

            end

        end

        labels = strjoin(condNames(condOpt), "_");

        disp("Mass-Uni Independent T-test for groups done");

        % Make results file and save
        results = struct();
        results.type.analysis = "muit-group";
        results.type.data = string(analDataOpt);
        results.type.substracted = doSubstractCond;
        if doSubstractCond, results.type.substractedCond = condNames{subCondOpt}; end
        results.stats.t = tVals;
        results.stats.p = pVals;
        results.stats.h = hVals;
        results.group.code = [];
        results.group.name = [];
        results.condition.code = [];
        results.condition.name = [];
        results.channel.code = [];
        results.channel.name = [];
        results.time = time;

        for grpIdx = 1:numel(groupNames)
            results.group.code = [results.group.code; grpIdx];
            results.group.name = [results.group.name; string(groupNames{grpIdx})];
        end

        results.condition.code = 1;
        results.condition.name = labels;

        for chanIdx = 1:nChan
            results.channel.code = [results.channel.code; chanIdx];
            results.channel.name = [results.channel.name; string(chanCombos{chanIdx})];
        end

        % export & save
        exportNIRS(results, saveDirCSV, time, doSubstractCond);

        if doSubstractCond

            if isTask
                save(fullfile(saveDir, "results_task_group_substracted.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_group_substracted.mat"), "results");
            end

        else

            if isTask
                save(fullfile(saveDir, "results_task_group.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_group.mat"), "results");
            end

        end

    end

    % mass uni dep T (mudt) - within groups
    % Compares within each group the effect of condition vs another
    if any(analOpt == 3)

        % preallocate
        nCond = numel(condNames(condOpt));
        if nCond == 2, nCondMax = 1; else, nCondMax = nCond; end
        nGrp = numel(groupNames);
        nChan = numel(chanCombos);
        nTime = length(time);

        tVals = zeros(nChan, nTime, nCondMax, nGrp);
        pVals = zeros(nChan, nTime, nCondMax, nGrp);
        hVals = zeros(nChan, nTime, nCondMax, nGrp);

        labels = strings(nCondMax, 1);

        for grpIdx = 1:nGrp
            group = groupNames{grpIdx};

            for condIdx = 1:nCondMax
                condIdx2 = condIdx + 1;
                if condIdx2 > nCond, condIdx2 = 1; end

                for chanIdx = 1:nChan
                    chan = chanCombos{chanIdx};

                    if doSubstractCond
                        groupOneData = AnalData.(group).(condNames{condOpt(condIdx)}).(chan) - AnalData.(group).(condNames{subCondOpt}).(chan);
                        groupTwoData = AnalData.(group).(condNames{condOpt(condIdx2)}).(chan) - AnalData.(group).(condNames{subCondOpt}).(chan);
                    else
                        groupOneData = AnalData.(group).(condNames{condOpt(condIdx)}).(chan);
                        groupTwoData = AnalData.(group).(condNames{condOpt(condIdx2)}).(chan);
                    end

                    labels(condIdx) = sprintf("%s_%s", condNames{condOpt(condIdx)}, condNames{condOpt(condIdx2)});

                    for timeIdx = 1:nTime

                        % t-test across subjects at this time point
                        [~, ~, ~, stats] = ttest(groupOneData(timeIdx, :), groupTwoData(timeIdx, :));

                        dataDiff = groupOneData(timeIdx, :) - groupTwoData(timeIdx, :); % compute diff for dep t (Winkler et al 2016)

                        permT = zeros(nPerm, 1);

                        if useParallel % run in parallel if toolbox installed for performance

                            % compute permuted t tests
                            parfor permIdx = 1:nPerm

                                flipSigns = (rand(size(dataDiff)) > 0.5) * 2 - 1; % *-1
                                permDiff = dataDiff .* flipSigns;

                                [~, ~, ~, permStats] = ttest(permDiff, 0);
                                permT(permIdx) = permStats.tstat;
                            end

                        else

                            for permIdx = 1:nPerm

                                flipSigns = (rand(size(dataDiff)) > 0.5) * 2 - 1; % *-1
                                permDiff = dataDiff .* flipSigns;

                                [~, ~, ~, permStats] = ttest(permDiff, 0);
                                permT(permIdx) = permStats.tstat;
                            end

                        end

                        tVals(chanIdx, timeIdx, condIdx, grpIdx) = stats.tstat;
                        pVals(chanIdx, timeIdx, condIdx, grpIdx) = mean(abs(permT) >= abs(stats.tstat));
                        hVals(chanIdx, timeIdx, condIdx, grpIdx) = pVals(chanIdx, timeIdx, condIdx, grpIdx) < p; % t-test across subjects at this time point

                    end

                end

                disp("Mass-Uni Dependent T-test for " + group + "-" + condNames{condOpt(condIdx)} + "_" + condNames{condOpt(condIdx2)} + " done");
            end

        end

        % Make results file and save
        results = struct();
        results.type.analysis = "mudt-within";
        results.type.data = string(analDataOpt);
        results.type.substracted = doSubstractCond;
        if doSubstractCond, results.type.substractedCond = condNames{subCondOpt}; end
        results.stats.t = tVals;
        results.stats.p = pVals;
        results.stats.h = hVals;
        results.group.code = [];
        results.group.name = [];
        results.condition.code = [];
        results.condition.name = [];
        results.channel.code = [];
        results.channel.name = [];
        results.time = time;

        for grpIdx = 1:numel(groupNames)
            results.group.code = [results.group.code; grpIdx];
            results.group.name = [results.group.name; string(groupNames{grpIdx})];
        end

        for condIdx = 1:nCondMax
            results.condition.code = [results.condition.code; condIdx];
            results.condition.name = labels;
        end

        for chanIdx = 1:nChan
            results.channel.code = [results.channel.code; chanIdx];
            results.channel.name = [results.channel.name; string(chanCombos{chanIdx})];
        end

        % export & save
        exportNIRS(results, saveDirCSV, time, doSubstractCond);

        if doSubstractCond

            if isTask
                save(fullfile(saveDir, "results_task_withinSub_substracted.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_withinSub_substracted.mat"), "results");
            end

        else

            if isTask
                save(fullfile(saveDir, "results_task_withinSub.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_withinSub.mat"), "results");
            end

        end

    end

    % mass uni dep T (mudt) - no groups
    % Ignores group separation and checks whether there's an effect condition-wise
    if any(analOpt == 4)

        % preallocate
        nCond = numel(condNames(condOpt));
        if nCond == 2, nCondMax = 1; else, nCondMax = nCond; end
        nChan = numel(chanCombos);
        nTime = length(time);

        tVals = zeros(nChan, nTime, nCondMax);
        pVals = zeros(nChan, nTime, nCondMax);
        hVals = zeros(nChan, nTime, nCondMax);

        labels = strings(nCondMax, 1);

        for condIdx = 1:nCondMax
            condIdx2 = condIdx + 1;
            if condIdx2 > nCond, condIdx2 = 1; end

            for chanIdx = 1:nChan
                chan = chanCombos{chanIdx};

                if doSubstractCond
                    groupOneData = [AnalData.(groupNames{1}).(condNames{condOpt(condIdx)}).(chan) - AnalData.(groupNames{1}).(condNames{subCondOpt}).(chan), AnalData.(groupNames{2}).(condNames{condOpt(condIdx)}).(chan) - AnalData.(groupNames{2}).(condNames{subCondOpt}).(chan)];
                    groupTwoData = [AnalData.(groupNames{1}).(condNames{condOpt(condIdx2)}).(chan) - AnalData.(groupNames{1}).(condNames{subCondOpt}).(chan), AnalData.(groupNames{2}).(condNames{condOpt(condIdx2)}).(chan) - AnalData.(groupNames{2}).(condNames{subCondOpt}).(chan)];

                else
                    groupOneData = [AnalData.(groupNames{1}).(condNames{condOpt(condIdx)}).(chan), AnalData.(groupNames{2}).(condNames{condOpt(condIdx)}).(chan)];
                    groupTwoData = [AnalData.(groupNames{1}).(condNames{condOpt(condIdx2)}).(chan), AnalData.(groupNames{2}).(condNames{condOpt(condIdx2)}).(chan)];
                end

                labels(condIdx) = sprintf("%s_%s", condNames{condOpt(condIdx)}, condNames{condOpt(condIdx2)});

                for timeIdx = 1:nTime

                    % t-test across subjects at this time point
                    [~, ~, ~, stats] = ttest(groupOneData(timeIdx, :), groupTwoData(timeIdx, :));

                    dataDiff = groupOneData(timeIdx, :) - groupTwoData(timeIdx, :); % compute diff for dep t (Winkler et al 2016)

                    permT = zeros(nPerm, 1);

                    if useParallel % run in parallel if toolbox installed for performance

                        % compute permuted t tests
                        parfor permIdx = 1:nPerm

                            flipSigns = (rand(size(dataDiff)) > 0.5) * 2 - 1; % *-1
                            permDiff = dataDiff .* flipSigns;

                            [~, ~, ~, permStats] = ttest(permDiff, 0);
                            permT(permIdx) = permStats.tstat;
                        end

                    else

                        for permIdx = 1:nPerm

                            flipSigns = (rand(size(dataDiff)) > 0.5) * 2 - 1; % *-1
                            permDiff = dataDiff .* flipSigns;

                            [~, ~, ~, permStats] = ttest(permDiff, 0);
                            permT(permIdx) = permStats.tstat;
                        end

                    end

                    tVals(chanIdx, timeIdx, condIdx) = stats.tstat;
                    pVals(chanIdx, timeIdx, condIdx) = mean(abs(permT) >= abs(stats.tstat));
                    hVals(chanIdx, timeIdx, condIdx) = pVals(chanIdx, timeIdx, condIdx) < p; % t-test across subjects at this time point

                end

            end

            disp("Mass-Uni Independent T-test for " + condNames{condOpt(condIdx)} + "_" + condNames{condOpt(condIdx2)} + " done");
        end

        % Make results file and save
        results = struct();
        results.type.analysis = "mudt-nogroup";
        results.type.data = string(analDataOpt);
        results.type.substracted = doSubstractCond;
        if doSubstractCond, results.type.substractedCond = condNames{subCondOpt}; end
        results.stats.t = tVals;
        results.stats.p = pVals;
        results.stats.h = hVals;
        results.condition.code = [];
        results.condition.name = [];
        results.channel.code = [];
        results.channel.name = [];
        results.time = time;

        for condIdx = 1:nCondMax
            results.condition.code = [results.condition.code; condIdx];
            results.condition.name = labels;
        end

        for chanIdx = 1:nChan
            results.channel.code = [results.channel.code; chanIdx];
            results.channel.name = [results.channel.name; string(chanCombos{chanIdx})];
        end

        % export & save
        exportNIRS(results, saveDirCSV, time, doSubstractCond);

        if doSubstractCond

            if isTask
                save(fullfile(saveDir, "results_task_nogroup_substracted.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_nogroup_substracted.mat"), "results");
            end

        else

            if isTask
                save(fullfile(saveDir, "results_task_nogroup.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_nogroup.mat"), "results");
            end

        end

    end

    %% Averages
    % AVERAGE mass uni ind T (muit) - condition
    % Compares condition effects between groups
    if any(analOpt == 5)

        % preallocate
        nCond = length(condOpt);
        nChan = numel(chanCombos);
        nTime = length(time);

        tVals = zeros(nTime, nCond);
        pVals = zeros(nTime, nCond);
        hVals = zeros(nTime, nCond);

        for condIdx = 1:nCond
            cond = condNames{condOpt(condIdx)};

            groupOneData = zeros([size(AnalData.(groupNames{1}).(cond).(chanCombos{1})) nChan]);
            groupTwoData = zeros([size(AnalData.(groupNames{2}).(cond).(chanCombos{1})) nChan]);

            for chanIdx = 1:nChan
                chan = chanCombos{chanIdx};

                groupOneData(:, :, chanIdx) = AnalData.(groupNames{1}).(cond).(chan);
                groupTwoData(:, :, chanIdx) = AnalData.(groupNames{2}).(cond).(chan);

            end

            avgGroupOneData = mean(groupOneData, 3);
            avgGroupTwoData = mean(groupTwoData, 3);

            if doSubstractCond
                groupOneDataSub = zeros([size(AnalData.(groupNames{1}).(condNames{subCondOpt}).(chanCombos{1})) nChan]);
                groupTwoDataSub = zeros([size(AnalData.(groupNames{2}).(condNames{subCondOpt}).(chanCombos{1})) nChan]);

                for chanIdx = 1:nChan
                    chan = chanCombos{chanIdx};

                    groupOneDataSub(:, :, chanIdx) = AnalData.(groupNames{1}).(condNames{subCondOpt}).(chan);
                    groupTwoDataSub(:, :, chanIdx) = AnalData.(groupNames{2}).(condNames{subCondOpt}).(chan);

                end

                avgGroupOneData = avgGroupOneData - mean(groupOneDataSub, 3);
                avgGroupTwoData = avgGroupTwoData - mean(groupTwoDataSub, 3);
            end

            for timeIdx = 1:nTime

                % t-test across subjects at this time point
                [~, ~, ~, stats] = ttest2(avgGroupOneData(timeIdx, :), avgGroupTwoData(timeIdx, :));

                allData = [avgGroupOneData(timeIdx, :), avgGroupTwoData(timeIdx, :)];

                nTmpts1 = length(avgGroupOneData(timeIdx, :));
                permT = zeros(nPerm, 1);

                if useParallel % run in parallel if toolbox installed for performance

                    % compute permuted t tests
                    parfor permIdx = 1:nPerm

                        permDataIdx = randperm(length(allData));
                        permGroupOneData = allData(permDataIdx(1:nTmpts1));
                        permGroupTwoData = allData(permDataIdx(nTmpts1 + 1:end));

                        [~, ~, ~, permStats] = ttest2(permGroupOneData, permGroupTwoData);
                        permT(permIdx) = permStats.tstat;
                    end

                else

                    for permIdx = 1:nPerm

                        permDataIdx = randperm(length(allData));
                        permGroupOneData = allData(permDataIdx(1:nTmpts1));
                        permGroupTwoData = allData(permDataIdx(nTmpts1 + 1:end));

                        [~, ~, ~, permStats] = ttest2(permGroupOneData, permGroupTwoData);
                        permT(permIdx) = permStats.tstat;
                    end

                end

                tVals(timeIdx, condIdx) = stats.tstat;
                pVals(timeIdx, condIdx) = mean(abs(permT) >= abs(stats.tstat));
                hVals(timeIdx, condIdx) = pVals(timeIdx, condIdx) < p;

            end

            disp("Average Mass-Uni Independent T-test for condition " + cond + " done");

        end

        % Make results file and save
        results = struct();
        results.type.analysis = "muit-condition-avg";
        results.type.data = string(analDataOpt);
        results.type.substracted = doSubstractCond;
        if doSubstractCond, results.type.substractedCond = condNames{subCondOpt}; end
        results.stats.t = tVals;
        results.stats.p = pVals;
        results.stats.h = hVals;
        results.group.code = [];
        results.group.name = [];
        results.condition.code = [];
        results.condition.name = [];
        results.channel.name = [];
        results.time = time;

        for grpIdx = 1:numel(groupNames)
            results.group.code = [results.group.code; grpIdx];
            results.group.name = [results.group.name; string(groupNames{grpIdx})];
        end

        for condIdx = 1:nCond
            results.condition.code = [results.condition.code; condIdx];
            results.condition.name = [results.condition.name; string(condNames{condOpt(condIdx)})];
        end

        results.channel.name = strjoin(chanCombos, ";");

        % export & save
        exportNIRS(results, saveDirCSV, time, doSubstractCond);

        if doSubstractCond

            if isTask
                save(fullfile(saveDir, "results_task_avg_condition_substracted.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_avg_condition_substracted.mat"), "results");
            end

        else

            if isTask
                save(fullfile(saveDir, "results_task_avg_condition.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_avg_condition.mat"), "results");
            end

        end

    end

    % AVERAGE mass uni ind T (muit) - group
    % Ignores condition separation and checks whether there's an effect group-wise
    if any(analOpt == 6)

        % preallocate
        nCond = numel(condNames(condOpt));
        nChan = numel(chanCombos);
        nTime = length(time);

        tVals = zeros(nTime, 1);
        pVals = zeros(nTime, 1);
        hVals = zeros(nTime, 1);

        groupOneData = [];
        groupTwoData = [];

        for chanIdx = 1:nChan
            chan = chanCombos{chanIdx};

            groupOneDataFused = [];
            groupTwoDataFused = [];

            for condIdx = 1:nCond

                if doSubstractCond

                    groupOneDataFused = [groupOneDataFused, AnalData.(groupNames{1}).(condNames{condOpt(condIdx)}).(chan) - AnalData.(groupNames{1}).(condNames{subCondOpt}).(chan)];
                    groupTwoDataFused = [groupTwoDataFused, AnalData.(groupNames{2}).(condNames{condOpt(condIdx)}).(chan) - AnalData.(groupNames{2}).(condNames{subCondOpt}).(chan)];

                else

                    groupOneDataFused = [groupOneDataFused, AnalData.(groupNames{1}).(condNames{condOpt(condIdx)}).(chan)];
                    groupTwoDataFused = [groupTwoDataFused, AnalData.(groupNames{2}).(condNames{condOpt(condIdx)}).(chan)];

                end

            end

            groupOneData(:, :, chanIdx) = groupOneDataFused;
            groupTwoData(:, :, chanIdx) = groupTwoDataFused;

        end

        avgGroupOneData = mean(groupOneData, 3);
        avgGroupTwoData = mean(groupTwoData, 3);

        for timeIdx = 1:nTime

            % t-test across subjects at this time point
            [~, ~, ~, stats] = ttest2(avgGroupOneData(timeIdx, :), avgGroupTwoData(timeIdx, :));

            allData = [avgGroupOneData(timeIdx, :), avgGroupTwoData(timeIdx, :)];

            nTmpts1 = length(avgGroupOneData(timeIdx, :));
            permT = zeros(nPerm, 1);

            if useParallel % run in parallel if toolbox installed for performance

                % compute permuted t tests
                parfor permIdx = 1:nPerm

                    permDataIdx = randperm(length(allData));
                    permGroupOneData = allData(permDataIdx(1:nTmpts1));
                    permGroupTwoData = allData(permDataIdx(nTmpts1 + 1:end));

                    [~, ~, ~, permStats] = ttest2(permGroupOneData, permGroupTwoData);
                    permT(permIdx) = permStats.tstat;
                end

            else

                for permIdx = 1:nPerm

                    permDataIdx = randperm(length(allData));
                    permGroupOneData = allData(permDataIdx(1:nTmpts1));
                    permGroupTwoData = allData(permDataIdx(nTmpts1 + 1:end));

                    [~, ~, ~, permStats] = ttest2(permGroupOneData, permGroupTwoData);
                    permT(permIdx) = permStats.tstat;
                end

            end

            tVals(timeIdx) = stats.tstat;
            pVals(timeIdx) = mean(abs(permT) >= abs(stats.tstat));
            hVals(timeIdx) = pVals(timeIdx) < p;

        end

        labels = string(strjoin(condNames(condOpt), "_"));

        disp("Average Mass-Uni Independent T-test for " + labels + " done");

        % Make results file and save
        results = struct();
        results.type.analysis = "muit-group-avg";
        results.type.data = string(analDataOpt);
        results.type.substracted = doSubstractCond;
        if doSubstractCond, results.type.substractedCond = condNames{subCondOpt}; end
        results.stats.t = tVals;
        results.stats.p = pVals;
        results.stats.h = hVals;
        results.group.code = [];
        results.group.name = [];
        results.condition.code = [];
        results.condition.name = [];
        results.channel.name = [];
        results.time = time;

        for grpIdx = 1:numel(groupNames)
            results.group.code = [results.group.code; grpIdx];
            results.group.name = [results.group.name; string(groupNames{grpIdx})];
        end

        results.condition.code = 1;
        results.condition.name = labels;

        results.channel.name = strjoin(chanCombos, ";");

        % export & save
        exportNIRS(results, saveDirCSV, time, doSubstractCond);

        if doSubstractCond

            if isTask
                save(fullfile(saveDir, "results_task_avg_group_substracted.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_avg_group_substracted.mat"), "results");
            end

        else

            if isTask
                save(fullfile(saveDir, "results_task_avg_group.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_avg_group.mat"), "results");
            end

        end

    end

    % AVERAGE mass uni dep T (mudt) - within groups
    % Compares within each group the effect of condition vs another
    if any(analOpt == 7)

        % preallocate
        nCond = numel(condNames(condOpt));
        if nCond == 2, nCondMax = 1; else, nCondMax = nCond; end
        nGrp = numel(groupNames);
        nChan = numel(chanCombos);
        nTime = length(time);

        tVals = zeros(nTime, nCondMax, nGrp);
        pVals = zeros(nTime, nCondMax, nGrp);
        hVals = zeros(nTime, nCondMax, nGrp);

        labels = strings(nCondMax, 1);

        for grpIdx = 1:nGrp
            group = groupNames{grpIdx};

            for condIdx = 1:nCondMax
                condIdx2 = condIdx + 1;
                if condIdx2 > nCond, condIdx2 = 1; end

                nRows = size(AnalData.(group).(condNames{condOpt(condIdx)}).(chanCombos{1}), 1);
                nCols = size(AnalData.(group).(condNames{condOpt(condIdx)}).(chanCombos{1}), 2);
                groupOneData = zeros(nRows, nCols, nChan);
                groupTwoData = zeros(nRows, nCols, nChan);

                if doSubstractCond

                    for chanIdx = 1:nChan
                        chan = chanCombos{chanIdx};

                        groupOneData(:, :, chanIdx) = AnalData.(group).(condNames{condOpt(condIdx)}).(chan) - AnalData.(group).(condNames{subCondOpt}).(chan);
                        groupTwoData(:, :, chanIdx) = AnalData.(group).(condNames{condOpt(condIdx2)}).(chan) - AnalData.(group).(condNames{subCondOpt}).(chan);

                    end

                else

                    for chanIdx = 1:nChan
                        chan = chanCombos{chanIdx};

                        groupOneData(:, :, chanIdx) = AnalData.(group).(condNames{condOpt(condIdx)}).(chan);
                        groupTwoData(:, :, chanIdx) = AnalData.(group).(condNames{condOpt(condIdx2)}).(chan);

                    end

                end

                avgGroupOneData = mean(groupOneData, 3);
                avgGroupTwoData = mean(groupTwoData, 3);

                labels(condIdx) = sprintf("%s_%s", condNames{condOpt(condIdx)}, condNames{condOpt(condIdx2)});

                for timeIdx = 1:nTime

                    % t-test across subjects at this time point
                    [~, ~, ~, stats] = ttest(avgGroupOneData(timeIdx, :), avgGroupTwoData(timeIdx, :));

                    dataDiff = avgGroupOneData(timeIdx, :) - avgGroupTwoData(timeIdx, :); % compute diff for dep t (Winkler et al 2016)

                    permT = zeros(nPerm, 1);

                    if useParallel % run in parallel if toolbox installed for performance

                        % compute permuted t tests
                        parfor permIdx = 1:nPerm

                            flipSigns = (rand(size(dataDiff)) > 0.5) * 2 - 1; % *-1
                            permDiff = dataDiff .* flipSigns;

                            [~, ~, ~, permStats] = ttest(permDiff, 0);
                            permT(permIdx) = permStats.tstat;
                        end

                    else

                        for permIdx = 1:nPerm

                            flipSigns = (rand(size(dataDiff)) > 0.5) * 2 - 1; % *-1
                            permDiff = dataDiff .* flipSigns;

                            [~, ~, ~, permStats] = ttest(permDiff, 0);
                            permT(permIdx) = permStats.tstat;
                        end

                    end

                    tVals(timeIdx, condIdx, grpIdx) = stats.tstat;
                    pVals(timeIdx, condIdx, grpIdx) = mean(abs(permT) >= abs(stats.tstat));
                    hVals(timeIdx, condIdx, grpIdx) = pVals(timeIdx, condIdx, grpIdx) < p; % t-test across subjects at this time point

                end

                disp("Average Mass-Uni Dependent T-test for " + group + "-" + condNames{condOpt(condIdx)} + "_" + condNames{condOpt(condIdx2)} + " done");
            end

        end

        % Make results file and save
        results = struct();
        results.type.analysis = "mudt-within-avg";
        results.type.data = string(analDataOpt);
        results.type.substracted = doSubstractCond;
        if doSubstractCond, results.type.substractedCond = condNames{subCondOpt}; end
        results.stats.t = tVals;
        results.stats.p = pVals;
        results.stats.h = hVals;
        results.group.code = [];
        results.group.name = [];
        results.condition.code = [];
        results.condition.name = [];
        results.channel.name = [];
        results.time = time;

        for grpIdx = 1:numel(groupNames)
            results.group.code = [results.group.code; grpIdx];
            results.group.name = [results.group.name; string(groupNames{grpIdx})];
        end

        for condIdx = 1:nCondMax
            results.condition.code = [results.condition.code; condIdx];
            results.condition.name = labels;
        end

        results.channel.name = strjoin(chanCombos, ";");

        % export & save
        exportNIRS(results, saveDirCSV, time, doSubstractCond);

        if doSubstractCond

            if isTask
                save(fullfile(saveDir, "results_task_avg_withinSub_substracted.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_avg_withinSub_substracted.mat"), "results");
            end

        else

            if isTask
                save(fullfile(saveDir, "results_task_avg_withinSub.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_avg_withinSub.mat"), "results");
            end

        end

    end

    % AVERAGE mass uni dep T (mudt) - no groups
    % Ignores group separation and checks whether there's an effect condition-wise
    if any(analOpt == 8)

        % preallocate
        nCond = numel(condNames(condOpt));
        if nCond == 2, nCondMax = 1; else, nCondMax = nCond; end
        nChan = numel(chanCombos);
        nTime = length(time);

        tVals = zeros(nTime, nCond);
        pVals = zeros(nTime, nCond);
        hVals = zeros(nTime, nCond);

        labels = strings(nCondMax, 1);

        for condIdx = 1:nCondMax
            condIdx2 = condIdx + 1;
            if condIdx2 > nCond, condIdx2 = 1; end

            nRows = size(AnalData.(groupNames{1}).(condNames{condOpt(condIdx)}).(chanCombos{1}), 1);
            nCols = size(AnalData.(groupNames{1}).(condNames{condOpt(condIdx)}).(chanCombos{1}), 2) + size(AnalData.(groupNames{2}).(condNames{condOpt(condIdx)}).(chanCombos{1}), 2);
            groupOneData = zeros(nRows, nCols, nChan);
            groupTwoData = zeros(nRows, nCols, nChan);

            if doSubstractCond

                for chanIdx = 1:nChan
                    chan = chanCombos{chanIdx};

                    groupOneData(:, :, chanIdx) = [AnalData.(groupNames{1}).(condNames{condOpt(condIdx)}).(chan) - AnalData.(groupNames{1}).(condNames{subCondOpt}).(chan), AnalData.(groupNames{2}).(condNames{condOpt(condIdx)}).(chan) - AnalData.(groupNames{2}).(condNames{subCondOpt}).(chan)];
                    groupTwoData(:, :, chanIdx) = [AnalData.(groupNames{1}).(condNames{condOpt(condIdx2)}).(chan) - AnalData.(groupNames{1}).(condNames{subCondOpt}).(chan), AnalData.(groupNames{2}).(condNames{condOpt(condIdx2)}).(chan) - AnalData.(groupNames{2}).(condNames{subCondOpt}).(chan)];

                end

            else

                for chanIdx = 1:nChan
                    chan = chanCombos{chanIdx};

                    groupOneData(:, :, chanIdx) = [AnalData.(groupNames{1}).(condNames{condOpt(condIdx)}).(chan), AnalData.(groupNames{2}).(condNames{condOpt(condIdx)}).(chan)];
                    groupTwoData(:, :, chanIdx) = [AnalData.(groupNames{1}).(condNames{condOpt(condIdx2)}).(chan), AnalData.(groupNames{2}).(condNames{condOpt(condIdx2)}).(chan)];

                end

            end

            avgGroupOneData = mean(groupOneData, 3);
            avgGroupTwoData = mean(groupTwoData, 3);

            labels(condIdx) = sprintf("%s_%s", condNames{condOpt(condIdx)}, condNames{condOpt(condIdx2)});

            for timeIdx = 1:nTime

                % t-test across subjects at this time point
                [~, ~, ~, stats] = ttest(avgGroupOneData(timeIdx, :), avgGroupTwoData(timeIdx, :));

                dataDiff = avgGroupOneData(timeIdx, :) - avgGroupTwoData(timeIdx, :); % compute diff for dep t (Winkler et al 2016)

                permT = zeros(nPerm, 1);

                if useParallel % run in parallel if toolbox installed for performance

                    % compute permuted t tests
                    parfor permIdx = 1:nPerm

                        flipSigns = (rand(size(dataDiff)) > 0.5) * 2 - 1; % *-1
                        permDiff = dataDiff .* flipSigns;

                        [~, ~, ~, permStats] = ttest(permDiff, 0);
                        permT(permIdx) = permStats.tstat;
                    end

                else

                    for permIdx = 1:nPerm

                        flipSigns = (rand(size(dataDiff)) > 0.5) * 2 - 1; % *-1
                        permDiff = dataDiff .* flipSigns;

                        [~, ~, ~, permStats] = ttest(permDiff, 0);
                        permT(permIdx) = permStats.tstat;
                    end

                end

                tVals(timeIdx, condIdx) = stats.tstat;
                pVals(timeIdx, condIdx) = mean(abs(permT) >= abs(stats.tstat));
                hVals(timeIdx, condIdx) = pVals(timeIdx, condIdx) < p; % t-test across subjects at this time point

            end

            disp("Average Mass-Uni Dependent T-test for " + condNames{condOpt(condIdx)} + "_" + condNames{condOpt(condIdx2)} + " done");
        end

        % Make results file and save
        results = struct();
        results.type.analysis = "mudt-nogroup-avg";
        results.type.data = string(analDataOpt);
        results.type.substracted = doSubstractCond;
        if doSubstractCond, results.type.substractedCond = condNames{subCondOpt}; end
        results.stats.t = tVals;
        results.stats.p = pVals;
        results.stats.h = hVals;
        results.condition.code = [];
        results.condition.name = [];
        results.channel.name = [];
        results.time = time;

        for condIdx = 1:nCondMax
            results.condition.code = [results.condition.code; condIdx];
            results.condition.name = labels;
        end

        results.channel.name = strjoin(chanCombos, ";");

        % export & save
        exportNIRS(results, saveDirCSV, time, doSubstractCond);

        if doSubstractCond

            if isTask
                save(fullfile(saveDir, "results_task_avg_nogroup_substracted.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_avg_nogroup_substracted.mat"), "results");
            end

        else

            if isTask
                save(fullfile(saveDir, "results_task_avg_nogroup.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_avg_nogroup.mat"), "results");
            end

        end

    end

    % save data if no save exists
    if nargin < 1
        ALLDATA.ALLDATATASK = ALLDATATASK;
        ALLDATA.ALLDATAREST = ALLDATAREST;
        ALLDATA.time = time;

        save(fullfile(saveDir, "ALLDATA.mat"), "ALLDATA");
    end

    % Display completion
    fprintf('\n\t-------Analyses completed successfully-------\n');
    fprintf('\n\t\t  /\\_/\\ \t  /\\_/\\ \n\t\t ( o.o )\t ( ^.^ )\n\t\t  > ^ <\t\t  > ^ <\n');

end
