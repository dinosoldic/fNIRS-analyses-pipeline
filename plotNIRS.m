%   Plot NIRS analysis results.
%
%   plotNIRS() prompts the user to select a result file and the corresponding
%   NIRS data file. It then generates plots of the selected statistical
%   analysis (e.g., mass univariate t-tests) across channels, conditions, and groups.
%
%   plotNIRS(results, data) plots the provided 'results' struct using the
%   provided NIRS 'data' struct.
%
%   Inputs:
%       results - Struct containing statistical analysis results
%       data    - Struct containing NIRS data (ALLDATATASK, ALLDATAREST, time)
%
%   The function supports plotting for:
%       - Mass Univariate Independent T-tests (group)
%       - Mass Univariate Independent T-tests (condition)
%       - Mass Univariate Dependent T-tests (within subjects)
%       - Average analyses for all of the above
%
%   Features:
%       - Automatically detects whether analysis is averaged ('avg') or
%         per-channel.
%       - Plots multiple channels in subplots, with automatic layout.
%       - Adds significance markers from results.stats.h.
%       - Legends, axis labels, and reference lines (x=0, y=0) are included.
%       - Clickable subplots expand to larger views using expandPlot callback.
%
%   Example usage:
%       % Load result and data manually
%       load('results_task_group.mat');
%       load('ALLDATA.mat');
%       plotNIRS(results, ALLDATA);
%
%       % Prompt user to select files
%       plotNIRS();
%
%   Notes:
%       - expandPlot function must be on the path to enable interactive subplot
%         expansion.
%
%   Author: Dino Soldic
%   Email: dino.soldic@urjc.es
%   Date: 2025-10-24
%
%   See also NIRSAnalysis

function plotNIRS(results, data)

    if nargin < 1 || ~isstruct(results)
        [file, path] = uigetfile(".mat", "Select result file");
        if file == 0, error("Operation Canceled"); end
        load(fullfile(path, file), "results");
    end

    if nargin < 2 || ~isstruct(data)
        [file, path] = uigetfile(".mat", "Select data file");
        if file == 0, error("Operation Canceled"); end
        load(fullfile(path, file), "results");
        data = ALLDATA;
    end

    time = data.time;

    if strcmpi(results.type.data, "task")
        plotData = data.ALLDATATASK;
    elseif strcmpi(results.type.data, "rest")
        plotData = data.ALLDATAREST;
    else
        error("Data type missing")
    end

    if contains(results.type.analysis, 'avg')
        isAvg = true;
    else
        isAvg = false;
    end

    %% plots
    % plot for muit
    if strcmpi(results.type.analysis, "muit-group") || strcmpi(results.type.analysis, "muit-group-avg")

        groups = string(results.group.name);
        conditions = string(results.condition.name);
        chanLabels = string(results.channel.name);

        colors = lines(length(groups)); % colors

        for condIdx = 1:length(conditions)

            nChans = length(chanLabels);
            nRows = ceil(sqrt(nChans));
            nCols = ceil(nChans / nRows);

            figure('Name', sprintf('%s', conditions(condIdx)), 'NumberTitle', 'off');

            for chanIdx = 1:length(chanLabels)

                ax = subplot(nRows, nCols, chanIdx);
                hold(ax, 'on');

                plotHandles = gobjects(1, numel(groups)); % store line handles

                for grpIdx = 1:length(groups)

                    if isAvg
                        chans = split(results.channel.name, ";");
                        dataToAverage = zeros([size(plotData.(groups{grpIdx}).(conditions{condIdx}).(chans{1})) length(chans)]);

                        for c = 1:length(chans)
                            dataToAverage(:, :, c) = plotData.(groups{grpIdx}).(conditions{condIdx}).(chans{c});
                        end

                        y = mean(mean(dataToAverage, 3), 2);

                    else
                        y = mean(plotData.(groups(grpIdx)).(conditions(condIdx)).(chanLabels(chanIdx)), 2);
                    end

                    % Plot and store handle
                    plotHandles(grpIdx) = plot(time, y, 'Color', colors(grpIdx, :), 'LineWidth', 1.5);

                end

                % Add significance markers
                sigY = ax.YAxis.Limits(1) * 1.01;

                if isAvg
                    sig = logical(results.stats.h(:, condIdx));
                else
                    sig = logical(results.stats.h(chanIdx, :, condIdx));
                end

                plot(ax, time(sig == 1), repmat(sigY, sum(sig == 1), 1), 'r.', 'MarkerSize', 5);

                % Axis formatting
                set(ax, 'FontName', 'Times New Roman', 'FontSize', 8);
                xlim(ax, [min(time) * 1.01, max(time) * 1.01]);
                ylim(ax, 'padded');
                xlabel(ax, 'Time (s)');
                ylabel(ax, 'Concentration');
                xticks(ax, floor(min(time)):ceil(max(time)));
                xline(ax, 0, '--', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.8);
                yline(ax, 0, '--', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.8);
                legend(ax, plotHandles, groups, 'FontSize', 8);

                title(ax, chanLabels(chanIdx), "FontSize", 12);

                % Add callback to expand subplot on click
                set(ax, 'ButtonDownFcn', @(src, evt) expandPlot(plotData, sig, groups, conditions(condIdx), chanLabels(chanIdx), time, colors, isAvg, "muit-group"));
                children = get(ax, 'Children');

                for c = 1:length(children)
                    set(children(c), 'HitTest', 'off'); % important!
                end

            end

        end

    end

    if strcmpi(results.type.analysis, "muit-condition") || strcmpi(results.type.analysis, "muit-condition-avg")

        conditions = split(string(results.condition.name), "_")';
        chanLabels = string(results.channel.name);
        nCond = size(conditions, 1);
        nSubCond = size(conditions, 2);

        colors = lines(nSubCond); % colors

        groupNames = fieldnames(plotData);

        for condIdx = 1:nCond

            nChans = length(chanLabels);
            nRows = ceil(sqrt(nChans));
            nCols = ceil(nChans / nRows);

            figure('Name', sprintf('%s-%s', conditions(condIdx, 1), conditions(condIdx, 2)), 'NumberTitle', 'off');

            for chanIdx = 1:length(chanLabels)

                ax = subplot(nRows, nCols, chanIdx);
                hold(ax, 'on');

                plotHandles = gobjects(1, nSubCond); % store line handles

                for subCondIdx = 1:nSubCond

                    if isAvg
                        chans = split(results.channel.name, ";");
                        dataToAverage = zeros([size(plotData.(groupNames{1}).(conditions{condIdx, subCondIdx}).(chans{1}), 1), size(plotData.(groupNames{1}).(conditions{condIdx, subCondIdx}).(chans{1}), 2) + size(plotData.(groupNames{2}).(conditions{condIdx}).(chans{1}), 2), length(chans)]);

                        for c = 1:length(chans)
                            dataToAverage(:, :, c) = [plotData.(groupNames{1}).(conditions{condIdx, subCondIdx}).(chans{c}), plotData.(groupNames{2}).(conditions{condIdx, subCondIdx}).(chans{c})];
                        end

                        y = mean(mean(dataToAverage, 3), 2);

                    else
                        y = mean([plotData.(groupNames{1}).(conditions{condIdx, subCondIdx}).(chanLabels{chanIdx}), plotData.(groupNames{2}).(conditions{condIdx, subCondIdx}).(chanLabels{chanIdx})], 2);
                    end

                    % Plot and store handle
                    plotHandles(subCondIdx) = plot(time, y, 'Color', colors(subCondIdx, :), 'LineWidth', 1.5);
                end

                % Add significance markers
                sigY = ax.YAxis.Limits(1) * 1.01;

                if isAvg
                    sig = logical(results.stats.h(:, condIdx));
                else
                    sig = logical(results.stats.h(chanIdx, :, condIdx));
                end

                plot(ax, time(sig == 1), repmat(sigY, sum(sig == 1), 1), 'r.', 'MarkerSize', 5);

                % Axis formatting
                set(ax, 'FontName', 'Times New Roman', 'FontSize', 8);
                xlim(ax, [min(time) * 1.01, max(time) * 1.01]);
                ylim(ax, 'padded');
                xlabel(ax, 'Time (s)');
                ylabel(ax, 'Concentration');
                xticks(ax, floor(min(time)):ceil(max(time)));
                xline(ax, 0, '--', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.8);
                yline(ax, 0, '--', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.8);
                legend(ax, plotHandles, string(conditions(condIdx, :)), 'FontSize', 8);

                title(ax, chanLabels(chanIdx), "FontSize", 12);

                % Add callback to expand subplot on click
                set(ax, 'ButtonDownFcn', @(src, evt) expandPlot(plotData, sig, groupNames, conditions(condIdx, :), chanLabels(chanIdx), time, colors, isAvg, "muit-condition"));
                children = get(ax, 'Children');

                for c = 1:length(children)
                    set(children(c), 'HitTest', 'off'); % important!
                end

            end

        end

    end

    if strcmpi(results.type.analysis, "mudt") || strcmpi(results.type.analysis, "mudt-avg")

        groups = string(results.group.name);
        conditions = split(string(results.condition.name), "_")';
        chanLabels = string(results.channel.name);
        nCond = size(conditions, 1);
        nSubCond = size(conditions, 2);

        colors = lines(nSubCond); % colors

        for grpIdx = 1:length(groups)

            for condIdx = 1:nCond

                nChans = length(chanLabels);
                nRows = ceil(sqrt(nChans));
                nCols = ceil(nChans / nRows);

                figure('Name', sprintf('%s_%s-%s', groups{grpIdx}, conditions(condIdx, 1), conditions(condIdx, 2)), 'NumberTitle', 'off');

                for chanIdx = 1:length(chanLabels)

                    ax = subplot(nRows, nCols, chanIdx);
                    hold(ax, 'on');

                    plotHandles = gobjects(1, nSubCond); % store line handles

                    for subCondIdx = 1:nSubCond

                        if isAvg
                            chans = split(results.channel.name, ";");
                            dataToAverage = zeros([size(plotData.(groups{grpIdx}).(conditions{condIdx, subCondIdx}).(chans{1})) length(chans)]);

                            for c = 1:length(chans)
                                dataToAverage(:, :, c) = plotData.(groups{grpIdx}).(conditions{condIdx, subCondIdx}).(chans{c});
                            end

                            y = mean(mean(dataToAverage, 3), 2);

                        else
                            y = mean(plotData.(groups{grpIdx}).(conditions{condIdx, subCondIdx}).(chanLabels{chanIdx}), 2);
                        end

                        % Plot and store handle
                        plotHandles(subCondIdx) = plot(time, y, 'Color', colors(subCondIdx, :), 'LineWidth', 1.5);
                    end

                    % Add significance markers
                    sigY = ax.YAxis.Limits(1) * 1.01;

                    if isAvg
                        sig = logical(results.stats.h(:, condIdx, grpIdx));
                    else
                        sig = logical(results.stats.h(chanIdx, :, condIdx, grpIdx));
                    end

                    plot(ax, time(sig == 1), repmat(sigY, sum(sig == 1), 1), 'r.', 'MarkerSize', 5);

                    % Axis formatting
                    set(ax, 'FontName', 'Times New Roman', 'FontSize', 8);
                    xlim(ax, [min(time) * 1.01, max(time) * 1.01]);
                    ylim(ax, 'padded');
                    xlabel(ax, 'Time (s)');
                    ylabel(ax, 'Concentration');
                    xticks(ax, floor(min(time)):ceil(max(time)));
                    xline(ax, 0, '--', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.8);
                    yline(ax, 0, '--', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.8);
                    legend(ax, plotHandles, string(conditions(condIdx, :)), 'FontSize', 8);

                    title(ax, chanLabels(chanIdx), "FontSize", 12);

                    % Add callback to expand subplot on click
                    set(ax, 'ButtonDownFcn', @(src, evt) expandPlot(plotData, sig, groups{grpIdx}, conditions(condIdx, :), chanLabels(chanIdx), time, colors, isAvg, "mudt"));
                    children = get(ax, 'Children');

                    for c = 1:length(children)
                        set(children(c), 'HitTest', 'off'); % important!
                    end

                end

            end

        end

    end

end
