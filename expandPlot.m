%   Helper function for plotNIRS
%
%   Author: Dino Soldic
%   Email: dino.soldic@urjc.es
%   Date: 2025-10-29
%
%   See also plotNIRS

function expandPlot(plotData, sig, groups, condition, channel, time, colors, isAvg, aType, isSubstracted, subCond)

    switch aType
        case "muit-condition"
            figure('Name', channel, 'NumberTitle', 'off');
            hold on;

            lineHandles = gobjects(1, length(groups));

            for grpIdx = 1:length(groups)

                if isAvg
                    chans = split(channel, ";");
                    dataToAverage = zeros([size(plotData.(groups{grpIdx}).(condition).(chans{1})) length(chans)]);

                    if isSubstracted

                        for c = 1:length(chans)
                            dataToAverage(:, :, c) = plotData.(groups{grpIdx}).(condition).(chans{c}) - plotData.(groups{grpIdx}).(subCond).(chans{c});
                        end

                    else

                        for c = 1:length(chans)
                            dataToAverage(:, :, c) = plotData.(groups{grpIdx}).(condition).(chans{c});
                        end

                    end

                    dataMat = mean(dataToAverage, 3);

                else

                    if isSubstracted
                        dataMat = plotData.(groups(grpIdx)).(condition).(channel) - plotData.(groups(grpIdx)).(subCond).(channel);
                    else
                        dataMat = plotData.(groups(grpIdx)).(condition).(channel);
                    end

                end

                % Compute mean and SEM across subjects
                y = mean(dataMat, 2);
                semY = std(dataMat, 0, 2) ./ sqrt(size(dataMat, 2));

                % Define shading area
                upper = y + semY;
                lower = y - semY;

                % Plot SEM as shaded region
                fill([time; flipud(time)], [upper; flipud(lower)], ...
                    colors(grpIdx, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');

                % Plot mean and store handle for legend
                lineHandles(grpIdx) = plot(time, y, 'Color', colors(grpIdx, :), 'LineWidth', 1.5);
            end

            % Add significance at bottom
            ylims = ylim(gca);
            sigY = ylims(1) * 1.01;
            plot(time(sig), repmat(sigY, sum(sig), 1), 'r.', 'MarkerSize', 10);

            set(gca, 'FontName', 'Times New Roman', 'FontSize', 30);
            xlabel('Time (s)', 'FontSize', 50);
            ylabel('Concentration (mol/L^{-1})', 'FontSize', 50);
            xlim([min(time) * 1.01, max(time) * 1.01]);
            ylim padded;
            xticks(floor(min(time)):ceil(max(time)));
            xline(0, '--', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.8);
            yline(0, '--', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.8);
            legend(lineHandles, groups, 'FontSize', 24, 'Location', 'northwest');
            title(sprintf('Channel: %s', channel));
            hold off;

        case "muit-group"
            figure('Name', channel, 'NumberTitle', 'off');
            hold on;

            lineHandles = gobjects(1, length(groups));

            for grpIdx = 1:length(groups)

                groupData = [];

                if isAvg
                    chans = split(channel, ";");
                    dataToAverage = [];

                    for c = 1:length(chans)
                        groupData = [];

                        for condIdx = 1:length(condition)

                            if isSubstracted
                                groupData = [groupData, plotData.(groups{grpIdx}).(condition{condIdx}).(chans{c}) - plotData.(groups{grpIdx}).(subCond).(chans{c})]; %#ok<*AGROW>
                            else
                                groupData = [groupData, plotData.(groups{grpIdx}).(condition{condIdx}).(chans{c})];
                            end

                        end

                        dataToAverage(:, :, c) = groupData;

                    end

                    dataMat = mean(dataToAverage, 3);

                else

                    for condIdx = 1:length(condition)

                        if isSubstracted
                            groupData = [groupData, plotData.(groups{grpIdx}).(condition{condIdx}).(channel) - plotData.(groups{grpIdx}).(subCond).(channel)];
                        else
                            groupData = [groupData, plotData.(groups{grpIdx}).(condition{condIdx}).(channel)];
                        end

                    end

                    dataMat = groupData;

                end

                % Compute mean and SEM across subjects
                y = mean(dataMat, 2);
                semY = std(dataMat, 0, 2) ./ sqrt(size(dataMat, 2));

                % Define shading area
                upper = y + semY;
                lower = y - semY;

                % Plot SEM as shaded region
                fill([time; flipud(time)], [upper; flipud(lower)], ...
                    colors(grpIdx, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');

                % Plot mean and store handle for legend
                lineHandles(grpIdx) = plot(time, y, 'Color', colors(grpIdx, :), 'LineWidth', 1.5);
            end

            % Add significance at bottom
            ylims = ylim(gca);
            sigY = ylims(1) * 1.01;
            plot(time(sig), repmat(sigY, sum(sig), 1), 'r.', 'MarkerSize', 10);

            set(gca, 'FontName', 'Times New Roman', 'FontSize', 30);
            xlabel('Time (s)', 'FontSize', 50);
            ylabel('Concentration (mol/L^{-1})', 'FontSize', 50);
            xlim([min(time) * 1.01, max(time) * 1.01]);
            ylim padded;
            xticks(floor(min(time)):ceil(max(time)));
            xline(0, '--', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.8);
            yline(0, '--', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.8);
            legend(lineHandles, groups, 'FontSize', 24, 'Location', 'northwest');
            title(sprintf('Channel: %s', channel), "FontSize", 50, "FontWeight", "bold");
            hold off;

        case "mudt-within"
            figure('Name', channel, 'NumberTitle', 'off');
            hold on;

            lineHandles = gobjects(1, length(condition));

            for condIdx = 1:length(condition)

                if isAvg
                    chans = split(channel, ";");
                    dataToAverage = zeros([size(plotData.(groups).(condition{condIdx}).(chans{1})), length(chans)]);

                    for c = 1:length(chans)

                        if isSubstracted
                            dataToAverage(:, :, c) = plotData.(groups).(condition{condIdx}).(chans{c}) - plotData.(groups).(subCond).(chans{c});
                        else
                            dataToAverage(:, :, c) = plotData.(groups).(condition{condIdx}).(chans{c});
                        end

                    end

                    dataMat = mean(dataToAverage, 3);

                else

                    if isSubstracted
                        dataMat = plotData.(groups).(condition{condIdx}).(channel) - plotData.(groups).(subCond).(channel);
                    else
                        dataMat = plotData.(groups).(condition{condIdx}).(channel);
                    end

                end

                % Compute mean and SEM across subjects
                y = mean(dataMat, 2);
                semY = std(dataMat, 0, 2) ./ sqrt(size(dataMat, 2));

                % Define shading area
                upper = y + semY;
                lower = y - semY;

                % Plot SEM as shaded region
                fill([time; flipud(time)], [upper; flipud(lower)], ...
                    colors(condIdx, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');

                % Plot mean and store handle for legend
                lineHandles(condIdx) = plot(time, y, 'Color', colors(condIdx, :), 'LineWidth', 1.5);
            end

            % Add significance at bottom
            ylims = ylim(gca);
            sigY = ylims(1) * 1.01;
            plot(time(sig), repmat(sigY, sum(sig), 1), 'r.', 'MarkerSize', 10);

            set(gca, 'FontName', 'Times New Roman', 'FontSize', 30);
            xlabel('Time (s)', 'FontSize', 50);
            ylabel('Concentration (mol/L^{-1})', 'FontSize', 50);
            xlim([min(time) * 1.01, max(time) * 1.01]);
            ylim padded;
            xticks(floor(min(time)):ceil(max(time)));
            xline(0, '--', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.8);
            yline(0, '--', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.8);
            legend(lineHandles, condition, 'FontSize', 24, 'Location', 'northwest');
            title(sprintf('Channel: %s', channel));
            hold off;

        case "mudt-nogroup"
            figure('Name', channel, 'NumberTitle', 'off');
            hold on;

            lineHandles = gobjects(1, length(condition));

            for condIdx = 1:length(condition)

                if isAvg
                    chans = split(channel, ";");
                    dataToAverage = zeros([size(plotData.(groups{1}).(condition{condIdx}).(chans{1}), 1), size(plotData.(groups{1}).(condition{condIdx}).(chans{1}), 2) + size(plotData.(groups{2}).(condition{condIdx}).(chans{1}), 2), length(chans)]);

                    for c = 1:length(chans)

                        if isSubstracted
                            dataToAverage(:, :, c) = [plotData.(groups{1}).(condition{condIdx}).(chans{c}) - plotData.(groups{1}).(subCond).(chans{c}), plotData.(groups{2}).(condition{condIdx}).(chans{c}) - plotData.(groups{2}).(subCond).(chans{c})];
                        else
                            dataToAverage(:, :, c) = [plotData.(groups{1}).(condition{condIdx}).(chans{c}), plotData.(groups{2}).(condition{condIdx}).(chans{c})];
                        end

                    end

                    dataMat = mean(dataToAverage, 3);

                else

                    if isSubstracted
                        dataMat = [plotData.(groups{1}).(condition{condIdx}).(channel) - plotData.(groups{1}).(subCond).(channel), plotData.(groups{2}).(condition{condIdx}).(channel) - plotData.(groups{2}).(subCond).(channel)];
                    else
                        dataMat = [plotData.(groups{1}).(condition{condIdx}).(channel), plotData.(groups{2}).(condition{condIdx}).(channel)];
                    end

                end

                % Compute mean and SEM across subjects
                y = mean(dataMat, 2);
                semY = std(dataMat, 0, 2) ./ sqrt(size(dataMat, 2));

                % Define shading area
                upper = y + semY;
                lower = y - semY;

                % Plot SEM as shaded region
                fill([time; flipud(time)], [upper; flipud(lower)], ...
                    colors(condIdx, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');

                % Plot mean and store handle for legend
                lineHandles(condIdx) = plot(time, y, 'Color', colors(condIdx, :), 'LineWidth', 1.5);
            end

            % Add significance at bottom
            ylims = ylim(gca);
            sigY = ylims(1) * 1.01;
            plot(time(sig), repmat(sigY, sum(sig), 1), 'r.', 'MarkerSize', 10);

            set(gca, 'FontName', 'Times New Roman', 'FontSize', 30);
            xlabel('Time (s)', 'FontSize', 50);
            ylabel('Concentration (mol/L^{-1})', 'FontSize', 50);
            xlim([min(time) * 1.01, max(time) * 1.01]);
            ylim padded;
            xticks(floor(min(time)):ceil(max(time)));
            xline(0, '--', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.8);
            yline(0, '--', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.8);
            legend(lineHandles, condition, 'FontSize', 24, 'Location', 'northwest');
            title(sprintf('Channel: %s', channel));
            hold off;
    end

end
