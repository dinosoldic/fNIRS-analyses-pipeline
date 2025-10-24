%   Helper function for plotNIRS
%
%   Author: Dino Soldic
%   Email: dino.soldic@urjc.es
%   Date: 2025-10-24
%
%   See also plotNIRS

function expandPlot(plotData, sig, groups, condition, channel, time, colors, isAvg, aType)

    switch aType
        case "muit-group"
            figure('Name', channel, 'NumberTitle', 'off');
            hold on;

            lineHandles = gobjects(1, length(groups));

            for grpIdx = 1:length(groups)

                if isAvg
                    chans = split(channel, ";");
                    dataToAverage = zeros([size(plotData.(groups{grpIdx}).(condition).(chans{1})) length(chans)]);

                    for c = 1:length(chans)
                        dataToAverage(:, :, c) = plotData.(groups{grpIdx}).(condition).(chans{c});
                    end

                    dataMat = mean(dataToAverage, 3);

                else
                    dataMat = plotData.(groups(grpIdx)).(condition).(channel);
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

            xlabel('Time (s)');
            ylabel('Concentration');
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 28);
            xlim([min(time) * 1.01, max(time) * 1.01]);
            ylim padded;
            xticks(floor(min(time)):ceil(max(time)));
            xline(0, '--', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.8);
            yline(0, '--', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.8);
            legend(lineHandles, groups, 'FontSize', 24, 'Location', 'northwest');
            title(sprintf('Channel: %s', channel));
            hold off;

        case "muit-condition"
            figure('Name', channel, 'NumberTitle', 'off');
            hold on;

            lineHandles = gobjects(1, length(condition));

            for condIdx = 1:length(condition)

                if isAvg
                    chans = split(channel, ";");
                    dataToAverage = zeros([size(plotData.(groups{1}).(condition{condIdx}).(chans{1}), 1), size(plotData.(groups{1}).(condition{condIdx}).(chans{1}), 2) + size(plotData.(groups{2}).(condition{condIdx}).(chans{1}), 2), length(chans)]);

                    for c = 1:length(chans)
                        dataToAverage(:, :, c) = [plotData.(groups{1}).(condition{condIdx}).(chans{c}), plotData.(groups{2}).(condition{condIdx}).(chans{c})];
                    end

                    dataMat = mean(dataToAverage, 3);

                else
                    dataMat = [plotData.(groups{1}).(condition{condIdx}).(channel), plotData.(groups{2}).(condition{condIdx}).(channel)];
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

            xlabel('Time (s)');
            ylabel('Concentration');
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 28);
            xlim([min(time) * 1.01, max(time) * 1.01]);
            ylim padded;
            xticks(floor(min(time)):ceil(max(time)));
            xline(0, '--', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.8);
            yline(0, '--', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.8);
            legend(lineHandles, condition, 'FontSize', 24, 'Location', 'northwest');
            title(sprintf('Channel: %s', channel));
            hold off;

        case "mudt"
            figure('Name', channel, 'NumberTitle', 'off');
            hold on;

            lineHandles = gobjects(1, length(condition));

            for condIdx = 1:length(condition)

                if isAvg
                    chans = split(channel, ";");
                    dataToAverage = zeros([size(plotData.(groups).(condition{condIdx}).(chans{1})), length(chans)]);

                    for c = 1:length(chans)
                        dataToAverage(:, :, c) = plotData.(groups).(condition{condIdx}).(chans{c});
                    end

                    dataMat = mean(dataToAverage, 3);

                else
                    dataMat = plotData.(groups).(condition{condIdx}).(channel);
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

            xlabel('Time (s)');
            ylabel('Concentration');
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 28);
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
