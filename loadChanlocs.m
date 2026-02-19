function chanlocs = loadChanlocs

    try
        eeglab('nogui');
    catch
        errordlg("EEGLAB toolbox missing or not working", "EEGLAB Error")
        return
    end

    %% load standard channels
    [chansFile, chansDir] = uigetfile("*", "Channel Coordinates", "Standard_Channels.txt", "MultiSelect", "off");
    if chansFile == 0, error("Channel coordinates need to be selected"), end
    chans = readcell(fullfile(chansDir, chansFile));

    shChans = [chans{:, 3}] < 8; % remove short CH
    chans = chans(shChans, :);

    nSens = numel(chans(:, 2));

    labels = cellfun(@(x, y) [num2str(x), '-', num2str(y)], chans(:, 2), chans(:, 3), 'UniformOutput', false);

    %% Adapt nirs to eeg chans for topoplot
    % temp xyz to load
    fidPath = sprintf("%s/chanlocs.xyz", pwd);
    fid = fopen(fidPath, 'w');

    for i = 1:nSens
        fprintf(fid, "%d %.4f %.4f %.4f %s\n", i, chans{i, 4}, chans{i, 5}, chans{i, 6}, labels{i});
    end

    fclose(fid);

    % load and delete temp
    chanlocs = struct();
    chanlocs = pop_chanedit(chanlocs, 'load', fidPath);
    delete(fidPath);
end
