function [hVals, pVals] = clusterCorrection(tReal, tNull, adj, p_crit)
    % tReal: [nChan x nTime]
    % tNull: [nPerm x nTime x nChan] or [nPerm x nTime] (for 1D)
    % adj:   [nChan x nChan] adjacency matrix

    [nChan, nTime] = size(tReal);
    nPerm = size(tNull, 1);

    % Initialize outputs
    hVals = zeros(nChan, nTime);
    pVals = ones(nChan, nTime);

    % Define two tails to check: positive (1) and negative (-1)
    tails = [1, -1];

    for tail = tails
        % 1. Thresholding: Use a fixed threshold based on the null distribution
        thresh = prctile(abs(tNull(:)), 100 * (1 - p_crit));

        if tail == 1
            real_sig_mask = tReal > thresh;
        else
            real_sig_mask = tReal < -thresh;
        end

        % 2. Find Real Clusters (Helper handles 1D and 2D)
        [real_labels, num_real] = findClusters(real_sig_mask, adj);
        real_masses = zeros(1, num_real);

        for i = 1:num_real
            real_masses(i) = sum(abs(tReal(real_labels == i)));
        end

        % 3. Build Null Distribution for this tail
        null_max_masses = zeros(nPerm, 1);

        for pIdx = 1:nPerm
            % Handle 1D vs 2D tNull dimensions
            if nChan > 1
                perm_t = squeeze(tNull(pIdx, :, :))'; % [nChan x nTime]
            else
                % Ensure 1D null is a row vector [1 x nTime]
                perm_t = reshape(tNull(pIdx, :), 1, nTime);
            end

            if tail == 1
                perm_sig_mask = perm_t > thresh;
            else
                perm_sig_mask = perm_t < -thresh;
            end

            [perm_labels, num_perm] = findClusters(perm_sig_mask, adj);

            if num_perm > 0
                m = zeros(1, num_perm);

                for j = 1:num_perm
                    m(j) = sum(abs(perm_t(perm_labels == j)));
                end

                null_max_masses(pIdx) = max(m);
            else
                null_max_masses(pIdx) = 0;
            end

        end

        % 4. Map p-values for this tail
        for i = 1:num_real
            % Statistical p-value calculation
            p_val = (sum(null_max_masses >= real_masses(i)) + 1) / (nPerm + 1);

            if p_val <= p_crit
                idx = (real_labels == i);
                hVals(idx) = 1;
                pVals(idx) = min(pVals(idx), p_val);
            end

        end

    end

end

function [labels, num] = findClusters(binaryMask, adj)
    [nChan, nTime] = size(binaryMask);

    % Check if adj matches nChan. If not, treat as no-adj (1D temporal)
    useAdj = ~isempty(adj) && size(adj, 1) == nChan && nChan > 1;

    visited = false(nChan, nTime);
    labels = zeros(nChan, nTime);
    num = 0;

    for c = 1:nChan

        for t = 1:nTime

            if binaryMask(c, t) && ~visited(c, t)
                num = num + 1;
                queue = [c, t];
                visited(c, t) = true;

                while ~isempty(queue)
                    curr = queue(1, :);
                    queue(1, :) = [];
                    labels(curr(1), curr(2)) = num;

                    % Temporal Neighbors
                    for dt = [-1, 1]
                        nt = curr(2) + dt;

                        if nt >= 1 && nt <= nTime && binaryMask(curr(1), nt) && ~visited(curr(1), nt)
                            visited(curr(1), nt) = true;
                            queue = [queue; curr(1), nt];
                        end

                    end

                    % Spatial Neighbors - only if dimensions match
                    if useAdj
                        neighbors = find(adj(curr(1), :));

                        for nc = neighbors

                            if binaryMask(nc, curr(2)) && ~visited(nc, curr(2))
                                visited(nc, curr(2)) = true;
                                queue = [queue; nc, curr(2)];
                            end

                        end

                    end

                end

            end

        end

    end

end
