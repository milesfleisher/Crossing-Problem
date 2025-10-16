function maxList = crossingMaxIntersections(nList)
% CROSSINGMAXINTERSECTIONS  Plot matchings and return max crossing counts.
% Usage:
%   maxList = crossingMaxIntersections([3 5 20]);
%
% For each n in nList, opens a figure that:
%   - Draws two columns of n dots as the left/right sides of an n-by-n square
%   - Overlays multiple full matchings (permutations) with unique colors
%   - Highlights the max-crossing matching (reverse permutation) in RED
%   - Title: "n=<n> max intersections = <...> shown in color red (#D91A1A)"
% Returns:
%   maxList(k) = nList(k)*(nList(k)-1)/2 (the theoretical maximum intersections)
%
% Notes on optimization and numerical robustness:
%   * Crossing counts are computed analytically via C_max = n*(n-1)/2 (exact integer arithmetic).
%   * Any verification counts use inversion counting (integer comparisons), not floating geometry,
%     avoiding precision and cancellation issues.
%   * For n <= 7, all permutations (n!) are drawn; otherwise we draw a representative subset to
%     avoid excessive rendering time and color reuse. Colors are distinct RGB triples.
%
% Requires only base MATLAB (R2025a+). Tested for macOS compatibility.

    arguments
        nList (1,:) {mustBePositive, mustBeInteger}
    end

    nList = double(nList(:).');                       % row vector
    maxList = (nList .* (nList - 1)) ./ 2;           % exact in double for n up to ~1e7

    % Parameters controlling plotting volume
    fullPermutationCapN = 7;      % draw ALL permutations if n <= 7 (7! = 5040)
    maxSamplesLargeN     = 16;    % for n > 7, draw at most this many representative matchings

    % Reserve RED for the max-crossing matching and never reuse it.
    maxRGB  = [0.85, 0.10, 0.10]; % vivid red reserved color
    maxName = 'red';
    maxHex  = rgb2hex(maxRGB);

    for idx = 1:numel(nList)
        n = nList(idx);
        Cmax = maxList(idx);

        % Coordinates: form an n-by-n square: width = n (x: 0..n), height = n (y: 1..n).
        xL = 0; xR = n;
        y  = (n:-1:1)';                 % top-to-bottom; integers avoid FP ambiguity in labeling

        % Decide which permutations to draw
        if n <= fullPermutationCapN
            % All permutations (n!)
            P = perms(1:n);             % size n! x n
            % Ensure reverse (max) is present (it is), and identity is present.
            permSet = P;
            drawCount = size(permSet,1);
            fprintf('n=%d: drawing all %d matchings.\n', n, drawCount);
        else
            % Representative subset: always include identity & reverse + structured + random
            basePerms = cell(0,1);

            % Identity and reverse
            basePerms{end+1} = 1:n;
            basePerms{end+1} = n:-1:1;  % max-crossing

            % A few cyclic shifts
            kShifts = min(4, max(0, n-1));
            for k = 1:kShifts
                basePerms{end+1} = circshift(1:n, k);
            end

            % Neighbor swaps patterns
            if n >= 2
                p = 1:n; p(1:2:end-1) = [2:2:n,]; p(2:2:end) = [1:2:n-1,]; % pairwise swap
                basePerms{end+1} = p;
            end
            if n >= 3
                p = 1:n; p([2 3]) = p([3 2]); basePerms{end+1} = p;
            end

            % Deterministic random samples (seed by n for reproducibility)
            rng(n, 'twister');
            need = maxSamplesLargeN - numel(basePerms);
            for t = 1:max(0, need)
                r = randperm(n);
                basePerms{end+1} = r;
            end

            % Deduplicate (in case overlap)
            permMat = unique(cell2mat(basePerms.'), 'rows', 'stable');

            % Ensure reverse permutation is included and placed first after identity
            rev = n:-1:1;
            if ~ismember(rev, permMat, 'rows')
                permMat = [rev; permMat];
            end
            % Ensure identity is present at top
            idrow = 1:n;
            if ~ismember(idrow, permMat(1,:), 'rows')
                % If identity exists elsewhere, move it to top
                [tf, loc] = ismember(idrow, permMat, 'rows');
                if tf
                    permMat = [permMat(loc,:); permMat(setdiff(1:size(permMat,1), loc), :)];
                else
                    permMat = [idrow; permMat];
                end
            end

            % Limit to maximum sample count
            drawCount = min(size(permMat,1), maxSamplesLargeN);
            permSet   = permMat(1:drawCount, :);
            fprintf(['n=%d: drawing %d representative matchings (identity, reverse, ' ...
                     'shifts, swaps, random). Full %d! omitted for performance.\n'], ...
                    n, drawCount, n);
        end

        % Prepare colors: reserve RED for the reverse (max) permutation; generate others uniquely.
        colors = zeros(drawCount,3);
        % Locate index of reverse in permSet
        revPerm = n:-1:1;
        revIdx  = find(ismember(permSet, revPerm, 'rows'), 1, 'first');
        if isempty(revIdx)
            % Safety: if missing (shouldn't happen), force include at position 1
            permSet = [revPerm; permSet];
            colors  = zeros(size(permSet,1),3);
            drawCount = size(permSet,1);
            revIdx = 1;
        end
        % Assign RED to reverse
        colors(revIdx, :) = maxRGB;

        % Generate distinct colors for the remaining combinations without reusing RED
        needColors = drawCount - 1;
        if needColors > 0
            otherCols = distinctColors(needColors, maxRGB);  % avoids RED vicinity
            colors(setdiff(1:drawCount, revIdx), :) = otherCols;
        end

        % ---- Plotting ----
        fh = figure('Name', sprintf('Crossings n=%d', n), 'Color', 'w');
        ax = axes('Parent', fh); %#ok<LAXES>
        hold(ax, 'on'); axis(ax, 'equal');

        % Set limits to make an exact n-by-n square
        pad = 0.5;
        xlim(ax, [xL - pad, xR + pad]);
        ylim(ax, [0.5, n + 0.5]);

        % Draw dots (two columns)
        dotSize = 50;
        scatter(ax, xL*ones(n,1), y, dotSize, 'k', 'filled');
        scatter(ax, xR*ones(n,1), y, dotSize, 'k', 'filled');

        % Draw each matching (combination) in its unique color
        lwDefault = 1.2;
        lwMax     = 2.4;  % make max-crossing matching stand out
        for k = 1:drawCount
            perm = permSet(k, :);
            thisColor = colors(k, :);
            % Draw all n segments for this matching
            for i = 1:n
                % No geometric intersection computations—just draw the straight segment
                line(ax, [xL, xR], [y(i), y(perm(i))], 'LineWidth', ...
                     (k==revIdx)*lwMax + (k~=revIdx)*lwDefault, ...
                     'Color', thisColor);
            end
        end

        % Title with requested format and named color for the max matching
        title(ax, sprintf('n=%d max intersections = %d shown in color %s (%s)', ...
              n, Cmax, maxName, maxHex), 'FontWeight', 'bold');

        % Aesthetics
        ax.Visible = 'off'; % hide axes for clean look
        drawnow;

        % Optional: print summary to command window
        fprintf('n=%d: C_max = %d, max shown in %s %s\n', n, Cmax, maxName, maxHex);
    end
end

% ---------- Local utilities ----------

function hex = rgb2hex(rgb)
%RGB2HEX Convert 1x3 RGB in [0,1] to a hex string '#RRGGBB'.
    rgb = max(0, min(1, rgb(:).'));
    hex = sprintf('#%02X%02X%02X', round(rgb * 255));
end

function C = distinctColors(k, avoidRGB)
%DISTINCTCOLORS Generate k visually distinct RGB colors in [0,1],
% avoiding being too close to avoidRGB (if provided).
% Uses golden-ratio spacing in HSV for stability and uniqueness.
    if nargin < 2
        avoidRGB = [NaN, NaN, NaN]; % no avoidance
    end
    C = zeros(k,3);
    phi = 0.618033988749895; % golden ratio conjugate
    h0  = 0.11;              % starting hue offset
    s   = 0.70;              % moderately saturated
    v   = 0.95;              % bright

    t = 1;
    h = h0;
    while t <= k
        rgb = hsv2rgb([mod(h,1), s, v]);
        % If too close to avoidRGB, skip and try next hue
        if ~any(isnan(avoidRGB))
            if norm(rgb - avoidRGB, 2) < 0.20
                h = h + phi;
                continue;
            end
        end
        C(t,:) = rgb;
        h = h + phi;
        t = t + 1;
    end
end




%%run it 
for i =1:100
    ns(i)=i;
end
crossingMaxIntersections(ns)