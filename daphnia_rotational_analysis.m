clear all;

% ============================================================
%   Daphnia Cumulative Angular Dynamics Analysis
%
%   This script loads trajectory data for multiple Daphnia,
%   filters it using a physics-based velocity threshold (v*),
%   computes incremental angular changes, and integrates them
%   into cumulative rotation angles. It then averages the
%   cumulative angles across all individuals and plots results.
%
%   HOW IT FITS IN THE PIPELINE:
%     Step 1: Convert TREX .npz tracking to CSV
%     Step 2: (THIS SCRIPT) Calculate cumulative angular dynamics
%     Step 3: CW/CCW classification & handedness (per individual)
% ============================================================


% ============ USER-DEFINED PARAMETERS =======================

% Set the base name of your dataset (e.g., 'Individual_6-18')
base = 'Individual_6-18';

% Input directory containing the CSV files
% NOTE: change 'path' and 'folder' to your actual location
inputDir = fullfile('path', 'to', 'your', 'csv', 'folder', base, '_csv');

% Number of tracked Daphnia in the dataset
N = 42;

% ============================================================



%% ========= Initialize Storage ==============================

% Cell arrays to store each Daphnia’s trajectory and fps
X    = cell(1, N);
Y    = cell(1, N);
Time = cell(1, N);
FPS  = cell(1, N);

% Progress bar (file loading step)
h = waitbar(0, 'Loading data and calculating cumulative angles...');



%% ========= Load CSV Data ===================================

for i = 0:N-1
    % Construct filename for each Daphnia
    csvFile = fullfile(inputDir, sprintf('Individual_daphnia%d.csv', i));

    % Read CSV file
    data = readtable(csvFile);

    % Store in arrays
    X{i+1}    = data.X;
    Y{i+1}    = data.Y;
    Time{i+1} = data.Time;
    FPS{i+1}  = data.fps(1);   % fps may repeat, use first value

    % Update progress bar
    waitbar(i / N, h, sprintf('Loading data... (%d/%d)', i, N));
end

% Define reference time axis (from first Daphnia)
ref_time   = Time{1};
timeLength = numel(ref_time);

% Average fps across all datasets
fps = mean(cell2mat(FPS));

% Preallocate matrix for storing cumulative angle trajectories
cumulativeAnglesMatrix = NaN(N, timeLength);



%% ========= Compute v* and Cumulative Angles ================

% Reset progress bar
waitbar(0, h, 'Calculating v_star and cumulative angles...');

for i = 1:N
    % Clean NaNs/Infs from trajectory
    X_clean    = X{i};
    Y_clean    = Y{i};
    Time_clean = Time{i};

    validMask = isfinite(X_clean) & isfinite(Y_clean) & isfinite(Time_clean);
    X_clean    = X_clean(validMask);
    Y_clean    = Y_clean(validMask);
    Time_clean = Time_clean(validMask);

    % Skip if too little data
    if numel(X_clean) < 3
        warning("⚠️ Daphnia #%d has too little data. Skipping...", i);
        continue;
    end

    % Instantaneous velocity components
    Vx = diff(X_clean) ./ diff(Time_clean);
    Vy = diff(Y_clean) ./ diff(Time_clean);

    % Speeds
    speeds = sqrt(Vx.^2 + Vy.^2);
    speeds = speeds(isfinite(speeds));

    if isempty(speeds)
        warning("⚠️ Daphnia #%d has empty speed. Skipping...", i);
        continue;
    end

    % ---- Velocity threshold v* estimation ----
    binSize = 0.5; % cm/s bin size for histogram
    [C, E] = histcounts(speeds, max(1, ceil(range(speeds) / binSize)));
    PROB   = C / sum(C);

    binCenters = (E(1:end-1) + E(2:end)) / 2;
    validIdx   = PROB > 0;

    logProb          = log(PROB(validIdx));
    binCentersValid  = binCenters(validIdx);

    % Linear fit → slope = -1/v*
    fitCoeffs = polyfit(binCentersValid, logProb, 1);
    slope     = fitCoeffs(1);
    v_star    = 1 / abs(slope);

    % ---- Apply v* threshold ----
    V_total   = sqrt(Vx.^2 + Vy.^2);
    valid_idx = (V_total >= v_star) & isfinite(V_total);

    Vx_f   = Vx(valid_idx);
    Vy_f   = Vy(valid_idx);
    Time_f = Time_clean(2:end);
    Time_f = Time_f(valid_idx);

    % Skip if too little post-threshold data
    if numel(Vx_f) < 2
        warning("⚠️ Daphnia #%d has insufficient post-threshold velocity. Skipping...", i);
        continue;
    end

    % ---- Angular increments ----
    delta_phi = zeros(1, numel(Vx_f) - 1);
    for a = 1:numel(delta_phi)
        dVx = Vx_f(a+1) - Vx_f(a);
        dVy = Vy_f(a+1) - Vy_f(a);

        numr = (Vx_f(a) * dVy) - (Vy_f(a) * dVx);
        denr = Vx_f(a)^2 + Vy_f(a)^2;

        delta_phi(a) = numr / denr;
    end

    % Normalize by 2π and integrate cumulatively
    norm_delta_phi   = delta_phi ./ (2 * pi);
    cumulative_angle = cumsum(norm_delta_phi);

    % Interpolate onto reference time grid
    interp_time     = Time_f(1:end-1);
    cumulative_interp = interp1(interp_time, cumulative_angle, ref_time, 'linear', NaN);

    % Store trajectory
    cumulativeAnglesMatrix(i, :) = cumulative_interp';

    % Update progress bar
    waitbar(i / N, h, sprintf('Processed daphnia %d/%d', i, N));
end

close(h);



%% ========= Average Across All Daphnia ======================

% Average cumulative angle (ignoring NaNs)
averageCumulativeAngle = nanmean(cumulativeAnglesMatrix, 1);

% Count number of valid contributors per time point
validCounts = sum(~isnan(cumulativeAnglesMatrix), 1);

% Discard averages where fewer than minDaphniaCount contributed
minDaphniaCount = 20;
averageCumulativeAngle(validCounts < minDaphniaCount) = NaN;



%% ========= Plot Results ====================================

% Plot # of valid Daphnia at each time point
figure;
plot(ref_time, validCounts);
xlabel('Time (s)');
ylabel('# of Daphnia contributing');
title('Daphnia count at each time point');

% Plot average cumulative angle across all Daphnia
figure;
plot(ref_time, averageCumulativeAngle, 'LineWidth', 1.4);
xlabel('Time (s)');
ylabel('Average Cumulative Angle');
title('Average Cumulative Angle of All Daphnia Over Time');
grid on;
