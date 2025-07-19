clear all;

% ============ USER-DEFINED PARAMETERS ============

% Set the base name of your dataset (e.g., 'Individual_6-18')
base = 'Individual_6-18';

% Set the input directory containing the CSV files
inputDir = fullfile('path', 'to', 'your', 'csv', 'folder', base, '_csv');  % <- Change this

% Number of tracked Daphnia (42 here)
N = 42;

% ================================================

% Initialize cell arrays to store the data
X = cell(1, N);
Y = cell(1, N);
Time = cell(1, N);
FPS = cell(1, N);

% Progress bar for file loading
h = waitbar(0, 'Loading data and calculating cumulative angles...');

% Load CSV files
for i = 0:N-1
    csvFile = fullfile(inputDir, sprintf('Individual_daphnia%d.csv', i));
    data = readtable(csvFile);
    X{i+1} = data.X;
    Y{i+1} = data.Y;
    Time{i+1} = data.Time;
    FPS{i+1} = data.fps(1);  % Use first value if repeated
    waitbar(i / N, h, sprintf('Loading data... (%d/%d)', i, N));
end

ref_time = Time{1};
timeLength = numel(ref_time);
fps = mean(cell2mat(FPS));
cumulativeAnglesMatrix = NaN(N, timeLength);

% Reset progress bar
waitbar(0, h, 'Calculating v_star and cumulative angles...');

% Loop through each Daphnia
for i = 1:N
    X_clean = X{i};
    Y_clean = Y{i};
    Time_clean = Time{i};

    validMask = isfinite(X_clean) & isfinite(Y_clean) & isfinite(Time_clean);
    X_clean = X_clean(validMask);
    Y_clean = Y_clean(validMask);
    Time_clean = Time_clean(validMask);

    if numel(X_clean) < 3
        warning("⚠️ Daphnia #%d has too little data. Skipping...", i);
        continue;
    end

    Vx = diff(X_clean) ./ diff(Time_clean);
    Vy = diff(Y_clean) ./ diff(Time_clean);
    speeds = sqrt(Vx.^2 + Vy.^2);
    speeds = speeds(isfinite(speeds));
    if isempty(speeds)
        warning("⚠️ Daphnia #%d has empty speed. Skipping...", i);
        continue;
    end

    binSize = 0.5;
    [C, E] = histcounts(speeds, max(1, ceil(range(speeds) / binSize)));
    PROB = C / sum(C);
    binCenters = (E(1:end-1) + E(2:end)) / 2;
    validIdx = PROB > 0;

    logProb = log(PROB(validIdx));
    binCentersValid = binCenters(validIdx);
    fitCoeffs = polyfit(binCentersValid, logProb, 1);
    slope = fitCoeffs(1);
    v_star = 1 / abs(slope);

    V_total = sqrt(Vx.^2 + Vy.^2);
    valid_idx = (V_total >= v_star) & isfinite(V_total);
    Vx_f = Vx(valid_idx);
    Vy_f = Vy(valid_idx);
    Time_f = Time_clean(2:end);
    Time_f = Time_f(valid_idx);

    if numel(Vx_f) < 2
        warning("⚠️ Daphnia #%d has insufficient post-threshold velocity. Skipping...", i);
        continue;
    end

    delta_phi = zeros(1, numel(Vx_f) - 1);
    for a = 1:numel(delta_phi)
        dVx = Vx_f(a+1) - Vx_f(a);
        dVy = Vy_f(a+1) - Vy_f(a);
        numr = (Vx_f(a) * dVy) - (Vy_f(a) * dVx);
        denr = Vx_f(a)^2 + Vy_f(a)^2;
        delta_phi(a) = numr / denr;
    end

    norm_delta_phi = delta_phi ./ (2 * pi);
    cumulative_angle = cumsum(norm_delta_phi);
    interp_time = Time_f(1:end-1);
    cumulative_interp = interp1(interp_time, cumulative_angle, ref_time, 'linear', NaN);
    cumulativeAnglesMatrix(i, :) = cumulative_interp';
    waitbar(i / N, h, sprintf('Processed daphnia %d/%d', i, N));
end

close(h);

averageCumulativeAngle = nanmean(cumulativeAnglesMatrix, 1);
validCounts = sum(~isnan(cumulativeAnglesMatrix), 1);
minDaphniaCount = 20;
averageCumulativeAngle(validCounts < minDaphniaCount) = NaN;

figure;
plot(ref_time, validCounts);
xlabel('Time (s)');
ylabel('# of Daphnia contributing');
title('Daphnia count at each time point');

figure;
plot(ref_time, averageCumulativeAngle, 'LineWidth', 1.4);
xlabel('Time (s)');
ylabel('Average Cumulative Angle');
title('Average Cumulative Angle of All Daphnia Over Time');
grid on;
