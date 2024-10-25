%Q1
clc;
clear;

data_files = {'NewData1.mat', 'NewData2.mat', 'NewData3.mat', 'NewData4.mat'};
signals = cell(1, numel(data_files));


for i = 1:numel(data_files)
    data = load(data_files{i});
    signals{i} = data.EEG_Sig;
end

for i = 1:numel(signals)
    plotEEG(signals{i})
    title(['EEG Signal ', num2str(i)], 'Fontsize', 14, 'Interpreter', 'latex');
end




%%
% Q3
[F_1, W_1, ~] = COM2R(signals{1}, 32);
components_signal1 = W_1 * signals{1};

[F_3, W_3, ~] = COM2R(signals{3}, 32);
components_signal3 = W_3 * signals{3};


%% Q4
clc;
fs=250;
Electrodes = load('Electrodes.mat');

% Plot Independent Components for signal 1
plotEEG(components_signal1);
title('Independent Components for signal 1', 'Fontsize', 14, 'Interpreter', 'latex');

% Plot Independent Components for signal 3
plotEEG(components_signal3);
title('Independent Components for signal 3', 'Fontsize', 14, 'Interpreter', 'latex');

% Plot Power Spectral Density (PSD) for signal 1
figure;
for i = 1:21
    subplot(7, 3, i)
    [pxx, f] = pwelch(components_signal1(i, :), [], [], [], fs);
    plot(f, pxx, 'Linewidth', 1);
    title(['PSD of ', num2str(i), 'th Source in Signal 1'], 'FontSize', 14, 'Interpreter', 'latex');
    xlabel('Frequency (Hz)', 'Interpreter', 'latex');
    ylabel('PSD (power/Hz)', 'Interpreter', 'latex');
end

% Plot Power Spectral Density (PSD) for signal 3
figure;
for i = 1:21
    subplot(7, 3, i)
    [pxx, f] = pwelch(components_signal3(i, :), [], [], fs);
    plot(f, pxx, 'Linewidth', 1);
    title(['PSD of ', num2str(i), 'th Source in Signal 3'], 'FontSize', 14, 'Interpreter', 'latex');
    xlabel('Frequency (Hz)', 'Interpreter', 'latex');
    ylabel('PSD (power/Hz)', 'Interpreter', 'latex');
end

% Plot topographic maps for signal 1
figure;
for i = 1:21
    subplot(4, 6, i)
    plottopomap(Electrodes.Electrodes.X(:, 1), Electrodes.Electrodes.Y(:, 1), Electrodes.Electrodes.labels(1, :), F_1(:, i));
    title([num2str(i), 'th Source in Signal 1'], 'FontSize', 14, 'Interpreter', 'latex');
end

% Plot topographic maps for signal 3
figure;
for i = 1:21
    subplot(4, 6, i)
    plottopomap(Electrodes.Electrodes.X(:, 1), Electrodes.Electrodes.Y(:, 1), Electrodes.Electrodes.labels(1, :), F_3(:, i));
    title([num2str(i), 'th Source in Signal 3'], 'FontSize', 14, 'Interpreter', 'latex');
end


%% Q5
clc;
SelSources_1 = [2,3,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20,21];
SelSources_3 = [2,6,9,10,14,17,18,19,20,21];
signal_1_den = F_1(:,SelSources_1)*components_signal1(SelSources_1,:);
signal_3_den = F_3(:,SelSources_3)*components_signal3(SelSources_3,:);
%% Q6
clc;

plotEEG(signal_1_den);
title('Denoised signal 1','Fontsize',14,'Interpreter','latex');
plotEEG(signal_3_den);
title('Denoised signal 3','Fontsize',14,'Interpreter','latex');

%% Functions

function plotEEG(X)
    load('Electrodes.mat');
    offset = max(abs(X(:)));
    feq = 250;
    ElecName = Electrodes.labels;
    disp_eeg(X, offset, feq, ElecName);
end

function plottopomap(elocsX, elocsY, elabels, data)
    % Define XY points for interpolation
    interp_detail = 100;
    interpX = linspace(min(elocsX)-.2, max(elocsX)+.25, interp_detail);
    interpY = linspace(min(elocsY), max(elocsY), interp_detail);

    % Create 2D grid locations based on 1D inputs
    [gridX, gridY] = meshgrid(interpX, interpY);

    % Interpolate the data on a 2D grid
    interpFunction = TriScatteredInterp(elocsY, elocsX, data);
    topodata = interpFunction(gridX, gridY);

    % Plot map
    contourf(interpY, interpX, topodata);
    hold on;
    scatter(elocsY, elocsX, 10, 'ro', 'filled');
    for i = 1:length(elocsX)
        text(elocsY(i), elocsX(i), elabels(i));
    end
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
end
