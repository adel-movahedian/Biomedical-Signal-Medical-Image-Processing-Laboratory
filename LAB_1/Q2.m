%% Q1 _1
ECG_data = load('ECG_sig.mat');
signal = ECG_data.Sig;
fs = ECG_data.sfreq;
t = 0:1/fs:((length(signal)-1))/fs;
clc; close all;
subplot(2,1,1)
plot(t,signal(:,1))
xlabel("Time(s)",'Interpreter','latex')
ylabel("Amplitude(mV)",'Interpreter','latex')
title('Channel One Signal','Interpreter','latex');
xlim([0 max(t)])
subplot(2,1,2)
plot(t,signal(:,2))
xlabel("Time(s)",'Interpreter','latex')
ylabel("Amplitude(mV)",'Interpreter','latex')
title('Channel Two Signal','Interpreter','latex');
xlim([0 max(t)])
%% Q1 - 2
time_1=[22*fs, 26*fs];
heartbeat_1 = signal(time_1(1):time_1(2),1);
time_2=[148*fs, 152*fs];
heartbeat_2 = signal(time_2(1):time_2(2),1);
subplot(2,1,1)
plot(t(time_1(1):time_1(2)),heartbeat_1);
xlabel("Time(s)",'Interpreter','latex')
ylabel("Amplitude(mV)",'Interpreter','latex')
title('Channel One Signal 22s-26s','Interpreter','latex');
subplot(2,1,2)
plot(t(time_2(1):time_2(2)),heartbeat_2);
xlabel("Time(s)",'Interpreter','latex')
ylabel("Amplitude(mV)",'Interpreter','latex')
title('Channel One Signal 148s - 152s','Interpreter','latex');

%% Q1 - 3
time_1=[42*fs, 44*fs];
heartbeat_1 = signal(time_1(1):time_1(2),1);
heartbeat_2 = signal(time_1(1):time_1(2),2);

subplot(2,1,1)
plot(t(time_1(1):time_1(2)),heartbeat_1,'red');
xlabel("Time(s)",'Interpreter','latex')
ylabel("Amplitude(mV)",'Interpreter','latex')
title('Channel One Signal 42s-44s','Interpreter','latex');

subplot(2,1,2)
plot(t(time_1(1):time_1(2)),heartbeat_2,'red');
xlabel("Time(s)",'Interpreter','latex')
ylabel("Amplitude(mV)",'Interpreter','latex')
title('Channel Two Signal 42s-44s','Interpreter','latex');

%% Q1 - 4
time_1=[round(143.6*fs), round(144.6*fs)];
heartbeat_1 = signal(time_1(1):time_1(2),1);

time_2=[round(23.9*fs), round(24.6*fs)];
heartbeat_2 = signal(time_2(1):time_2(2),1);

subplot(2,1,1)
plot(t(time_1(1):time_1(2)),heartbeat_1);
xlabel("Time(s)",'Interpreter','latex')
ylabel("Amplitude(mV)",'Interpreter','latex')
title('Channel One Signal 143.6s-144.6s','Interpreter','latex');
xlim([143.6 144.6])
hold on
x = [144.206,144.144,144.011, 144.278,144.439];
y = [0.875,-0.355,-0.27,-0.35,-0.15];
plot(x,y,'O');

subplot(2,1,2)
plot(t(time_2(1):time_2(2)),heartbeat_2);
xlabel("Time(s)",'Interpreter','latex')
ylabel("Amplitude(mV)",'Interpreter','latex')
title('Channel One Signal 23.9s - 24.4s','Interpreter','latex');
xlim([23.9 24.4])
hold on
x = [23.99,24.017,24.075,24.18,24.3];
y = [-0.14,-0.197,0.89,-0.289,-0.07];
plot(x,y,'O');

%% Q2 
R_points = ECG_data.ATRTIMED; 
R_numbers = ECG_data.ANNOTD; 
anomalyLabels = [        
        "NOTQRS", "LBBB", "RBBB", "ABERR", "PVC", "FUSION", "NPC", ...
        "APC", "SVPB", "VESC", "NESC", "PACE", "UNKNOWN", "NOISE", "", "ARFCT", ...
        "", "STCH", "TCH", "SYSTOLE", "DIASTOLE", "NOTE", "MEASURE", "PWAVE", "BBB", ...
        "PACESP", "TWAVE", "RHYTHM", "UWAVE", "LEARN", "FLWAV", "VFON", "VFOFF", ...
        "AESC", "SVESC", "LINK", "NAPC", "PFUS", "WFON", "WFOFF", "RONT"];
labelMap = [
        "NOTQRS", "NORMAL", "LBBB", "RBBB", "ABERR", "PVC", "FUSION", "NPC", ...
        "APC", "SVPB", "VESC", "NESC", "PACE", "UNKNOWN", "NOISE", "", "ARFCT", ...
        "", "STCH", "TCH", "SYSTOLE", "DIASTOLE", "NOTE", "MEASURE", "PWAVE", "BBB", ...
        "PACESP", "TWAVE", "RHYTHM", "UWAVE", "LEARN", "FLWAV", "VFON", "VFOFF", ...
        "AESC", "SVESC", "LINK", "NAPC", "PFUS", "WFON", "WFOFF", "RONT"
    ];

R_labels = strings(length(R_numbers),1); 
for i = 1:length(R_points)
    R_labels(i,1) = labelMap(R_numbers(i) + 1); 
end

% because we can't plot all the data in one plot, we just select part of it
time_period = [1790, 1805]; 
R_points_selected = R_points(R_points >= time_period(1) & R_points<=time_period(2));
R_numbers_selected = R_numbers(R_points >= time_period(1) & R_points<=time_period(2));
R_labels_selected = R_labels(R_points >= time_period(1) & R_points<=time_period(2));

t_selected = t(time_period(1)*fs : (time_period(2))*fs);
samples_selected = time_period(1)*fs : (time_period(2))*fs;


subplot(2,1,1)
plot(t_selected, signal(samples_selected,1));

for i = 1:length(R_points_selected)
    % Check if the current label is in the anomaly list
    if any(strcmp(R_labels_selected(i), anomalyLabels))
        text(R_points_selected(i), 1, R_labels_selected(i), 'FontSize', 7, 'Color', 'red', 'HorizontalAlignment', 'center');
    else
        text(R_points_selected(i), 1, R_labels_selected(i), 'FontSize', 7, 'HorizontalAlignment', 'center');
    end
end
xlabel("Time(s)",'Interpreter','latex')
ylabel("Amplitude(mV)",'Interpreter','latex')
ylim([-1, 1.5])
title('Channel One Signal','Interpreter','latex');

subplot(2,1,2)
plot(t_selected, signal(samples_selected,2));

for i = 1:length(R_points_selected)
    % Check if the current label is in the anomaly list
    if any(strcmp(R_labels_selected(i), anomalyLabels))
        text(R_points_selected(i), -1.2, R_labels_selected(i), 'FontSize', 7, 'Color', 'red', 'HorizontalAlignment', 'center');
    else
        text(R_points_selected(i), -1.2, R_labels_selected(i), 'FontSize', 7, 'HorizontalAlignment', 'center');
    end
end
xlabel("Time(s)",'Interpreter','latex')
ylabel("Amplitude(mV)",'Interpreter','latex')
ylim([-1.5, 1])
title('Channel Two Signal ','Interpreter','latex');


%% Q3
[existing_numbers, ~, idx] = unique(R_numbers);
for i = 1:length(existing_numbers)
    label_idx = find(idx == i);
    if ~isempty(label_idx)
        selected_idx = label_idx(randi(length(label_idx)));
        window_start = round(max((R_points(selected_idx) - 2) * fs, 1));
        window_end = round(min((R_points(selected_idx) + 2) * fs, length(t)));
        window_data = signal(window_start:window_end, 1);
        window_time = t(window_start:window_end);
        subplot(length(existing_numbers)/2, 2, i);
        plot(window_time, window_data, 'LineWidth', 1);
        hold on; 
        plot(R_points(selected_idx), 0, 'o', 'MarkerSize', 20);
        title(labelMap(existing_numbers(i) + 1), 'Interpreter', 'latex');
        xlabel("Time (s)", 'Interpreter', 'latex');
        ylabel("Amplitude (mV)", 'Interpreter', 'latex');
        xlim([t(window_start), t(window_end)]);
    end
end

%% Q4
clc;
abnormal_beats = find(R_numbers ~= 1);
abnormal_index = [];
for i = 600:(length(R_numbers)-360)
    if all(ismember(i + [-1, 0, 1], abnormal_beats))
        abnormal_index = i;
        break;
    end
end

abnormal_time = [R_points(abnormal_index - 1) - 0.5, R_points(abnormal_index + 1) + 0.5];
normal_time = [53.4, 55.4];

t_abnormal = round(abnormal_time * fs);
signal_abnormal = signal(t_abnormal(1):t_abnormal(2), 1);

t_normal = round(normal_time * fs);
signal_normal = signal(t_normal(1):t_normal(2), 1);

subplot(3, 2, 1);
plot(t(t_abnormal(1):t_abnormal(2)), signal_abnormal, 'LineWidth', 1);
xlabel("Time (s)", 'Interpreter', 'latex');
ylabel("Amplitude (mV)", 'Interpreter', 'latex');
title('3 Abnormal Heart Beats', 'Interpreter', 'latex');

subplot(3, 2, 2);
plot(t(t_normal(1):t_normal(2)), signal_normal, 'LineWidth', 1);
xlabel("Time (s)", 'Interpreter', 'latex');
ylabel("Amplitude (mV)", 'Interpreter', 'latex');
title('3 Normal Heart Beats', 'Interpreter', 'latex');

L_abnormal = length(signal_abnormal);
fft_signal_abnormal = fft(signal_abnormal);
fft_shifted_signal_abnormal = fft_signal_abnormal / L_abnormal; 
f_abnormal = (0:(L_abnormal/2)) * (fs/L_abnormal);

subplot(3, 2, 3);
plot(f_abnormal, 2 * abs(fft_shifted_signal_abnormal(1:L_abnormal/2 + 1))); 
title('Fourier Transform of Abnormal Signal', 'Interpreter', 'latex');
ylabel('Magnitude', 'Interpreter', 'latex');
xlabel('Frequency (Hz)', 'Interpreter', 'latex');
xlim([0 max(f_abnormal)]); 
grid on;
L_normal = length(signal_normal);
fft_signal_normal = fft(signal_normal);
fft_shifted_signal_normal = fft_signal_normal / L_normal; 
f_normal = (0:(L_normal/2)) * (fs/L_normal);

subplot(3, 2, 4);
plot(f_normal, 2 * abs(fft_shifted_signal_normal(1:L_normal/2 + 1))); 
title('Fourier Transform of Normal Signal', 'Interpreter', 'latex');
ylabel('Magnitude', 'Interpreter', 'latex');
xlabel('Frequency (Hz)', 'Interpreter', 'latex');
xlim([0 max(f_normal)]); 
grid on;

subplot(3, 2, 5);
spectrogram(signal_abnormal, hamming(128), 64, 128, fs, 'yaxis');
title('Spectrogram of Abnormal Signal', 'Interpreter', 'latex');
subplot(3, 2, 6);
spectrogram(signal_normal, hamming(128), 64, 128, fs, 'yaxis');
title('Spectrogram of Normal Signal', 'Interpreter', 'latex');

