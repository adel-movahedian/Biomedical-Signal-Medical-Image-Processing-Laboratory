%% part1
load('EEG_sig.mat');
channel_5_data = Z(5, :);  
fs = des.samplingfreq; 
time = (0:length(channel_5_data) - 1) / fs;
figure;
plot(time, channel_5_data, 'r', 'LineWidth', 1);
xlabel('Time (seconds)');
ylabel(['EEG Signal - Channel 5 (' des.channelnames{5} ')']);
title('EEG Signal of Channel 5 Over Time');
set(gcf, 'Position', [100, 100, 1000, 400]); 
%% part2
clc;
time_1=[1,15*fs];
time_2=[18*fs, 40*fs];
time_3=[45*fs, 50*fs];
time_4=[50*fs, length(channel_5_data)];
figure()
subplot(2,2,1)
plot(time(time_1(1):time_1(2)),channel_5_data(time_1(1):time_1(2)));
xlabel('t');
ylabel(channel_names(5));
title('Channel 5 EEG signal for t = 0 : 15s','Interpreter','latex')
grid on
subplot(2,2,2)
plot(time(time_2(1):time_2(2)),channel_5_data(time_2(1):time_2(2)));
xlabel('t');
ylabel(channel_names(5));
title('Channel 5 EEG signal for t = 18 : 40s','Interpreter','latex')
grid on
subplot(2,2,3)
plot(t(time_3(1):time_3(2)),channel_5_data(time_3(1):time_3(2)));
xlabel('t');
ylabel(channel_names(5));
title('Channel 5 EEG signal for t = 45 : 50s','Interpreter','latex')
grid on
subplot(2,2,4)
plot(time(time_4(1):time_4(2)),channel_5_data(time_4(1):time_4(2)));
xlabel('t');
ylabel(channel_names(5));
title('Channel 5 EEG signal for t = 50 : end','Interpreter','latex')
grid on

%% part3
channel_10_data = Z(10, :); 
figure;
subplot(2, 1, 1);
plot(time, channel_10_data, 'r');
title(['EEG Signal - Channel 10 (' des.channelnames{10} ')']);
xlabel('Time (seconds)');
ylabel('Amplitude');
subplot(2, 1, 2);
plot(time, channel_5_data, 'b');
title(['EEG Signal - Channel 5 (' des.channelnames{5} ')']);
xlabel('Time (seconds)');
ylabel('Amplitude');


%% part4
load('EEG_sig.mat');
Z = Z;  
fs = des.samplingfreq;
ElecName = des.channelnames;  
offset = max(max(abs(Z)))/5; 
disp_eeg(Z, offset, fs, ElecName);
%% part6
channel_C3 = Z(5, :);  
fs = des.samplingfreq;  
N = length(channel_C3);  
time = (0:N-1) / fs;
intervals = [2 7; 30 35; 42 47; 50 55];
for i = 1:size(intervals, 1)
    start_time = intervals(i, 1);
    end_time = intervals(i, 2);
    start_idx = round(start_time * fs);
    end_idx = round(end_time * fs);
    segment_time = time(start_idx:end_idx);
    segment_signal = channel_C3(start_idx:end_idx);
    figure;
    subplot(2, 1, 1);
    plot(segment_time, segment_signal);
    title(['EEG Signal in Time Domain - Interval ' num2str(start_time) ' to ' num2str(end_time) ' seconds']);
    xlabel('Time (seconds)');
    ylabel('Amplitude');
    L = length(segment_signal); 
    Y = fft(segment_signal);
    P2 = abs(Y / L);  
    P1 = P2(1:L/2+1); 
    P1(2:end-1) = 2 * P1(2:end-1); 
    f = fs * (0:(L/2)) / L;
    subplot(2, 1, 2);
    plot(f, P1);
    xlim([0 50]); 
    title(['EEG Signal in Frequency Domain - Interval ' num2str(start_time) ' to ' num2str(end_time) ' seconds']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
end

%% part7
channel_C3 = Z(5, :);  
fs = des.samplingfreq; 
N = length(channel_C3); 
time = (0:N-1) / fs;
intervals = [2 7; 30 35; 42 47; 50 55];
for i = 1:size(intervals, 1)
    start_time = intervals(i, 1);
    end_time = intervals(i, 2);
    start_idx = round(start_time * fs);
    end_idx = round(end_time * fs);
    segment_signal = channel_C3(start_idx:end_idx);
    figure;
    subplot(2, 1, 1);
    segment_time = time(start_idx:end_idx);
    plot(segment_time, segment_signal);
    title(['EEG Signal in Time Domain - Interval ' num2str(start_time) ' to ' num2str(end_time) ' seconds']);
    xlabel('Time (seconds)');
    ylabel('Amplitude');
    subplot(2, 1, 2);
    [pxx, f] = pwelch(segment_signal, [], [], [], fs); 
    plot(f, 10*log10(pxx));
    xlim([0 50]); 
    title(['Power Spectral Density (PSD) - Interval ' num2str(start_time) ' to ' num2str(end_time) ' seconds']);
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
end

%%
%q1_p8
% Load the .mat file
load('EEG_sig.mat');

% Extract the 5th channel data (C3 is the 5th channel)
channel_C3 = Z(5, :);  

% Extract the sampling frequency from 'des'
fs = des.samplingfreq;  % Sampling frequency
N = length(channel_C3);  % Length of the signal

% Define time intervals in seconds
intervals = [2 7; 30 35; 42 47; 50 55];

% Loop over the intervals to compute and plot the spectrogram
for i = 1:size(intervals, 1)
    % Define the start and end times in samples
    start_time = intervals(i, 1);
    end_time = intervals(i, 2);
    
    % Convert time to sample indices
    start_idx = round(start_time * fs);
    end_idx = round(end_time * fs);
    
    % Extract the segment of the signal
    segment_signal = channel_C3(start_idx:end_idx);
    
    % Create the spectrogram
    figure;
    window_length = 128;  % Length of the window
    overlap_length = 64;   % Length of overlap
    nfft = 128;            % Number of points for FFT
    spectrogram(segment_signal, hamming(window_length), overlap_length, nfft, fs);
    
    % Plot the spectrogram
    %imagesc(T, F, 10*log10(abs(S)));  % Convert to dB
    axis xy;
    colorbar;
    xlabel('Frequency (Hz)');
    Time (s)
    ylabel('Time (s)');
    title(['Spectrogram of EEG Signal - Interval ' num2str(start_time) ' to ' num2str(end_time) ' seconds']);

end
%% part8
load('EEG_sig.mat');
channel_C3 = Z(5, :);  
fs = des.samplingfreq; 
N = length(channel_C3); 

intervals = [2 7; 30 35; 42 47; 50 55];
windowLength = 128; 
overlapLength = 64; 
nfft = 128;    
for i = 1:size(intervals, 1)
    start_time = intervals(i, 1);
    end_time = intervals(i, 2);
    start_idx = round(start_time * fs);
    end_idx = round(end_time * fs);
    segment_signal = channel_C3(start_idx:end_idx);
    figure;
    spectrogram(segment_signal, hamming(windowLength), overlapLength, nfft, fs, 'yaxis');
    title(['Spectrogram of EEG Signal - Interval ' num2str(start_time) ' to ' num2str(end_time) ' seconds']);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    colorbar; 
end

%% part9
channel_C3 = Z(5, :);  

% Extract the sampling frequency from 'des'
fs = des.samplingfreq;  % Original sampling frequency

% Define the time intervals for the second segment
start_time = 30;  % start at 30 seconds
end_time = 35;    % end at 35 seconds

% Convert time to sample indices
start_idx = round(start_time * fs);
end_idx = round(end_time * fs);

% Extract the segment of the signal
segment_signal = channel_C3(start_idx:end_idx);

% Design a low-pass filter
new_fs = fs / 2;  % New sampling frequency
cutoff_freq = new_fs / 2; % Cutoff frequency for the low-pass filter
order = 4; % Order of the filter
lpFilt = designfilt('lowpassiir', 'FilterOrder', order, ...
                     'HalfPowerFrequency', cutoff_freq, ...
                     'DesignMethod', 'butter', ...
                     'SampleRate', fs);

% Apply the low-pass filter
filtered_signal = filtfilt(lpFilt, segment_signal);

% Downsample the signal
downsample_factor = 2; % Downsample by a factor of 2
downsampled_signal = downsample(filtered_signal, downsample_factor);

% Calculate the new sampling frequency
fs_new = fs / downsample_factor;

% Plot the original segment signal
figure;
subplot(3, 1, 1);
t_original = (0:length(segment_signal)-1) / fs;
plot(t_original, segment_signal);
title('Original EEG Signal (Segment)');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot the filtered and downsampled signal
t_downsampled = (0:length(downsampled_signal)-1) / fs_new;
subplot(3, 1, 2);
plot(t_downsampled, downsampled_signal);
title('Filtered and Downsampled EEG Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Compute DFT of the original and downsampled signal
N_original = length(segment_signal);
N_downsampled = length(downsampled_signal);

% Compute DFT
dft_original = abs(fft(segment_signal, N_original));
dft_downsampled = abs(fft(downsampled_signal, N_downsampled));

% Frequency axis
f_original = (0:N_original-1) * fs / N_original;
f_downsampled = (0:N_downsampled-1) * fs_new / N_downsampled;

% Plot DFT
subplot(3, 1, 3);
hold on;
plot(f_original(1:N_original/2), dft_original(1:N_original/2), 'b', 'DisplayName', 'Original DFT');
plot(f_downsampled(1:N_downsampled/2), dft_downsampled(1:N_downsampled/2), 'r', 'DisplayName', 'Downsampled DFT');
title('DFT of Original and Downsampled Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend;
hold off;

% Compute STFT for original and downsampled signals
windowLength = 128; % Length of the window
overlapLength = 64;  % Length of overlap
nfft = 128;         % Number of FFT points

% Original STFT
figure;
subplot(2, 1, 1);
spectrogram(segment_signal, hamming(windowLength), overlapLength, nfft, fs, 'yaxis');
title('STFT of Original Signal');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

% Downsampled STFT
subplot(2, 1, 2);
spectrogram(downsampled_signal, hamming(windowLength), overlapLength, nfft, fs_new, 'yaxis');
title('STFT of Downsampled Signal');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

