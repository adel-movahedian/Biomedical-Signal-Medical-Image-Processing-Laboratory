%% 1
load('EMG_sig.mat');
t_healthym = (0:length(emg_healthym)-1)/fs;
t_myopathym = (0:length(emg_myopathym)-1)/fs;
t_neuropathym = (0:length(emg_neuropathym)-1)/fs;

figure;
subplot(3,1,1);
plot(t_healthym, emg_healthym);
title('Healthy EMG Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(3,1,2);
plot(t_myopathym, emg_myopathym);
title('Myopathy EMG Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(3,1,3);
plot(t_neuropathym, emg_neuropathym);
title('Neuropathy EMG Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% 2
figure;
subplot(3, 1, 1);
L = length(emg_healthym);
Y = fft(emg_healthym);
P2 = abs(Y/L);
P1 = P2(1:L/2+1); 
P1(2:end-1) = 2*P1(2:end-1); 
f = fs*(0:(L/2))/L;
plot(f, P1);
title('Healthy EMG - Frequency Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;

subplot(3, 1, 2);
L = length(emg_myopathym);
Y = fft(emg_myopathym);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;
plot(f, P1);
title('Myopathy EMG - Frequency Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;

subplot(3, 1, 3);
L = length(emg_neuropathym);
Y = fft(emg_neuropathym);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;
plot(f, P1);
title('Neuropathy EMG - Frequency Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;

figure;
subplot(3, 1, 1);
spectrogram(emg_healthym, 256, 250, 256, fs, 'yaxis'); 
title('Healthy EMG - Time-Frequency Representation');
colorbar;
xlabel('Time (s)');
ylabel('Frequency (Hz)');

subplot(3, 1, 2);
spectrogram(emg_myopathym, 256, 250, 256, fs, 'yaxis');
title('Myopathy EMG - Time-Frequency Representation');
colorbar;
xlabel('Time (s)');
ylabel('Frequency (Hz)');

subplot(3, 1, 3);
spectrogram(emg_neuropathym, 256, 250, 256, fs, 'yaxis');
title('Neuropathy EMG - Time-Frequency Representation');
colorbar;
xlabel('Time (s)');
ylabel('Frequency (Hz)');


