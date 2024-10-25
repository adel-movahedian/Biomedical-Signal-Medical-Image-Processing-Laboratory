%% 1
load('EOG_sig.mat');
t = (0:length(Sig)-1)/fs;

figure;

plot(t,Sig);
title('Eye signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
legend('Left eye','Right eye');

%% 2
figure;
subplot(2, 1, 1);
L = length(Sig(1,:));
Y = fft(Sig(1,:));
P2 = abs(Y/L);
P1 = P2(1:L/2+1); 
P1(2:end-1) = 2*P1(2:end-1); 
f = fs*(0:(L/2))/L;
plot(f, P1);
title('left eye - Frequency Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;

subplot(2, 1, 2);
L = length(Sig(2,:));
Y = fft(Sig(2,:));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;
plot(f, P1);
title('Myopathy EMG - Frequency Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;

figure;
subplot(2, 1, 1);
spectrogram(Sig(1,:), 256, 250, 256, fs, 'yaxis'); 
title('Left Eye EOG - Time-Frequency Representation');
colorbar;
xlabel('Time (s)');
ylabel('Frequency (Hz)');

subplot(2, 1, 2);
spectrogram(Sig(2,:), 256, 250, 256, fs, 'yaxis');
title('Right Eye EOG - Time-Frequency Representation');
colorbar;
xlabel('Time (s)');
ylabel('Frequency (Hz)');

