%% Part 1: Q1
clc;

ECG_mother = load('mecg1.dat');
ECG_fetus = load('fecg1.dat');
Noise = load('noise1.dat');

Mixed_signal = ECG_mother + ECG_fetus + Noise ;
fs = 256;
L_signal = length(ECG_fetus)/fs;
t = 0:1/fs:L_signal-1/fs;

subplot(4,1,1);
plot(t,ECG_mother);
title('ECG for mother','Fontsize',14,'Interpreter','latex');
grid on;
xlabel('time(s)', 'Interpreter','latex');
ylabel('Amplitude(mV)', 'Interpreter','latex');

subplot(4,1,2);
plot(t,ECG_fetus);
title('ECG for fetus','Fontsize',14,'Interpreter','latex');
grid on;
xlabel('time(s)', 'Interpreter','latex');
ylabel('Amplitude(mV)', 'Interpreter','latex');

subplot(4,1,3);
plot(t,Noise);
title('noise','Fontsize',14,'Interpreter','latex');
grid on;
xlabel('time(s)', 'Interpreter','latex');
ylabel('Amplitude(mV)', 'Interpreter','latex');

subplot(4,1,4) ;
plot(t,Mixed_signal);
title('mixed signal','Fontsize',14,'Interpreter','latex');
grid on;
xlabel('time(s)', 'Interpreter','latex');
ylabel('Amplitude(mV)', 'Interpreter','latex');

%% Q2
clc;

subplot(3,1,1) ;
[pxx,f] = pwelch(ECG_mother,[],[],[],fs);
plot(f,pxx,'Linewidth',1);
title('power spectrum of the ECG for mother','FontSize',14,'Interpreter','latex');
grid on;
xlabel('Frequency (Hz)', 'Interpreter','latex');
ylabel('Power', 'Interpreter','latex');

subplot(3,1,2) ;
[pxx,f] = pwelch(ECG_fetus,[],[],[],fs);
plot(f,pxx,'Linewidth',1);
title('power spectrum of the ECG for fetus','FontSize',14,'Interpreter','latex');
grid on;
xlabel('Frequency (Hz)', 'Interpreter','latex');
ylabel('Power', 'Interpreter','latex');

subplot(3,1,3) ;
[pxx,f] = pwelch(Noise,[],[],[],fs);
plot(f,pxx,'Linewidth',1);
title('power spectrum of noise','FontSize',14,'Interpreter','latex');
grid on;
xlabel('Frequency (Hz)', 'Interpreter','latex');
ylabel('Power', 'Interpreter','latex');

%% Q3
clc;

mean_ECG_mother = mean(ECG_mother);
variance_ECG_mother = var(ECG_mother);
output = ['mean ECG mother = ',num2str(mean_ECG_mother),' --------- varince ECG mother = ',num2str(variance_ECG_mother)];
disp(output);

mean_ECG_fetus = mean(ECG_fetus);
variance_ECG_fetus = var(ECG_fetus);
output = ['mean ECG fetus = ',num2str(mean_ECG_fetus),' --------- varince ECG fetus = ',num2str(variance_ECG_fetus)];
disp(output);

mean_Noise = mean(Noise);
variance_Noise = var(Noise);
output = ['mean Noise = ',num2str(mean_Noise),' --------- varince Noise = ',num2str(variance_Noise)];
disp(output);

%% Q4
clc;

subplot(3,1,1) ;
histogram(ECG_mother,150)
title('histogram plot of the ECG for mother','Fontsize',14,'Interpreter','latex');
xlim([-3 3]);
grid on;

subplot(3,1,2) ;
histogram(ECG_fetus,150)
title('histogram plot of the ECG for fetus','Fontsize',14,'Interpreter','latex');
xlim([-3 3]);
grid on;

subplot(3,1,3) ;
histogram(Noise,150)
title('histogram plot of the noise','Fontsize',14,'Interpreter','latex');
xlim([-3 3]);
grid on;

kurtosis_ECG_mother = kurtosis(ECG_mother);
disp(['kurtosis of ECG mother = ',num2str(kurtosis_ECG_mother)]);
kurtosis_ECG_fetus = kurtosis(ECG_fetus);
disp(['kurtosis of ECG fetus = ',num2str(kurtosis_ECG_fetus)]);
kurtosis_Noise = kurtosis(Noise);
disp(['kurtosis of noise = ',num2str(kurtosis_Noise)]);


%% Part 2: Q1
clc;

X = load('X.dat');
fs = 256;
plot3ch(X,fs,'Main Signal');

[U,S,V] = svd(X);

%% Q2
clc;

plot3dv(V(:,1),S(1,1),'r')
plot3dv(V(:,2),S(2,2),'g')
plot3dv(V(:,3),S(3,3),'b')

%% Q3
clc;

plot3ch(U(:,1:3),fs,'Seperated signals after SVD');

figure;
stem([S(1,1) S(2,2) S(3,3)],'LineWidth',2);
xlim([0 4]);
title('eigenspectrum','Fontsize',14,'Interpreter','latex');
grid on;

%% Q4
clc;

S_selected = zeros(2560,3);
S_selected (2,2) = S(2,2);

Y_reconstructed = U * S_selected * transpose(V);

plot3ch(Y_reconstructed(:,1:3),fs,'Reconstructed ECG signal of fetus');

%% Part 3: Q1
clc; 

X = load('X.dat');
fs=256;
[W,Z_hat_T] = ica(X');
A = inv(W);

%% Q2
clc;

plot3ch(X);

figure;
plot3dv(A(:,1),[],'r')
plot3dv(A(:,2),[],'g')
plot3dv(A(:,3),[],'b')
title('W inverse vectors','Fontsize',14,'Interpreter','latex');

%% Q3
clc;

Z_hat = transpose(Z_hat_T);
plot3ch(Z_hat,fs,'Z hat columns');
Z_selected = zeros(2560,3);
Z_selected(:,3) = Z_hat(:,3); 

X_reconstructed = transpose(A * transpose(Z_selected));

%% Q4
clc;

plot3ch(X_reconstructed,fs);
title('Recuntructed ECG using ICA');

%% Part 4: Q1
clc; close all;

scatter3(Y_reconstructed(:,1), Y_reconstructed(:,2), Y_reconstructed(:,3),'filled','MarkerFaceColor',[0.9290 0.6940 0.1250]);
hold on;
scatter3(X(:,1), X(:,2), X(:,3),'filled','MarkerFaceColor',[0 0.4470 0.7410]);
hold on;
scatter3(X_reconstructed(:,1), X_reconstructed(:,2), X_reconstructed(:,3),'filled','MarkerFaceColor',[0.8500 0.3250 0.0980]);

hold on;
plot3dv(V(:,1).',S(:,1),'r')
plot3dv(V(:,2).',S(:,2),'r')
plot3dv(V(:,3).',S(:,3),'r')

hold on;
plot3dv(A(:,1),[],'black')
plot3dv(A(:,2),[],'black')
plot3dv(A(:,3),[],'black')

legend('Using SVD','X','Using ICA','V directions','A directions');

figure;

plot3dv(V(:,1).',S(:,1),'r')
plot3dv(V(:,2).',S(:,2),'r')
plot3dv(V(:,3).',S(:,3),'r')

hold on;
plot3dv(A(:,1),[],'black')
plot3dv(A(:,2),[],'black')
plot3dv(A(:,3),[],'black')

legend('V');

deg1 = acosd(dot(V(:, 1), V(:, 2)) / (sum(abs(V(:, 1))) * sum(abs(V(:, 2)))));
deg2 = acosd(dot(V(:, 2), V(:, 3)) / (sum(abs(V(:, 2))) * sum(abs(V(:, 3)))));
deg3 = acosd(dot(V(:, 1), V(:, 3)) / (sum(abs(V(:, 1))) * sum(abs(V(:, 3)))));
disp('angles between V = ');
disp([deg1,deg2,deg3]);

norm1 = norm(V(:,1));
norm2 = norm(V(:,2));
norm3 = norm(V(:,3));
disp('norms V = ');
disp([norm1,norm2,norm3]);

deg1 = acosd(dot(A(:, 1), A(:, 2)) / (sum(abs(A(:, 1))) * sum(abs(A(:, 2)))));
deg2 = acosd(dot(A(:, 2), A(:, 3)) / (sum(abs(A(:, 2))) * sum(abs(A(:, 3)))));
deg3 = acosd(dot(A(:, 1), A(:, 3)) / (sum(abs(A(:, 1))) * sum(abs(A(:, 3)))));
disp('angles between A = ');
disp([deg1,deg2,deg3]);

norm1 = norm(A(:,1));
norm2 = norm(A(:,2));
norm3 = norm(A(:,3));
disp('norms A  = ');
disp([norm1,norm2,norm3]);


%% Q2
clc; close all;

fecg2 = load('fecg2.dat');

t = 0:1/fs:10 - 1/fs;

Y_reconstructed_normal = Y_reconstructed(:,1)/norm(Y_reconstructed(:,1));
X_reconstructed_normal = X_reconstructed(:,1)/norm(X_reconstructed(:,1)); 

figure;

subplot(3,1,1);
plot(t,fecg2/norm(fecg2));
title('Normalized Fetal ECG','Fontsize',14,'Interpreter','latex');
xlabel('time(s)','Interpreter','latex');
ylabel('Normalized','Interpreter','latex');

subplot(3,1,2);
plot(t,Y_reconstructed_normal);
title('Using SVD','Fontsize',14,'Interpreter','latex');
xlabel('time(s)','Interpreter','latex');
ylabel('Normalized','Interpreter','latex');

subplot(3,1,3);
plot(t,X_reconstructed_normal);
title('Using ICA','Fontsize',14,'Interpreter','latex');
xlabel('time(s)','Interpreter','latex');
ylabel('Normalized','Interpreter','latex');

%% Q3
corr_coef_SVD = corrcoef(fecg2,Y_reconstructed_normal);
disp(['correlation coef Using SVD = ',num2str(corr_coef_SVD(1,2))]);

corr_coef_ICA = corrcoef(fecg2,X_reconstructed_normal);
disp(['correlation coef Using ICA = ',num2str(corr_coef_ICA(1,2))]);
