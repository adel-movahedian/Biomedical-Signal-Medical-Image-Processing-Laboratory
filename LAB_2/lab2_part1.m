%% Q1
clc; clear;
load('X_org.mat')
load('X_noise.mat')
fs = 250;
plotEEG(X_org);
title('Original signal','Fontsize',14,'Interpreter','latex');























%% Q2
clc;
load('X_noise.mat')
fs = 250;
plotEEG(X_noise)
title('Noise signal','Fontsize',14,'Interpreter','latex')

%% Q3
clc;

P_X = sum(sum(X_org.^2));
P_N = sum(sum(X_noise.^2));

SNR_5 = -5;
SNR_15 = -15;

sigma_1 = sqrt((P_X/P_N)*10^(-1*(SNR_5)/10));
sigma_2 = sqrt((P_X/P_N)*10^(-1*(SNR_15)/10));

X_SNR_5 = X_org + sigma_1.*(X_noise);
X_SNR_15 = X_org + sigma_2.*(X_noise);

plotEEG(X_SNR_5);
title('-5 db noisy signal','Fontsize',14,'Interpreter','latex');


plotEEG(X_SNR_15);
title('-15 db noisy signal','Fontsize',14,'Interpreter','latex');

%% Q4
clc;
[F_5,W_5,~] = COM2R(X_SNR_5,32);
ICA_5_components = W_5*X_SNR_5;

[F_15,W_15,~]=COM2R(X_SNR_15,32);
ICA_15_components = W_15*X_SNR_15;

plotEEG(ICA_5_components);
title('Sources for SNR = -5dB','Fontsize',14,'Interpreter','latex');
plotEEG(ICA_15_components);
title('Sources for SNR = -15dB','Fontsize',14,'Interpreter','latex');

%% Q5 
clc; 

desired_components_5 = [2,5,12,23];
desired_components_15  = [15,19,20,30];

%% Q6
clc;

X_denoised_5 = F_5(:,desired_components_5)*ICA_5_components(desired_components_5,:);
X_denoised_15 = F_15(:,desired_components_15)*ICA_15_components(desired_components_15,:);

plotEEG(X_denoised_5);
title('denoised signal with SNR = -5dB','Fontsize',14,'Interpreter','latex');

plotEEG(X_denoised_15);
title('denoised signal with SNR = -15dB','Fontsize',14,'Interpreter','latex');

%% Q7
clc;
figure()
subplot(5,2,1)
plot(X_org(13,:));
xlim([0 10240])
title('Original signal of channel 13','Fontsize',14,'Interpreter','latex');
subplot(5,2,3)
plot(X_SNR_5(13,:));
xlim([0 10240])
title('Noisy signal of channel 13 with SNR = -5dB','Fontsize',14,'Interpreter','latex');
subplot(5,2,5)
plot(X_denoised_5(13,:));
xlim([0 10240])
title('Denoised signal of channel 13 Using ICA','Fontsize',14,'Interpreter','latex');
subplot(5,2,7)
plot(X_SNR_15(13,:));
xlim([0 10240])
title('Noisy signal of channel 13 with SNR = -15dB','Fontsize',14,'Interpreter','latex');
subplot(5,2,9)
plot(X_denoised_15(13,:));
xlim([0 10240])
title('Denoised signal of channel 13 with ICA','Fontsize',14,'Interpreter','latex');

subplot(5,2,2)
plot(X_org(24,:));
xlim([0 10240])
title('Original signal of channel 24','Fontsize',14,'Interpreter','latex');
subplot(5,2,4)
plot(X_SNR_5(24,:));
xlim([0 10240])
title('Noisy signal of channel 24 with SNR = -5dB','Fontsize',14,'Interpreter','latex');
subplot(5,2,6)
plot(X_denoised_5(24,:));
xlim([0 10240])
title('Denoised signal of channel 24 with ICA','Fontsize',14,'Interpreter','latex');
subplot(5,2,8)
plot(X_SNR_15(24,:));
xlim([0 10240])
title('Noisy signal of channel 24 with SNR = -15dB','Fontsize',14,'Interpreter','latex');
subplot(5,2,10)
plot(X_denoised_15(24,:));
xlim([0 10240])
title('Denoised signal of channel 24 with ICA','Fontsize',14,'Interpreter','latex');


%% Q8
clc;

RRMSE_5 = RRMSE(X_org,X_denoised_5);
RRMSE_15 = RRMSE(X_org,X_denoised_15);

disp('RRMSE for snr = -5: ');
disp(RRMSE_5);
disp('RRMSE for snr = -15: ');
disp(RRMSE_15);

%% Functions
function plotEEG(X) 
load('Electrodes.mat') ;
offset = max(abs(X(:))) ;
feq = 250 ;
ElecName = Electrodes.labels ;
disp_eeg(X,offset,feq,ElecName);
end

function result = RRMSE(X_org,X_den)

temp_num = sum(sum((X_org-X_den).^2,2));
temp_den = sum(sum((X_org.^2),2));

result = sqrt(temp_num)/sqrt(temp_den);

end

function [F,W,K]=COM2R(Y,Pest)
disp('COM2')

[N,TT]=size(Y);T=max(N,TT);N=min(N,TT);
if TT==N, Y=Y';[N,T]=size(Y);end; 
[U,S,V]=svd(Y',0);tol=max(size(S))*norm(S)*eps;
s=diag(S);I=find(s<tol);

%--- modif de Laurent le 03/02/2009
r = min(Pest,N);
U=U(:,1:r);
S=S(1:r,1:r);
V=V(:,1:r);
%---

Z=U'*sqrt(T);L=V*S'/sqrt(T);F=L; %%%%%% on a Y=L*Z;
%%%%%% INITIAL CONTRAST
T=length(Z);contraste=0;
for i=1:r,
 gii=Z(i,:)*Z(i,:)'/T;Z2i=Z(i,:).^2;;giiii=Z2i*Z2i'/T;
 qiiii=giiii/gii/gii-3;contraste=contraste+qiiii*qiiii;
end;
%%%% STEPS 3 & 4 & 5: Unitary transform
S=Z;
if N==2,K=1;else,K=1+round(sqrt(N));end;  % K= max number of sweeps
Rot=eye(r);
for k=1:K,                           %%%%%% strating sweeps
Q=eye(r);
  for i=1:r-1,
  for j= i+1:r,
    S1ij=[S(i,:);S(j,:)];
    [Sij,qij] = tfuni4(S1ij);    %%%%%% processing a pair
    S(i,:)=Sij(1,:);S(j,:)=Sij(2,:);
    Qij=eye(r);Qij(i,i)=qij(1,1);Qij(i,j)=qij(1,2);
    Qij(j,i)=qij(2,1);Qij(j,j)=qij(2,2);
    Q=Qij*Q;
  end;
  end;
Rot=Rot*Q';
end;                                    %%%%%% end sweeps
F=F*Rot;
%%%%%% FINAL CONTRAST
S=Rot'*Z;
T=length(S);contraste=0;
for i=1:r,
 gii=S(i,:)*S(i,:)'/T;S2i=S(i,:).^2;;giiii=S2i*S2i'/T;
 qiiii=giiii/gii/gii-3;contraste=contraste+qiiii*qiiii;
end;
%%%% STEP 6: Norming columns
delta=diag(sqrt(sum(F.*conj(F))));
%%%% STEP 7: Sorting
[d,I]=sort(-diag(delta));E=eye(r);P=E(:,I)';delta=P*delta*P';F=F*P';
%%%% STEP 8: Norming
F=F*inv(delta);
%%%% STEP 9: Phase of columns
[y,I]=max(abs(F));
for i=1:r,Lambda(i)=conj(F(I(i),i));end;Lambda=Lambda./abs(Lambda);
F=F*diag(Lambda);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCUL DE LA MATRICE DE FILTRAGE
%---------------------------------
W = pinv(F);
end


