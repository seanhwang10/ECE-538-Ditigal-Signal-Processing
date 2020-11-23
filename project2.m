%-----------------------------%
%         Sean Hwang          %
% ECE 538 - MATLAB Project 2  %
%     Source Code File        %
%-----------------------------%

clc
clf 
clear all
close all 

% A) h_0 = {1/sqrt(2), 1/sqrt(2)} ,, h_1 = {1/sqrt(2), -1/sqrt(2)} 
x=randn(1,128);
M=8;
h0=[1 1]; h1=[1 -1];
h00=[1 0 1];   h10=[1 0 -1];
h000=[1 0 0 0 1]; h100=[1 0 0 0 -1];
H_tmp=conv(h0,h00); H(1,:)=conv(H_tmp,h000); G(1,:)=H(1,:);
H_tmp=conv(h0,h00); H(2,:)=conv(H_tmp,h100); G(2,:)=-H(2,:);
H_tmp=conv(h0,h10); H(3,:)=conv(H_tmp,h000); G(3,:)=-H(3,:);
H_tmp=conv(h0,h10); H(4,:)=conv(H_tmp,h100); G(4,:)=H(4,:);
H_tmp=conv(h1,h00); H(5,:)=conv(H_tmp,h000); G(5,:)=-H(5,:);
H_tmp=conv(h1,h00); H(6,:)=conv(H_tmp,h100); G(6,:)=H(6,:);
H_tmp=conv(h1,h10); H(7,:)=conv(H_tmp,h000); G(7,:)=H(7,:);
H_tmp=conv(h1,h10); H(8,:)=conv(H_tmp,h100); G(8,:)=-H(8,:);

% i)   Please see attached for figure 1a)
for ind=1:M
hf(ind,:)=abs(fftshift(fft(H(ind,:),512)));
end
domega=2*pi/512;
omega=-pi:domega:pi-domega;
figure(1);
plot(omega,hf(1,:),omega,hf(2,:), ...
     omega,hf(3,:),omega,hf(4,:), ...
     omega,hf(5,:),omega,hf(6,:), ...
     omega,hf(7,:),omega,hf(8,:));
axis([-pi pi 0 max(hf(1,:))])
legend('H0(w)','H1(w)','H2(w)','H3(w)','H4(w)','H5(w)','H6(w)','H7(w)');
title('fig 1a: Frequency Response of h0[n] through h7[n]');
xlabel('\omega (radians/s)');
% ii)  Please see attached for figure 1b)
for ind=1:M
gf(ind,:)=abs(fftshift(fft(G(ind,:),512)));
end
domega=2*pi/512;
omega=-pi:domega:pi-domega;
figure(2);
plot(omega,gf(1,:),omega,gf(2,:), ...
     omega,gf(3,:),omega,gf(4,:), ...
     omega,gf(5,:),omega,gf(6,:), ...
     omega,gf(7,:),omega,gf(8,:));
axis([-pi pi 0 max(gf(1,:))])
legend('G0(w)','G1(w)','G2(w)','G3(w)','G4(w)','G5(w)','G6(w)','G7(w)');
title('fig 1b: Frequency Response of g0[n] through g7[n]');
xlabel('\omega (radians/s)');

% iii)Table 1
% H*H' =
%   8     0     0     0     0     0     0     0
%   0     8     0     0     0     0     0     0
%   0     0     8     0     0     0     0     0
%   0     0     0     8     0     0     0     0 
%   0     0     0     0     8     0     0     0
%   0     0     0     0     0     8     0     0
%   0     0     0     0     0     0     8     0
%   0     0     0     8     0     0     0     8 
% iv) Please see attached figure 1c)
% v)  Please see attached figure 1d)
for m=1:M;
W(m,:)=conv(x,H(m,:));
X(m,:)=W(m,1:M:length(W(m,:)));
end
for m=1:M;
Z(m,:)=zeros(1,M*length(X(m,:)));
Z(m,1:M:length(Z(m,:)))=X(m,:);
Y(m,:)=conv(Z(m,:),G(m,:));
end
y=zeros(1,length(Y(1,:)));
for m=1:M;
y=y+Y(m,:);
end
domega=2*pi/1024;
omega=-pi:domega:pi-domega;
yf1=abs(fftshift(fft(x,1024)));
yf2=abs(fftshift(fft(y,1024)));
figure(3);
plot(omega,yf1,'Linewidth',1)
axis([-pi pi 0 max(yf1)])
title('fig 1c: DTFT of Gaussian Random Process (mean = 0 with unit power)');
xlabel('\omega (radians/s)');
figure(4);
plot(omega,yf2,'Linewidth',1)
axis([-pi pi 0 max(yf2)])
title('fig 1d: DTFT of output of Gaussian Random Process');
xlabel('\omega (radians/s)');
 
% B)
N=16; beta=0.35;
n=-N:(N-1); %이렇게 해야 n 이 0 ... 31 으로 됨 
n=n+0.5;
h=2*beta*cos((1+beta)*pi*n/2)./(pi*(1-4*beta^2*n.^2));
h=h+sin((1-beta)*pi*n/2)./(pi*(n-4*beta^2*n.^3));
M=8;
h0=h;
h1=(-1).^(0:(length(n)-1)).*h; 
h00=zeros(1,2*length(h)); h10=h00;
h00(1,1:2:length(h00))=h0;
h10(1,1:2:length(h10))=h1;
h000=zeros(1,4*length(h)); h100=h000;
h000(1,1:4:length(h000))=h0;
h100(1,1:4:length(h100))=h1;
H_tmp=conv(h0,h00); H_b(1,:)=conv(H_tmp,h000); G_b(1,:)=H_b(1,:);
H_tmp=conv(h0,h00); H_b(2,:)=conv(H_tmp,h100); G_b(2,:)=-H_b(2,:);
H_tmp=conv(h0,h10); H_b(3,:)=conv(H_tmp,h000); G_b(3,:)=-H_b(3,:);
H_tmp=conv(h0,h10); H_b(4,:)=conv(H_tmp,h100); G_b(4,:)=H_b(4,:);
H_tmp=conv(h1,h00); H_b(5,:)=conv(H_tmp,h000); G_b(5,:)=-H_b(5,:);
H_tmp=conv(h1,h00); H_b(6,:)=conv(H_tmp,h100); G_b(6,:)=H_b(6,:);
H_tmp=conv(h1,h10); H_b(7,:)=conv(H_tmp,h000); G_b(7,:)=H_b(7,:);
H_tmp=conv(h1,h10); H_b(8,:)=conv(H_tmp,h100); G_b(8,:)=-H_b(8,:);
% i)   Please see attached for figure 2a)
for ind=1:M
hf(ind,:)=abs(fftshift(fft(H_b(ind,:),512)));
end
domega=2*pi/512;
omega=-pi:domega:pi-domega;
figure(5);
plot(omega,hf(1,:),omega,hf(2,:), ...
     omega,hf(3,:),omega,hf(4,:), ...
     omega,hf(5,:),omega,hf(6,:), ...
     omega,hf(7,:),omega,hf(8,:));
axis([-pi pi 0 max(hf(1,:))])
legend('H0(w)','H1(w)','H2(w)','H3(w)','H4(w)','H5(w)','H6(w)','H7(w)');
title('fig 2a: Frequency Response of h0[n] through h7[n]');
xlabel('\omega (radians/s)');
% ii)  Please see attached for figure 2b)
for ind=1:M
gf(ind,:)=abs(fftshift(fft(G_b(ind,:),512)));
end
domega=2*pi/512;
omega=-pi:domega:pi-domega;
figure(6);
plot(omega,gf(1,:),omega,gf(2,:), ...
     omega,gf(3,:),omega,gf(4,:), ...
     omega,gf(5,:),omega,gf(6,:), ...
     omega,gf(7,:),omega,gf(8,:));
axis([-pi pi 0 max(gf(1,:))])
legend('G0(w)','G1(w)','G2(w)','G3(w)','G4(w)','G5(w)','G6(w)','G7(w)');
title('fig 2b: Frequency Response of g0[n] through g7[n]');
xlabel('\omega (radians/s)');
% iii) Table 2
% H_b*H_b' =
%    0.1250    0.0000    0.0000    0.0000    0.0000    0.0000   -0.0000   -0.0000
%    0.0000    0.1250   -0.0000    0.0000    0.0000    0.0000    0.0000    0.0000
%    0.0000   -0.0000    0.1250    0.0000   -0.0000    0.0000    0.0000    0.0000
%    0.0000    0.0000    0.0000    0.1249   -0.0000    0.0000    0.0000    0.0000
%    0.0000    0.0000   -0.0000   -0.0000    0.1250    0.0000    0.0000    0.0000
%    0.0000    0.0000    0.0000    0.0000    0.0000    0.1250   -0.0000    0.0000
%   -0.0000    0.0000    0.0000    0.0000    0.0000   -0.0000    0.1250    0.0000
%   -0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.1249
%                                       
% iv) Please see attached figure 2c)
% v)  Please see attached figure 2d)
for m=1:M;
W_b(m,:)=conv(x,H_b(m,:));
X_b(m,:)=W_b(m,1:M:length(W_b(m,:)));
end
for m=1:M;
Z_b(m,:)=zeros(1,M*length(X_b(m,:)));
Z_b(m,1:M:length(Z_b(m,:)))=X_b(m,:);
Y_b(m,:)=conv(Z_b(m,:),G_b(m,:));
end
y_b=zeros(1,length(Y_b(1,:)));
for m=1:M;
y_b=y_b+Y_b(m,:);
end
domega=2*pi/1024;
omega=-pi:domega:pi-domega;
yf1=abs(fftshift(fft(x,1024)));
yf2=M*abs(fftshift(fft(y_b,1024)));
figure(7);
plot(omega,yf1,'Linewidth',1)
axis([-pi pi 0 max(yf1)])
title('fig 2c: DTFT of Gaussian Random Process (mean = 0 with unit power)');
xlabel('\omega (radians/s)');
figure(8);
plot(omega,yf2,'Linewidth',1)
axis([-pi pi 0 max(yf2)])
title('fig 2d: DTFT of output of Gaussian Random Process');
xlabel('\omega (radians/s)');
 
% C)
N=24; beta=0.1;
n=-N:(N-1);
n=n+0.5;
h=2*beta*cos((1+beta)*pi*n/2)./(pi*(1-4*beta^2*n.^2));
h=h+sin((1-beta)*pi*n/2)./(pi*(n-4*beta^2*n.^3));
M=8;
h0=h;
h1=(-1).^(0:(length(n)-1)).*h; 
h00=zeros(1,2*length(h)); h10=h00;
h00(1,1:2:length(h00))=h0;
h10(1,1:2:length(h10))=h1;
h000=zeros(1,4*length(h)); h100=h000;
h000(1,1:4:length(h000))=h0;
h100(1,1:4:length(h100))=h1;
H_tmp=conv(h0,h00); H_c(1,:)=conv(H_tmp,h000); G_c(1,:)=H_c(1,:);
H_tmp=conv(h0,h00); H_c(2,:)=conv(H_tmp,h100); G_c(2,:)=-H_c(2,:);
H_tmp=conv(h0,h10); H_c(3,:)=conv(H_tmp,h000); G_c(3,:)=-H_c(3,:);
H_tmp=conv(h0,h10); H_c(4,:)=conv(H_tmp,h100); G_c(4,:)=H_c(4,:);
H_tmp=conv(h1,h00); H_c(5,:)=conv(H_tmp,h000); G_c(5,:)=-H_c(5,:);
H_tmp=conv(h1,h00); H_c(6,:)=conv(H_tmp,h100); G_c(6,:)=H_c(6,:);
H_tmp=conv(h1,h10); H_c(7,:)=conv(H_tmp,h000); G_c(7,:)=H_c(7,:);
H_tmp=conv(h1,h10); H_c(8,:)=conv(H_tmp,h100); G_c(8,:)=-H_c(8,:);
% i)   Please see attached for figure 3a)
for ind=1:M
hf(ind,:)=abs(fftshift(fft(H_c(ind,:),512)));
end
domega=2*pi/512;
omega=-pi:domega:pi-domega;
figure(9);
plot(omega,hf(1,:),omega,hf(2,:), ...
     omega,hf(3,:),omega,hf(4,:), ...
     omega,hf(5,:),omega,hf(6,:), ...
     omega,hf(7,:),omega,hf(8,:));
axis([-pi pi 0 max(hf(1,:))])
legend('H0(w)','H1(w)','H2(w)','H3(w)','H4(w)','H5(w)','H6(w)','H7(w)');
title('fig 3a: Frequency Response of h0[n] through h7[n]');
xlabel('\omega (radians/s)');
% ii)  Please see attached for figure 3b)
for ind=1:M
gf(ind,:)=abs(fftshift(fft(G_c(ind,:),512)));
end
domega=2*pi/512;
omega=-pi:domega:pi-domega;
figure(10);
plot(omega,gf(1,:),omega,gf(2,:), ...
     omega,gf(3,:),omega,gf(4,:), ...
     omega,gf(5,:),omega,gf(6,:), ...
     omega,gf(7,:),omega,gf(8,:));
axis([-pi pi 0 max(gf(1,:))])
legend('G0(w)','G1(w)','G2(w)','G3(w)','G4(w)','G5(w)','G6(w)','G7(w)');
title('fig 3b: Frequency Response of g0[n] through g7[n]');
xlabel('\omega (radians/s)');
% iii) Table 3
% H_c*H_c' =
%    0.1250   -0.0000   -0.0000    0.0000    0.0000    0.0000   -0.0000   -0.0000
%   -0.0000    0.1250   -0.0000   -0.0000    0.0000    0.0000    0.0000    0.0000
%   -0.0000   -0.0000    0.1250   -0.0000   -0.0000    0.0000   -0.0000   -0.0000
%    0.0000   -0.0000   -0.0000    0.1249   -0.0000    0.0000   -0.0000    0.0000
%    0.0000    0.0000   -0.0000   -0.0000    0.1250   -0.0000   -0.0000    0.0000
%    0.0000    0.0000    0.0000    0.0000   -0.0000    0.1250   -0.0000   -0.0000
%   -0.0000    0.0000   -0.0000   -0.0000   -0.0000   -0.0000    0.1250   -0.0000
%   -0.0000    0.0000   -0.0000    0.0000    0.0000   -0.0000   -0.0000    0.1249
%                                    
% iv) Please see attached figure 2c)
% v)  Please see attached figure 2d)
for m=1:M;
W_c(m,:)=conv(x,H_c(m,:));
X_c(m,:)=W_c(m,1:M:length(W_c(m,:)));
end
for m=1:M;
Z_c(m,:)=zeros(1,M*length(X_c(m,:)));
Z_c(m,1:M:length(Z_c(m,:)))=X_c(m,:);
Y_c(m,:)=conv(Z_c(m,:),G_c(m,:));
end
y_c=zeros(1,length(Y_c(1,:)));
for m=1:M;
y_c=y_c+Y_c(m,:);
end
domega=2*pi/1024;
omega=-pi:domega:pi-domega;
yf1=abs(fftshift(fft(x,1024)));
yf2=M*abs(fftshift(fft(y_c,1024)));
figure(11);
plot(omega,yf1,'Linewidth',1)
axis([-pi pi 0 max(yf1)])
title('fig 3c: DTFT of Gaussian Random Process (mean = 0 with unit power)');
xlabel('\omega (radians/s)');
figure(12);
plot(omega,yf2,'Linewidth',1)
axis([-pi pi 0 max(yf2)])
title('fig 3d: DTFT of output of Gaussian Random Process');
xlabel('\omega (radians/s)');
 
% D)
h0_d=[1 i]; h1_d=[1 -i];
h00_d=[1 0 i];   h10_d=[1 0 -i];
h000_d=[1 0 0 0 i]; h100_d=[1 0 0 0 -i];
H_tmp=conv(h0_d,h00_d); H_d(1,:)=conv(H_tmp,h000_d); G_d(1,:)=H_d(1,:);
H_tmp=conv(h0_d,h00_d); H_d(2,:)=conv(H_tmp,h100_d); G_d(2,:)=-H_d(2,:);
H_tmp=conv(h0_d,h10_d); H_d(3,:)=conv(H_tmp,h000_d); G_d(3,:)=-H_d(3,:);
H_tmp=conv(h0_d,h10_d); H_d(4,:)=conv(H_tmp,h100_d); G_d(4,:)=H_d(4,:);
H_tmp=conv(h1_d,h00_d); H_d(5,:)=conv(H_tmp,h000_d); G_d(5,:)=-H_d(5,:);
H_tmp=conv(h1_d,h00_d); H_d(6,:)=conv(H_tmp,h100_d); G_d(6,:)=H_d(6,:);
H_tmp=conv(h1_d,h10_d); H_d(7,:)=conv(H_tmp,h000_d); G_d(7,:)=H_d(7,:);
H_tmp=conv(h1_d,h10_d); H_d(8,:)=conv(H_tmp,h100_d); G_d(8,:)=-H_d(8,:);
% i)   Please see attached for figure 4a)
for ind=1:M
hf(ind,:)=abs(fftshift(fft(H_d(ind,:),512)));
end
domega=2*pi/512;
omega=-pi:domega:pi-domega;
figure(13);
plot(omega,hf(1,:),omega,hf(2,:), ...
     omega,hf(3,:),omega,hf(4,:), ...
     omega,hf(5,:),omega,hf(6,:), ...
     omega,hf(7,:),omega,hf(8,:));
axis([-pi pi 0 max(hf(1,:))])
legend('H0(w)','H1(w)','H2(w)','H3(w)','H4(w)','H5(w)','H6(w)','H7(w)');
title('fig 4a: Frequency Response of h0[n] through h7[n]');
xlabel('\omega (radians/s)');
% ii)  Please see attached for figure 4b)
for ind=1:M
gf(ind,:)=abs(fftshift(fft(G_d(ind,:),512)));
end
domega=2*pi/512;
omega=-pi:domega:pi-domega;
figure(14);
plot(omega,gf(1,:),omega,gf(2,:), ...
     omega,gf(3,:),omega,gf(4,:), ...
     omega,gf(5,:),omega,gf(6,:), ...
     omega,gf(7,:),omega,gf(8,:));
axis([-pi pi 0 max(gf(1,:))])
legend('G0(w)','G1(w)','G2(w)','G3(w)','G4(w)','G5(w)','G6(w)','G7(w)');
title('fig 4b: Frequency Response of g0[n] through g7[n]');
xlabel('\omega (radians/s)');
% iii) Table 4
% H*H' =
%       8     0     0     0     0     0     0     0
%       0     8     0     0     0     0     0     0
%       0     0     8     0     0     0     0     0
%       0     0     0     8     0     0     0     0 
%       0     0     0     0     8     0     0     0
%       0     0     0     0     0     8     0     0
%       0     0     0     0     0     0     8     0
%       0     0     0     8     0     0     0     8     
% iv) Please see attached figure 4c)
% v)  Please see attached figure 4d)
for m=1:M;
W_d(m,:)=conv(x,H_d(m,:));
X_d(m,:)=W_d(m,1:M:length(W_d(m,:)));
end
for m=1:M;
Z_d(m,:)=zeros(1,M*length(X_d(m,:)));
Z_d(m,1:M:length(Z_d(m,:)))=X_d(m,:);
Y_d(m,:)=conv(Z_d(m,:),G_d(m,:));
end
y_d=zeros(1,length(Y_d(1,:)));
for m=1:M;
y_d=y_d+Y_d(m,:);
end
domega=2*pi/1024;
omega=-pi:domega:pi-domega;
yf1=abs(fftshift(fft(x,1024)));
yf2=abs(fftshift(fft(y_d,1024)));
figure(15);
plot(omega,yf1,'Linewidth',1)
axis([-pi pi 0 max(yf1)])
title('fig 4c: DTFT of Gaussian Random Process (mean = 0 with unit power)');
xlabel('\omega (radians/s)');
figure(16);
plot(omega,yf2,'Linewidth',1)
axis([-pi pi 0 max(yf2)])
title('fig 4d: DTFT of output of Gaussian Random Process');
xlabel('\omega (radians/s)');
