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

% Setting up M=8 channel DFT filter bank 
M = 8; 
h0 = [1 1] / sqrt(2); 
h1 = [1 -1] / sqrt(2); 
h00 = [1 0 1] / sqrt(2); 
h10 = [1 0 -1] / sqrt(2); 
h000 = [1 0 0 0 1] / sqrt(2); 
h100 = [1 0 0 0 -1] / sqrt(2); 


H_temp = conv(h0, h00); 
H(1,:) = conv(H_temp, h000); G(1,:) = H(1, :); 
H(2,:) = conv(H_temp, h100); G(2,:) = -H(2, :); 

H_temp = conv(h0, h10);
H(3,:) = conv(H_temp, h000); G(3,:) = -H(3, :); 
H(4,:) = conv(H_temp, h100); G(4,:) = H(4, :); 

H_temp = conv(h1, h00);
H(5,:) = conv(H_temp, h000); G(5,:) = -H(5, :); 
H(6,:) = conv(H_temp, h100); G(6,:) = H(6, :); 

H_temp = conv(h1, h10);
H(7,:) = conv(H_temp, h000); G(7,:) = H(7, :); 
H(8,:) = conv(H_temp, h100); G(8,:) = -H(8, :); 

% i) All corresponding DTFT's H_m(w) 

h_m(1,:)=abs(fftshift(fft(H(1,:),512)));
h_m(2,:)=abs(fftshift(fft(H(2,:),512)));
h_m(3,:)=abs(fftshift(fft(H(3,:),512)));
h_m(4,:)=abs(fftshift(fft(H(4,:),512)));
h_m(5,:)=abs(fftshift(fft(H(5,:),512)));
h_m(6,:)=abs(fftshift(fft(H(6,:),512)));
h_m(7,:)=abs(fftshift(fft(H(7,:),512)));
h_m(8,:)=abs(fftshift(fft(H(8,:),512)));

domega = 2 * pi / 512; 
omega = -pi:domega:pi-domega; 

figure(1) %Figure 1(a)  
plot(omega, h_m(1,:), ... 
     omega, h_m(2,:), ... 
     omega, h_m(3,:), ... 
     omega, h_m(4,:), ... 
     omega, h_m(5,:), ... 
     omega, h_m(6,:), ... 
     omega, h_m(7,:), ... 
     omega, h_m(8,:)); 
axis([-pi pi 0 3]); 
title('Figura 1(a): All corresponding DTFT''s H_m(\omega)')
ylabel('H_m(\omega)'); 
xlabel('Omega, \omega (rad/sec)');
legend('H_1(\omega)', 'H_2(\omega)', ...
       'H_3(\omega)', 'H_4(\omega)', ...
       'H_5(\omega)', 'H_6(\omega)', ...
       'H_7(\omega)', 'H_8(\omega)'); 
grid on

%ii) 8x8 matrix HH^H 
table1 = H*H' 

