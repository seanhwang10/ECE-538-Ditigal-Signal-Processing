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




