%-----------------------------%
%         Sean Hwang          %
% ECE 538 - MATLAB Project 1 %
%     Source Code File        %
%-----------------------------%

% Parameters to be changed 

a2 = -1; 
D2 = 21; 

v = randn(1,200) * 0.7; %variance == std.div, since 1. 

% Generating sequence for M=127 using given M=15 file
register=[1 0 0 0 0 0 0];
for ri=1:127,
m127(ri)=register(1,7);
register(2:7)=register(1:6);
register(1,1)=rem((register(1,1)+m127(1,ri)),2);
end
m127=2*m127-1; 

x = [m127 zeros(1, 200-127)];
x20 = [zeros(1,20) m127 zeros(1, 200-127-20)]; 
%x[n-20] 
xD2 = [zeros(1,D2) m127 zeros(1, 200-127-D2)];  
%x[n-D] 

% Parameters not to be changed 

n = 0:199; 

y = x20 + a2*xD2 + v; %y[n] = x[n-20] + a2x[n-D2] + v[n]  

ryx = conv(y,x(end:-1:1)); 
[a,bound] = xcorr(y,x); 

%Plot generation 

%Plot i
subplot(3,1,1); 
plot (m127); 
xlim([0 127]); 
grid on; 
xlabel('n'); 
ylabel('x[n]'); 
title('Plot (i): Plot of values of x[n], for n = 0:127'); 
 
%Plot ii 
subplot(3,1,2); 
plot (y); 
grid on; 
xlabel('n'); 
ylabel('y[n]');
title('Plot (ii): Plot of values of y[n], for n = 0:199'); 

%plot iii 
subplot(3,1,3); 
plot (bound,ryx); 
xlim([0 59]);
grid on; 
xlabel ('Lag'); 
ylabel ('Cross Correlation'); 
title ('Plot (iii): Plot of Cross correlation ryx(ell)'); 



