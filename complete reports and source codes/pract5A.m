
%----------Practical No.5 A------------------------------------%
%---------To perform linear convolution of given sequences---------%
clc
close all
clear all


x= [exp(j*pi*0) exp(j*pi*1) exp(j*pi*2) exp(j*pi*3) exp(j*pi*4) exp(j*pi*5) exp(j*pi*6) exp(j*pi*7)]
h= [(-1)^0 (-1)^1 (-1)^2 (-1)^3 ]

m=length(x);
n=length(h);
subplot(4,1,1)
stem(1:m,x,'fill','r')
grid on;
title('input sequence, x(n)=')
xlabel('time n------>')
ylabel('amplitude----->')

subplot(4,1,2)
stem(1:n,h,'fill','r')
grid on;
title('impulse sequence, h(n)=')
xlabel('time n------>')
ylabel('amplitude----->')
%------linear convolution using inbuilt command------------------------%
v=conv(x,h)
l=m+n-1;
subplot(4,1,3)
stem(1:l,v,'fill','r')
grid on;
title('output sequence using inbuilt command, v(n)=')
xlabel('time n------>')
ylabel('amplitude----->')
%--------linear convolution using 'for' loop------------------------------%
X=zeros(1,l);
H=zeros(1,l);
X(1:m)=x;
H(1:n)=h;
for i=1:l
    Y(i)=0;
    for j=1:i
        Y(i)=Y(i)+X(j)*H(i-j+1);
    end
end
Y
subplot(4,1,4)
stem(1:l,Y,'fill','r')
grid on;
title('output sequence using loop, Y(n)=')
xlabel('time n------>')
ylabel('amplitude----->')

