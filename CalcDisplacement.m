clear; close all;

V0 = 0;
g = 9.8;
T = 0.1:0.01:1;
omega = 2*pi./T;
A = sqrt((V0./omega).^2 + g^2./omega.^4);
plot(T,A)