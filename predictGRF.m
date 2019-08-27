%for test purpose
%predicting GRF from the vertical dispalcement
clear;close all;
%X(t) = A*sin(omega*t)+B*cos(omega*t)+C was assumed with the contraints 
%V(0) = V0 & A(0) = -g
T = 0.1;%in sec
g = 9.80;
m = 66.4;%kg from Takeshita et al.
b_a = 3;%lever arm ratio
dt = T/100;
t = 0:dt:1.1*T;
omega = 2*pi/T;
V0 = -1;
x = V0/omega*sin(omega*t) + g/omega^2 *cos(omega*t);
v = V0*cos(omega*t) -g/omega*sin(omega*t);
v2 = diff(x)/dt;
a = -V0*omega*sin(omega*t) - g*cos(omega*t);
th = g/omega^2;
[tTakeoff, idxTakeoff] = findTcontact(t,x,T,th)

a2 = diff(v)/dt;
GRF = m*(a + g);
F_tend = GRF*b_a;

nfig = 5;
subplot(nfig,1,1)
plot(t,x)
hold on
plot(tTakeoff, x(idxTakeoff),'ro')
subplot(nfig,1,2)
plot(t,v)
hold on
plot(t(2:end),v2,'rx')
subplot(nfig,1,3)
plot(t,a)
hold on
plot(t(2:end),a2,'bx')
subplot(nfig,1,4)
plot(t,GRF)
subplot(nfig,1,5)
plot(t,F_tend)



