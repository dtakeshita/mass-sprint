clear; 
close all;
Tc = 0.15;
m = 66.4;
M = 520.38;%effective mass, see CalcEffectiveMass.m
K = 1.657e5;%N/m - 165.7 N/mm
V0 = 0.6;%m/sec - initial velocity of CM
g = 9.80;
%
R0 = 1/3;%R = r_b/r_a (=Ankle joint to COP/)
V0 = 0.6;
Tc = 0.15;
f = @(R)abs(V0*R*cos(R*sqrt(K/m)*Tc/2)+g/(R*sqrt(K/m))*sin(R*sqrt(K/m)*Tc/2));
R = fminsearch(f,R0);
1/R*4.5