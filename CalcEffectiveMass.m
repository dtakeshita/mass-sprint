clear;
f_n = 2.84;%Hz
omega_n = 2*pi*f_n;
K_mm = 165.7;%N/mm
g = 9.8;
K = K_mm*1e3;%N/m
M = K/omega_n^2
%% 
F0 = 1800;%N, from Takeshita et al. 
m = 66.4;%kg
M2 = F0^2/(m*g^2)