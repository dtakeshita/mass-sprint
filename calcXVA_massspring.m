function [x,v,a] = calcXVA_massspring(k, m, v0)
    % calculating x, v, a for "hopping" of mass-spring system
    %ma = mg - ky, then y(0) = 0 (no force from the spring at the landing)
    % & v(0) = v0
    omega = sqrt(k/m);
    x = v0/omega*sin(omega*t) - m*g/k*cos(omega*t) + m*g/k;
    v = v0*cos(omega*t) + g/omega*sin(omega*t);
    a = -v0*omega*sin(omega*t) + g*cos(omega*t);

end