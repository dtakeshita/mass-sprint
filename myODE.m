function dydt = myode(t,y,ft,f)
f = interp1(ft,f,t);
dydt = [y(2);-10*y(1) - 5*y(2) + f];