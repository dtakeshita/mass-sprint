function dydt = focrceoscifcn(t,y,ft,X_CC,M,K,C,X_SEC0,g)
X_CCnew = interp1(ft,X_CC,t);
dydt = [y(2);g-K/M*(y(1)-X_SEC0-X_CCnew) - C/M*y(2)];
