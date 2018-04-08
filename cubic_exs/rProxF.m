function [x1,u1] = rProxF(x,u,tau,b,sigma)
[x1,u1] = ProxF(x,u,tau,b,sigma);
x1 = 2*x1-x;
u1 = 2*u1-u;
end