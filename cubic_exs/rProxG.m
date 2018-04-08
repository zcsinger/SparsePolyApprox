function [x1,u1] = rProxG(x,u,A,b,pAlower)
[x1,u1] = ProxG(x,u,A,b,pAlower);
x1 = 2*x1-x;
u1 = 2*u1-u;
end