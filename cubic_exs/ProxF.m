function [x1,u1] = ProxF(x,u,tau,b,sigma)
x1= max(abs(x)-tau,0) .* sign(x);
u1 = b + (u-b) * min(sigma/norm(u-b,2),1);
end