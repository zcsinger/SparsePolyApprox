function [x1,u1] = ProxG(x,u,A,b,pAlower)
y = x + A'*u;
x1 = pAlower'\(pAlower\y);
%x1 = pAlower * (b + A'*u);
u1 = A*x1;
end