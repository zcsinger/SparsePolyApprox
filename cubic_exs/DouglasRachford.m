function x = DouglasRachford(A,b,sigma,tau,mu,MaxIt,tol)

% ============================================================
% Reference:
%   http://www.numerical-tours.com/matlab/sparsity_5_dictionary_learning_denoising/
%
% Inputs:
%   A = m*n matrix
%   b = m*1 vector
%   sigma = constrain parameter
%   tau = step size
%   mu = convergence parameter
%   MaxIt = maximum number of iterations allowed
%   tol = tolerance
% Outputs:
%   u = min_{x} ||x||_1
%       sub to  ||Ax-b||_2 <= sigma
%
% Date: Wednesday, July 19, 2017
% ============================================================

% Initialization
T = 1;
N = size(A,2);
M = length(b);
x = zeros(N,1);
x1 = zeros(N,1);
u1 = zeros(M,1);
error = tol+1;

pAlower = chol(eye(N) + A'*A,'lower');
mu1 = 1-mu;

while T<=MaxIt %&& error>=tol
    T = T+1;
    % first step
    [p,q] = rProxF(x1,u1,tau,b,sigma);
    [p,q] = rProxG(p,q,A,b,pAlower);
    x1 = mu1*x1 + mu*p;
    u1 = mu1*u1 + mu*q;
   
    % second step
    [xnew,~] = ProxF(x1,u1,tau,b,sigma);
    
    % update
    error = norm(x-xnew);
    x = xnew;
    
end

end