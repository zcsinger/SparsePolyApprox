function [phiX, stanX, M, max_phi] = dictionary_quadratic(U, SizeOfBurst, a, b)
%% Documentation
% Goal: Construct the quadratic dictionary matrix phiX containing 
% all multivariate monomials up to degree two for the Lorenz 96 system

% Input: U = [x1(t1) x2(t1) .... xn(t1)
%             x1(t2) x2(t2) .... xn(t2)
%                    ......
%             x1(tm) x2(tm) .... xn(tm)]

%        SizeOfBurst : number of bursts
%        a, b : domain interval

% Output: 
%   phiX : dictionary matrix D(X) with respect to orthonormal poly BOS
%   stanX : dictionary matrix with respect to standard polys
%   M : transformation matrix between orthonormal and standard polys
%   max_phi : calculated K boundedness constant for orthonormal poly BOS

% Copyrights: 
% Zachary Singer, Hayden Schaeffer, Department of Mathematical Sciences,
% Carnegie Mellon University, Pittsburgh PA

% Last Edited: 5/3/18

%% Data Initialization

m_dim = size(U, 1); % number of measurements
n = size(U, 2); % dimension of the ODE

N1 = n + 1;
N2 = (n^2 + 3*n + 2) / 2;

% Calculate central moments and mean through data

mu_1 = mean(U);
mu_2 = moment(U, 2);
mu_3 = moment(U, 3);
mu_4 = moment(U, 4);

alpha = mu_1.^2 + mu_2;
beta = mu_1.^3 + 3*mu_1.*mu_2 + mu_3;
gamma = mu_1.^4 + 4*mu_1.*mu_3 + 6*mu_1.^2.*mu_2 + mu_4;
        
kappa = gamma - beta.^2 ./ alpha  - alpha.^2 + 2*beta.*mu_1; % denom const
        
%% Matrix of Linear Transformation

M = zeros(N2);

% Indexing functions
g = @(i) 0.5 * (-i^2 + i + 2);
f = @(i, j) i*n + (j + g(i));

for i = 1 : n
    for j = i : n
        ij_idx = f(i, j);
        if i == j
            M(ij_idx, ij_idx) = 1 / sqrt(kappa(i));
            M(ij_idx, i + 1) = -beta(i) / (alpha(i) * sqrt(kappa(i)));
            M(ij_idx, 1) = -alpha(i) / sqrt(kappa(i));
        else
            M(ij_idx, ij_idx) = 1 / sqrt(mu_2(i)*mu_2(j));
            M(ij_idx, i + 1) = -mu_1(i) / sqrt(mu_2(i)*mu_2(j));
            M(ij_idx, j + 1) = -mu_1(j) / sqrt(mu_2(i)*mu_2(j));
            M(ij_idx, 1) = mu_1(i)*mu_1(j) / sqrt(mu_2(i)*mu_2(j));
        end
    end
end


for i = 1:n
    M(i + 1, i + 1) = 1 / sqrt(mu_2(i));
    M(i + 1, 1) = -mu_1(i) / sqrt(mu_2(i)); % fixed mistake here
end

M(1, 1) = 1; 

%% Monomial Term Calculation

stanX = zeros(m_dim, N2);

stanX(:, 1) = ones(m_dim, 1);
stanX(:, 2 : n+1) = U;

% Quadratic Terms
quad_idx = N1 + 1;
lin_mat = stanX(:, 2:n+1);

for i = 1 : n
    new_quad_idx = quad_idx + (n-i+1);
    i_mat = repmat(U(:, i), 1, n-i+1); % x_i repeated
    i_quad_mat = i_mat .* lin_mat(:, i:n);
    stanX(:, quad_idx : new_quad_idx-1) = i_quad_mat;
    quad_idx = new_quad_idx;
end

phiX = stanX * M'; % OUTPUT

%% K Boundedness Constant Calculations

phi2_r1 = (mu_1-a).^2 ./ mu_2;
phi2_r2 = (mu_1-b).^2 ./ mu_2;
phi2_r3 = abs((a-mu_1).*(b-mu_1)) ./ mu_2;

phi3_r1 = abs(a^2 - beta./alpha * a - alpha) ./ sqrt(kappa);
phi3_r2 = abs(b^2 - beta./alpha * b - alpha) ./ sqrt(kappa);
phi3_r3 = (beta.^2 / (4*alpha.^2) + alpha) ./ sqrt(kappa);

crits = [phi2_r1 phi2_r2 phi2_r3 phi3_r1 phi3_r2 phi3_r3];

max_phi = max(max(crits));