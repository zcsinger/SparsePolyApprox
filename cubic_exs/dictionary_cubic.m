function [phiX, stanX, M, max_phi] = dictionary_cubic(U, SizeOfBurst, a, b)
% Description: Construct the cubic dictionary matrix phiX containing 
% all multivariate monomials up to degree three for the Lorenz 96
% Input: U = [x1(t1) x2(t1) .... xn(t1)
%             x1(t2) x2(t2) .... xn(t2)
%                    ......
%             x1(tm) x2(tm) .... xn(tm)]
%        option = [] (monomial) or 'legendre'

%        SizeOfBurst : number of bursts
%        a, b : domain interval

m_dim = size(U,1); % number of measurements
n = size(U,2); % dimension of the ODE

N1 = n + 1;
N2 = (n^2 + 3*n + 2) / 2;
N3 = (n^3 + 6*n^2 + 11*n + 6) / 6;

% indexing function (for quadratic and cubic)
g = @(i) 0.5 * (-i^2 + i + 2);
f = @(i, j) i*n + (j + g(i));

%%%%% DATA AND CONSTANT INITIALIZATION %%%%%

% retrieve data from U (XFull)

t0_data = U(1:SizeOfBurst:end, :); % take every SizeOfBurst-st row
t0_data_vec = t0_data(:);

% Calculate central moments and mean through data

mu_1 = mean(t0_data_vec);
mu_2 = moment(t0_data_vec, 2);
mu_3 = moment(t0_data_vec, 3);
mu_4 = moment(t0_data_vec, 4);
mu_5 = moment(t0_data_vec, 5);
mu_6 = moment(t0_data_vec, 6);

alpha = mu_1^2 + mu_2;
beta = mu_1^3 + 3*mu_1*mu_2 + mu_3;
gamma = mu_1^4 + 4*mu_1*mu_3 + 6*mu_1^2*mu_2 + mu_4;
eta = mu_5 + 5*mu_1*mu_4 + 10*mu_1^2*mu_3 + 10*mu_1^3*mu_2 + mu_1^5;
lambda = mu_6 + 6*mu_1*mu_5 + 15*mu_1^2*mu_4 + 20*mu_1^3*mu_3 + ...
            15*mu_1^4*mu_2 + mu_1^6;
     
% Denominator monster constants
        
kappa = gamma - beta^2 / alpha  - alpha^2 + 2*beta*mu_1;       
        
m = lambda - eta^2 / gamma - gamma^2 / alpha - beta^2 + 2*beta / alpha ...
     * (eta + gamma * mu_1) + 2*alpha*beta*eta / gamma;

%%%%% MATRIX OF LINEAR TRANS CONSTRUCTION %%%%% 

M = zeros(N3);

% Create index matrix (3 dimensions)
% M3 = zeros(n, n, n); % uncomment if you want the index matrix

ctr = N2; %ctr adds one initially
for i = 1:n
    for j = i:n
        for k = j:n
            ctr = ctr + 1;
            % M3(k, i, j) = ctr; % uncomment if you want the index matrix
            ij_idx = f(i, j);
            jk_idx = f(j, k);
            ik_idx = f(i, k);
            if i == j
                if j == k
                    % x_i^3
                    assert(ij_idx == ik_idx);
                    M(ctr, ctr) = 1/sqrt(m);
                    M(ctr, 1) = -beta/sqrt(m);
                    M(ctr, i + 1) = -gamma/(alpha * sqrt(m));
                    M(ctr, jk_idx) = - eta / (gamma * sqrt(m));
                else
                % x_i^2 x_k
                assert(ik_idx == jk_idx);
                M(ctr, ctr) = 1/sqrt(mu_2 * kappa);
                M(ctr, 1) = alpha * mu_1/sqrt(mu_2 * kappa);
                M(ctr, i + 1) = beta * mu_1/(alpha * sqrt(mu_2 * kappa));
                M(ctr, k + 1) = -alpha / sqrt(mu_2 * kappa);
                M(ctr, ij_idx) = -mu_1 / sqrt(mu_2 * kappa);
                M(ctr, jk_idx) = - beta / (alpha * sqrt(mu_2 * kappa));
                end
                
            % i ~= j
            elseif j == k
                % x_i x_j^2
                assert(ij_idx == ik_idx);
                M(ctr, ctr) = 1/sqrt(mu_2 * kappa);
                M(ctr, 1) = alpha * mu_1/sqrt(mu_2 * kappa);
                M(ctr, i + 1) = -alpha / sqrt(mu_2 * kappa);
                M(ctr, k + 1) = beta * mu_1 / (alpha * sqrt(mu_2 * kappa));
                M(ctr, jk_idx) = -mu_1 / sqrt(mu_2 * kappa);
                M(ctr, ij_idx) = -beta / (alpha * sqrt(mu_2 * kappa));
                
            else
                % x_i x_j x_k
                M(ctr, ctr) = 1/mu_2^(3/2);
                M(ctr, 1) = -mu_1^3 / mu_2^(3/2);
                M(ctr, i + 1) = mu_1^2 / mu_2^(3/2);
                M(ctr, j + 1) = mu_1^2 / mu_2^(3/2);
                M(ctr, k + 1) = mu_1^2 / mu_2^(3/2);
                M(ctr, ij_idx) = -mu_1 / mu_2^(3/2);
                M(ctr, jk_idx) = -mu_1 / mu_2^(3/2);
                M(ctr, ik_idx) = -mu_1 / mu_2^(3/2);
            end

                
        end
    end
end

for i = 1 : n
    for j = i : n
        ij_idx = f(i, j);
        if i == j
            M(ij_idx, ij_idx) = 1 / sqrt(kappa);
            M(ij_idx, i + 1) = -beta / (alpha * sqrt(kappa));
            M(ij_idx, 1) = -alpha / sqrt(kappa);
        else
            M(ij_idx, ij_idx) = 1 / mu_2;
            M(ij_idx, i + 1) = -mu_1 / mu_2;
            M(ij_idx, j + 1) = -mu_1 / mu_2;
            M(ij_idx, 1) = mu_1^2 / mu_2;
        end
    end
end


for i = 1:n
    M(i + 1, i + 1) = 1 / sqrt(mu_2);
    M(i + 1, 1) = -mu_1 / sqrt(mu_2); % fixed mistake here
end

M(1, 1) = 1; 

%%%%%% MONOMIAL TERM CALCULATION %%%%% 

stanX = zeros(m_dim, N3);

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

% Cubic Terms
quad_len = (n^2 + n) / 2;
quad_start = n + 2;
quad_stop = quad_start + quad_len;
quad_mat = stanX(:, quad_start : quad_stop); 
% g = @(i) 0.5 * (-i^2 + i + 2);

cubic_idx = N2 + 1;
for i = 1 : n
    i_len = (i-1)*n + g(i-1); % tells you which index to start
    % Example: For i = 1, i_len = 1. For i = 4, i_len = 3n-2.
    % In general: i_len(1) = 1, i_len(i) = i_len(i-1) + (n-i-2)
    new_cubic_idx = cubic_idx + (quad_len - i_len + 1);
    i_mat = repmat(U(:, i), 1, quad_len - i_len + 1);
    cubic_mat = i_mat .* quad_mat(:, i_len : quad_len);
    stanX(:, cubic_idx : new_cubic_idx-1) = cubic_mat;
    cubic_idx = new_cubic_idx;
end

%Mnormalized = M ./ repmat(sqrt(sum(M.^2,1)),size(M,1),1);
%M = Mnormalized;
phiX = stanX * M'; 

% phi calculations 

phi2_r1 = (mu_1-a)^2 / mu_2;
phi2_r2 = (mu_1-b)^2 / mu_2;
phi2_r3 = abs((a-mu_1)*(b-mu_1)) / mu_2;

phi3_r1 = abs(a^2 - beta/alpha * a - alpha) / sqrt(kappa);
phi3_r2 = abs(b^2 - beta/alpha * b - alpha) / sqrt(kappa);
phi3_r3 = (beta^2 / (4*alpha^2) + alpha) / sqrt(kappa);

% find roots for cubic, only evaluate if real
p4_discrm = eta^2/gamma^2 + 3*gamma/alpha;
phi4_r1 = 0;
phi4_r2 = 0;
if (p4_discrm >= 0)
    p4_root1 = eta/(3 * gamma) - sqrt(p4_discrm) / 3;
    p4_root2 = eta/(3 * gamma) + sqrt(p4_discrm) / 3;
    phi4_r1 = (p4_root1^3 - eta/gamma * p4_root1^2 - ...
            gamma/alpha * p4_root1 - beta) / sqrt(m);
    phi4_r2 = (p4_root2^3 - eta/gamma * p4_root2^2 - ...
            gamma/alpha * p4_root2 - beta) / sqrt(m);
end

phi32_r1 = 1/sqrt(mu_2*kappa) * (a^2 - beta/alpha * a - alpha) * (a - mu_1);
phi32_r2 = 1/sqrt(mu_2*kappa) * (b^2 - beta/alpha * b - alpha) * (b - mu_1);
phi22_r1 = 1/mu_2^1.5 * (a - mu_1)^3;
phi22_r2 = 1/mu_2^1.5 * (b - mu_1)^3;

crits = [phi2_r1 phi2_r2 phi2_r3 phi3_r1 phi3_r2 phi3_r3 ...
         phi4_r1 phi4_r2 phi32_r1 phi32_r2 phi22_r1 phi22_r2];

max_phi = max(crits);