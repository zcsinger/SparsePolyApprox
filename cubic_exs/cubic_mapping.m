% Function Mapping Test
% (CUBIC DISCRETE TIME DYNAMICAL SYSTEM RECOVERY)

% Sample points and iteratively apply a function mapping
% (discrete time dynamical system) before creating dictionary matrix.
% Uses cubic model for tensorized basis functions.

% Zachary Singer, 4/8/17

%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS / SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clf; clear all;

% dimension and number of initializations
n = 15;
optEquation = 5;
NumIC = 150;
SizeOfBurst = 4;
m = SizeOfBurst;    
M = NumIC * SizeOfBurst;

% initial distribution
a = -1;
b = 1;
anor = a*ones(NumIC*SizeOfBurst, n);
bnor = b*ones(NumIC*SizeOfBurst, n);
X_0 = (b-a)*rand(NumIC, n) + a;

        rng default  % For reproducibility
        p = haltonset(n, 'Skip',1e3,'Leap', 1e2); % retains every 101st val
        p = scramble(p, 'RR2'); % reverse radix scrambling
        
        X_0 = (b-a)*net(p, NumIC) + a;

auto_choose = "off";

if auto_choose == "on"

int_range = input(['Range (Type): 1 : [-1, 1], 2 : [-2, 2], 3 : [0, 2],', ...
            '\n 4 : [-2, 0], 5 : [-4, 2], 6 : [-2, 4], 7  : [-4, 4]', ...
            '\n 8 : [-6, -2], 9 : [2, 6], 10 : [-0.5, 0.5]: ']);
dist = input(['Distribution (Type): 1 : Uniform, 2 : Normal,', ...
                '\n 3 : Beta, 4 : Quasi-Uniform: ']);
trans = input(['Transformation (Type): 1 : None, 2 : Squared, ', ...
         '\n 3 : Cubed, 4 : Cosine, 5 : Sine: ']);
switch int_range
    case 1
        a = -1;
        b = 1;
    case 2
        a = -2;
        b = 2;
    case 3
        a = 0;
        b = 2;
    case 4
        a = -2;
        b = 0;
    case 5
        a = -4;
        b = 2;
    case 6
        a = -2;
        b = 4;
    case 7
        a = -4;
        b = 4;
    case 8
        a = -6;
        b = -2;
    case 9
        a = 2;
        b = 6;
    case 10
        a = -0.5;
        b = 0.5;
    otherwise
        disp('Error: Input 1-10.')
end

switch dist
    case 1
        X_0 = (b-a)*rand(NumIC, n) + a;
    case 2
        X_0 = (b-a)*randn(NumIC, n) + a;
    case 3
        alp = input('Alpha (> 0): ');
        bet = input('Beta (> 0): ');
        X_0 = (b-a)*betapdf(rand(NumIC, n), alp, bet) + a;
    case 4
        rng default  % For reproducibility
        p = haltonset(n, 'Skip',1e3,'Leap', 1e2); % retains every 101st val
        p = scramble(p, 'RR2'); % reverse radix scrambling
        
        X_0 = (b-a)*net(p, NumIC) + a;
    otherwise
        disp('Error: Input 1-4.')
end

switch trans
    case 1
        X_0 = X_0;
    case 2
        X_0 = X_0.^2;
    case 3
        X_0 = X_0.^3;
    case 4
        X_0 = cos(X_0);
    case 5
        X_0 = sin(X_0);
    otherwise
        disp('Error. Input 1-5.')
end

end

% Coefficients for mapping (defined directly in for loop later)

c1 = .1;
c2 = .1;
c3 = .1;
c4 = .1;
k = optEquation;
% X(k)
map_fun = @(X) - c1*X(k).*X(k+1).^2 + c2*X(k).^3 + c3*X(k+1).^3;

% true matrix
N2 = (n^2 + 3*n + 2) / 2;
N3 = (n^3 + 6*n^2 + 11*n + 6) / 6;

%%%%%%%%%% GENERATE INDEXING MATRIX AND TRUE COEFFICIENTS MATRIX %%%%%%%%%%

c_true_mat = zeros(N3, n);

% indexing function (for quadratic and cubic)
g = @(i) 0.5 * (-i^2 + i + 2);
h = @(i, j) i*n + (j + g(i));

coef = zeros(n, n, n);

ctr = N2;
for i = 1 : n
    for j = i : n
        for k = j : n
            ctr = ctr + 1;
            coef(i, j, k) = ctr;
        end
    end
end

% Generate true matrix

for i = 1 : n-1
    cubic1 = coef(i, i, i);
    cubic2 = coef(i, i, i+1);
    cubic3 = coef(i, i+1, i+1);
    cubic4 = coef(i+1, i+1, i+1);
    c_true_mat(i+1, i) = 1; % linear
    c_true_mat(cubic1, i) = c1;
    c_true_mat(cubic2, i) = -c2;
    c_true_mat(cubic3, i) = -c3;
    c_true_mat(cubic4, i) = c4;
end

% boundary case : i = n
n_cubic1 = coef(n, n, n);
n_cubic2 = coef(1, 1, n);
n_cubic3 = coef(1, 1, 1);
c_true_mat(n+1, n) = 1;
c_true_mat(n_cubic1, n) = -c1;
c_true_mat(n_cubic2, n) = c2;
c_true_mat(n_cubic3, n) = c3;

% Function Mapping
Xfull = []; % 0 to n-1
Xsome = []; % 1 to n
X_new = zeros(m, n);

for i = 1 : NumIC
    X_dat = X_0(i, :);
    Xfull = [Xfull ; X_dat];
    for j = 1 : SizeOfBurst 
        % apply function map
        for k = 1 : n - 1
              X_new(j, k) = X_dat(k) + c1*X_dat(k)^3  - c2*X_dat(k)^2*X_dat(k+1) ...
                 - c3*X_dat(k)*X_dat(k+1)^2 + c4*X_dat(k+1)^3;
             %  X_new(j, k) = 0.1*X_dat(k)^2 + .01*sin(X_dat(k));
        end
             X_new(j, n) = X_dat(n) + c1*X_dat(n)^2  - c2*X_dat(n)^2*X_dat(1) ...
                 - c3*X_dat(n)*X_dat(1)^2 + c4*X_dat(1)^3;
        X_dat = X_new(j, :);
    end
    Xfull = [Xfull ; X_new(1:m-1, :)];
    Xsome = [Xsome ; X_new];
    
end

%%%%%%%%%% GENERATE DICTIONARY AND SOLVE OPTIMIZATION PROBLEM %%%%%%%%%%

[D, S, M, max_phi] = dictionary_cubic(Xfull, SizeOfBurst, a, b);       

noise = 0.01;
b1 = Xsome(:, optEquation) + noise*max(Xsome(:, optEquation)) ...
     * randn(SizeOfBurst*NumIC, 1);

sigma = 0.2; 

Dnormalized = D ./ repmat(sqrt(sum(D.^2,1)),size(D,1),1);

c2 = DouglasRachford(Dnormalized, b1, sigma, 0.1, 0.1, 5000, 0.01);
c2 = c2 ./ (sqrt(sum(D.^2,1)))';

soln2 = M' * c2;

c = c2;
soln = soln2;
soln = soln .* (abs(soln) > 0.001);
debiased_soln = S(:, find(soln)) \ b1;


%%%%%%%%%%% PLOT AND DISPLAY RESULTS %%%%%%%%%%

set(gcf, 'Position', [100, 800, 600, 600])
set(0, 'defaulttextinterpreter', 'latex')

plot(soln2, 'o', 'MarkerSize', 8, 'LineWidth', 3)
hold on
plot(sparse(c_true_mat(:, optEquation)), 'x', 'MarkerSize', 8, 'LineWidth', 3, 'Color', 'red')

legend({'Recovered Solutions', 'True Solutions'}, 'Location', 'northeast', 'FontSize', 24)

display(['True coefficients using component ', num2str(optEquation)])
sparse(c_true_mat(:,optEquation))

display(['Recovered coefficients using component  ', num2str(optEquation)])
sparse(soln)

display(['Debiased coefficients using component  ', num2str(optEquation)])
sparse(debiased_soln)

display(['K boundedness coefficient:  ', num2str(max_phi)])    

    
