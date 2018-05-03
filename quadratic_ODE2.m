%% Documentation

% Description: 

% Component-wise sparse recovery of the Lorenz96 system using 
% l1 minimization techniques. Uses quadratic tensorized generalized 
% Legendre polynomials as the bounded orthonormal system for recovery.
% Randomly sampled data and subsequent "bursts" used to construct 
% a dictionary data matrix (D), which is used to solve for a sparse 
% coefficient matrix (alpha) through the minimization problem
%
%   min |alpha|_1 s.t. |D*alpha - V|_2 <= sigma
%
% where V is a velocity matrix calculated using finite differences.

% Instructions:
% - To have initial distribution / parameters chosen from defaults,
%   change auto_choose from "off" to "on" in Line 59.
% - To change the noise in X_0, go to Line 53. To change noise in V,
%   go to Line 150.
% - To change display parameters, go to bottom.

% Requires: 
% dictionary_quadratic.m
% Lorenz96_XV.m (and subsequent subscripts)
% DouglasRachford.m (and corresponding subscripts)

% Copyrights: 
% Zachary Singer, Hayden Schaeffer, Department of Mathematical Sciences,
% Carnegie Mellon University, Pittsburgh PA

% Last Edited: 5/3/18

%% Setup

close all; clear all; clc; clf;

% Display options
auto_choose = "off";    % auto-choose distribution
plotLorenz = "off";     % plot pca in 3d
disp_res = "off";        % display results for coefficients
reg_plot = "on";        % display plots of recovered vs. true coeff.
plotTraj = "off";        % display trajectories
cond_num = "on";         % display condition number

% Tuning parameters
n = 50;                     % Dimension (R^n)
N2 = (n^2 + 3*n + 2) / 2;   % Number of basis functions
optEquation = 20;           % Component of system to recover
NumIC = 1;                  % Number of bursts
F = 8.0;                    % Constant in Lorenz 96 System

% Other Parameters

dt = .025; % time step
SizeOfBurst = 1000; % size of each burst
K = SizeOfBurst;

% Distribution Parameters + Initialization

noise_x = 0.00; % X noise parameter
a = 2;
b = 4;
X_0 = (b-a)*rand(n, NumIC) + a;
X_0 = X_0 + noise_x*max(abs(X_0(:)))*randn(n, NumIC); % add X noise

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

% Generate X data, V data, and true V matrix
[Xfull, Vapproximate, Vexact] =  Lorenz96_XV(F, X_0, dt, SizeOfBurst); 

% Generate dictionary matrix D, standard monomial matrix S,
% moment transformation matrix S, and BOS constant max_phi
[D, S, M, max_phi] = dictionary_quadratic(Xfull, K, a, b);

%% Basis Pursuit Denoising Problem

noise_v = 0.00; % V noise
sig_err = 0.95; % 0.95 before


% measurement vector ("b" in "Ax = b")
Vtest = Vapproximate(:, optEquation) + noise_v*...
        max(abs(Vapproximate(:)))*randn(SizeOfBurst*NumIC, 1);

sigma = sig_err.*norm(Vtest-Vexact(:,optEquation), 2);

Dnormalized = D ./ repmat(sqrt(sum(D.^2,1)),size(D,1),1);

% (Parameter Format) DouglasRachford(A, b, sigma, tau, mu, MaxIt, tol)
c = DouglasRachford(Dnormalized, Vtest, sigma, 0.5, 0.5, 5000, 0.01);
c = c ./ (sqrt(sum(D.^2,1)))';

% Transform back to monomial basis
soln = M' * c;

% Tolerance for solution, debiased solution
tol = 0.01;
soln = soln .* (abs(soln) > tol);
debiased_soln = S(:, find(soln)) \ Vtest; % find(soln)

% Line up nonzero indices with debiased solution
disp_debiased = zeros(size(soln));
nonzero_idx = find(soln);
d_len = size(nonzero_idx, 1);
d_idx = 1;
for i = 1 : d_len
    j = nonzero_idx(i);
    disp_debiased(j) = debiased_soln(i);
end

% True Coefficients
c_true_mat = Lorenz96_true_coefficients(n, F);


%% PLOT RESULTS

set(gcf, 'Position', [100, 800, 800, 800])
set(0, 'defaulttextinterpreter', 'latex')

if reg_plot == "on"
    plot(sparse(c_true_mat(:, optEquation)), 'x', 'MarkerSize', 20, ...
         'LineWidth', 16, 'Color', 'red')
    hold on
    plot(soln, 'o', 'MarkerSize', 20, 'LineWidth', 16, 'Color', 'blue')
    set(gca,'fontsize', 24)
    legend({'True Coefficients', 'Recovered Coefficients'}, ...
            'Location', 'northeast', 'FontSize', 24)
end

% Plot Lorenz Trajectory in 3D using PCA
if plotLorenz == "on"
    pcaRes = pca(Xfull');
    figure
    plot3(pcaRes(:, 1), pcaRes(:, 2), pcaRes(:, 3), 'LineWidth', 3)
    grid on
    set(gca,'fontsize', 24)
end

if plotTraj == "on"
    % which components to plot in plane
    comp1 = 1;
    comp2 = 10; 
    wh = b-a; % width and height of rectangle
    
    cm = colormap(jet(NumIC));
    or = [1 0.6 0];
    pu = [0.4 0 1];
    gr = [0 0.6 0.05];
    re = [1 0 0.2];
    bl = [0 0.2 1];
    c_small = [or; pu; gr; re; bl];
    
    figure
    % plot rectangle
    rectangle('Position', [a a wh wh], 'LineWidth', 4)
    hold on
    % plot trajectories
    for i = 1 : NumIC
        start_idx = (i-1)*SizeOfBurst + 1;
        end_idx = i*SizeOfBurst;
        x1 = Xfull(start_idx:end_idx, comp1);
        x2 = Xfull(start_idx:end_idx, comp2);
        plot(x1, x2, 'LineWidth', 4, 'color', c_small(i, :))
        plot(x1(1), x2(1), 'o', 'Markersize', 8, 'LineWidth', 4, 'color', 'red')
    end
    grid on
    set(gca,'fontsize', 18)
    xlabel('$$x_1$$', 'FontSize', 30)
    ylabel('$$x_{10}$$', 'FontSize', 30)
end
    

%% OUTPUT RESULTS

if cond_num == "on"
    display(['Condition number of transformation matrix M is: ', num2str(cond(M))])
end

if disp_res == "on"
    display(['True coefficients using component ', num2str(optEquation)])
    sparse(c_true_mat(:,optEquation))

    display(['Recovered coefficients using component  ', num2str(optEquation)])
    sparse(soln)

    display(['Debiased coefficients using component  ', num2str(optEquation)])
    sparse(disp_debiased)

    display(['K boundedness coefficient:  ', num2str(max_phi)])   
 end
