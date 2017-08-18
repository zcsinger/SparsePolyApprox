% Solve a particular vegetative system through continuation methods

% Zachary Singer, University of Minnesota Twin Cities, 8/9/17

% Differences between continuation.m and continuation_h2.m : 
%
% - Looking for other heteroclinic solution. Boundary conditions change.
% - Format is (previous boundary condition), current boundary condition :
%   - (b(N) - lamb_neg*v(N)), b(1) - lamb_pos*v(1)
%   - (dot([b(1) v(1)] - [b_plus s*b_plus], stab_eval)), 
%       dot([b(N) v(N)] - [b_minus s*b_minus], unstab_eval)
%
% - Different initial interfaces (opposite in sign)

%%%%%%%%%%%%%%%%%%% Model and Parameters: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scaled system: 
% b_x = v - s*b
% v_x = -(theta - v)*b^2 + b

% theta - varied parameter (want to find theta in terms of s)

%%%%%%%%%%%%%%%%%%%% Specifications and Guide %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Created to use tangent arclength continuation to find a more accurate
% relationship between front speed s and parameter theta in model above.

% Relies on cont_df.m for function and jacobian calculations, and
% cont_dz.m for continuation method calculations

% How to use:

% There are modes you can choose for particular tasks:
% plotting  - see plots of space x vs. biomass b and space x vs. water w
%

% Mode 1 (Debug):
%   Debugs by testing to make sure that the Jacobian expression is correct.
%   Tests different magnitudes to show quadratic error.

% Mode 2 (Initial Newton):
%   A single newton iteration to find the first solution to the nonlinear
%   equation given in cont_df.m.
%   F takes in u = [b; v; s] (2*N+1 variables). Has (N-1) + (N-1) eqns.
%   related to b and v, as well as 2 boundary cond. and 1 phase cond.
%   Plots Space x vs. Biomass b just to see interface, not too instructive. 

% Mode 3 ("Baby Continuation" in Theta):
%   Sets up initial newton iteration as in Mode 2, but linearly varies
%   Theta in one direction and plots the interface as well as 
%   Theta vs. Front Speed s 

% Mode 4 (Tangent Arclength Continuation):
%   First will find single solution using fsolve, and then will use 
%   continuation in theta with initial direction dz_old (line ~225).
%   Then will use tangent arclength continuation to find new directions
%   and new solutions. Plots space x vs. biomass b (not too useful) and
%   theta vs. front speed s. Writes data to file.

% Mode 5 (Data Transformation and Plotting w_plus vs s):
%   Given data found in Mode 4, transforms the data by using scalings
%   from our model back to modified Klausmeier, and plots the curves for
%   the height of the water front w_plus vs. front speed s.
%   Also can use data from veg_model_80317.m to compare curves.

% Suggested Ranges of Variables:
% dx : 0.1 to 0.5 (smaller better)
% theta / s : a list of converging initial conditions given (picky)
% b_old / b_oldx : phase condition, looks like step function 
% dz_old (in mode 4) : initial direction solely in direction of theta
%                      (i.e. (0, 0, ..., 0, 1))

%%%%%%%%%%%%%%%%%%% Initializing Constants: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clf
clear all
set(gcf, 'Position', [100, 800, 800, 800])
set(0, 'defaulttextinterpreter', 'latex')

% Mode Input
m = input(['Enter a mode (1 = debug, 2 = single comp, ' ...
            '3 = init. cont, 4 = tang. cont, 5 = w+ vs. s): ']);

% Options for turning around or slowing down (turn "on" or "off")
        
slow_down = "off";  % turn on for slowdown near 0 discriminant
discrm_tol = .05;     % tolerance for determinant (small, > 0)
ds_turnt = 0;       % whether we have slowed down our search arclength step
ds_slowdown = 0.2;  % percentage of ds after slowdown

turn_around = "off"; % turn on for turnaround after dir_switch % of dir_num
dir_switch = 0.05;   % switch direction after what percentage
turnt = 0;          % whether direction has switched or not


% !!! FILENAME TO SAVE DATA TO !!!
filename = 's_theta_81017_h2_06.txt';
        
% Step Size Constants
L = 50;             % conservative placeholder
dx = 0.05;           % change in space
ds = .5;            % change in direction
dir_num = 20;      % number of direction searches
x = (-L: dx: L)';   % spacing of dx up to L, transpose
N_tot = size(x');
N = N_tot(2);       % how many time steps

ts_data = (1 : 0.05 : 4);

% Initial theta and s parameter values

% Some initial values for dx / L / theta / s that work : 
% dx = 0.2, L = 25, theta = 1.7, s = 0.5
% dx = 0.05, L = 25, theta = 1.8, s = 0.5

theta = 1.8;
s = 0.5;


% Set up Derivative and Averaging Matrix (Finite Differences)

e = ones(N, 1);
D = spdiags([-e e], 0: 1, N-1, N) / dx; % first order derivative, upwind

M = spdiags([e e], 0: 1, N-1, N) / 2;   % averaging matrix

% Phase Condition

cnst = 5;   % steepness

b_amp = theta/(2*s) + sqrt(theta^2/(4*s^2) - 1/s); % eigenval (CHANGE??)

% b_old = -b_amp ./ (1+exp(cnst*x));
b_old = (b_amp / 2) * tanh(cnst*x) + (b_amp / 2);
% b_oldx =  b_amp*cnst * exp(cnst*x) ./ (exp(cnst*x) + 1).^2;
b_oldx = (b_amp / 2) * cnst * (sech(cnst*x)).^2;

u0 = [b_old; b_oldx; s];

% Kinetics (uses cont_df.m)

F = @(u) cont_df_h2(u, N, dx, theta, b_old, b_oldx, M, D);

options = optimset('Jacobian', 'on', 'Display', 'iter', 'TolFun', 1e-8, ...
          'TolX',1e-8,'MaxIter',50,'Algorithm','trust-region-reflective');

% reference if you want to turn Display off (or other)
% options = optimset('Jacobian', 'on','Display','iter','TolFun',1e-8, ...
%          'TolX',1e-8,'MaxIter',50,'Algorithm','trust-region-reflective');

%%%%%%%%%%%%%%%%%%%%%%% DEBUGGING/TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch m % mode
    
    case 1

    ord = input('Enter an order of magnitude (negative): ');
    disp(['Testing Jacobian with order 10^' num2str(ord) '...'])
    du = 10^ord * rand(2*N + 1, 1); % random points, order 'ord'
    [F0, J0] = F(u0);
    [F1, J1] = F(u0 + du);
    err_u = (F1 - F0 - J0*du);
    disp(['Error: ' num2str(norm(err_u))]);
    plot(err_u)


%%%%%%%%%%%%%%%%%%% SINGLE USE OF FSOLVE: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 2

    u_new = fsolve(F, u0, options);

    b = u_new(1: N);
    v = u_new(N+1: 2*N);
    s = u_new(2*N+1);

    plot(x, b, 'b') % space x vs biomass b
    % axis([0 L 0 1])
    xlabel('Space x')
    ylabel('Biomass b')
    title(['Space x vs. Biomass b with s = ' num2str(s) ', theta = ' ...
            num2str(theta)])
    
       
%%%%%%%%%%%%%%%%%%% BABY CONTINUATION IN THETA: %%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    case 3

    u_new = fsolve(F, u0, options);

    b = u_new(1: N);
    v = u_new(N+1: 2*N);
    s = u_new(2*N+1);
    
    % update phase condition
    
    b_old = b;
    b_oldx = v - s*b;
    
    % set up theta arrays for plotting later
    
    theta_start = theta;
    d_theta = 0.05;
    theta_stop = theta - 1.5;
    % set theta_go/stop so direction of continuation d/n matter
    if theta_start > theta_stop
        d_theta = -d_theta;
    end
    theta_go = min(theta_start, theta_stop);
    theta_end = max(theta_start, theta_stop);
    % set range of theta for plotting
    theta_range = theta_start : d_theta : theta_stop;
    theta_vals = NaN(size(theta_range)); % data to be filled with speeds
    theta_ctr = 1; % where in tau range are we

    % initial theta plotting
    
    subplot(2, 1, 2)
    plot(x, b, 'b') % space x vs biomass b
    % axis([0 L 0 1])
    xlabel('Space x')
    ylabel('Biomass b')
    title(['Space x vs. Biomass b with s = ' num2str(s) ', theta = ' ...
            num2str(theta)])
    subplot(2, 1, 2)
    plot(theta_range, theta_vals)
    axis([theta_go theta_end 0 2])
    xlabel('Theta')
    ylabel('Speed s')
    title('Theta vs. Speed s')
 
    drawnow
    
    pause
        
    % vary theta
    
    for new_theta = theta_start : d_theta : theta_stop
    
    % Change F based on new_theta value
    F = @(u) cont_df_h2(u, N, dx, new_theta, b_old, b_oldx, M, D);
    u_new = fsolve(F, u_new, options);
    
    b = u_new(1: N);
    v = u_new(N+1: 2*N);
    s = u_new(2*N+1);
    
    b_old = b;
    b_oldx = v - s*b;
    
    theta_vals(theta_ctr) = s;
    theta_ctr = theta_ctr + 1;
    
    %%%%% Plotting %%%%%
    
    % Plot Space x vs. Biomass b
    subplot(2, 1, 1)
    plot(x, b, 'b') % space x vs biomass b
    % axis([0 L 0 4])
    xlabel('Space x')
    ylabel('Biomass b')
    title(['Space x vs. Biomass b with theta = ' num2str(new_theta) ...
            ', s = ' num2str(s)])
    % Plot theta vs. speed s
    subplot(2, 1, 2)
    plot(theta_range, theta_vals)
    axis([theta_go theta_end 0 2])
    xlabel('Theta')
    ylabel('Speed s')
    title('Theta vs. Speed s')
 
    drawnow
   
    end
  
%%%%%%%%%%%%%%%%%%%%%% TANGENT ARCLENGTH CONTINUATION %%%%%%%%%%%%%%%%%%%%%

    case 4
        
    % number of direction iterations
    theta_data = [];
    s_data = []; % want to plot theta vs. s
    
        
    u_new = fsolve(F, u0, options);

    b = u_new(1: N);
    v = u_new(N+1: 2*N);
    s = u_new(2*N+1);
    
    z_new = [u_new; theta]; 
    
    % Set up for initial cont_dz 
    
    z_old = z_new; % initial z_old?
    dz_old = zeros(1, 2*N+2);
    dz_old(2*N+2) = 1; % initial direction in theta
    
    F = @(z) cont_dz_h2(z, z_old, N, dx, ds, dz_old, b_old, b_oldx, M, D);
    
    % run initial cont_dz fsolve
    [z_new, fval, exitflag, output, jacobian] = fsolve(F, z_old, options);
    
    b = z_new(1: N);
    v = z_new(N+1: 2*N);
    s = z_new(2*N+1);
    theta = z_new(2*N+2);
    theta_data = [theta];
    s_data = [s];
    
    discrm = theta^2 / (4*s^2) - 1 / s;
    discrm_data = [discrm];
    
    z_old = z_new;
    b_old = b;
    b_oldx = v-s*b;
    
    % Update dz 
 
    Z = @(dz_var) dz_solver(dz_var, dz_old, jacobian(1:2*N+1, :), N);

    dz_new = fsolve(Z, dz_old);

    dz_new = dz_new / norm(dz_new) * ds; % normalize
    
    subplot(2, 1, 1)
    plot(x, b, 'g') % space x vs biomass b
    % axis([0 L 0 1])
    xlabel('Space x')
    ylabel('Biomass b')
    title(['Space x vs. Biomass b with s = ' num2str(s) ', theta = ' ...
            num2str(theta)])
    
    subplot(2, 1, 2)
    plot(s_data, theta_data, 'g')
    xlabel('Front speed s')
    ylabel('Theta')
    title(['Front speed s vs. Theta with theta = ' num2str(theta) ...
             ', s = ' num2str(s)])
             
    drawnow
    % pause
    
    plot(ts_data, 0.25*ts_data.^2, 'b')
    drawnow
    hold on
    
       
    cm = colormap(jet(dir_num)); % rainbow coloring
    
    for i = 1: dir_num
    
    F = @(z) cont_dz_h2(z, z_old, N, dx, ds, dz_old, b_old, b_oldx, M, D);
        
    [z_new, fval, exitflag, output, jacobian] = fsolve(F, z_old, options);
    
    b = z_new(1: N);
    v = z_new(N+1: 2*N);
    s = z_new(2*N+1);
    theta = z_new(2*N+2);
    
    discrm = theta^2 / (4*s^2) - 1 / s;
    
    %%%% TEMPORARY %%%%
    
    if discrm < 0.0071
        break
    end
    
    disp(['discrim = ' num2str(discrm) ', theta = ' num2str(theta) ...
            ', s = ' num2str(s)]);
    
    if slow_down == "on"
        % if getting close to 0, smaller arclength size
        if (discrm < discrm_tol) && (ds_turnt == 0)
            disp(["We are in there with discriminant = " num2str(discrm)]);
            ds = ds_slowdown * ds;
            ds_turnt = 1;
        end
    end
    
    % stop if imaginary
    if isreal(theta) == 0 
        disp(['Theta (complex) = ' num2str(theta)])
        disp(['s (complex) = ' num2str(s)])
        break
    end   
    
    theta_data = [theta_data theta];
    s_data = [s_data s];
    discrm_data = [discrm_data discrm];
    
    b_old = b;
    b_oldx = v-s*b;
    
    % Update dz 
 
    Z = @(dz_var) dz_solver(dz_var, dz_old, jacobian(1:2*N+1, :), N);

    dz_new = fsolve(Z, dz_old);
    
    dz_new = dz_new / norm(dz_new) * ds; % normalize, this is our dz
    
    % switch direction after a certain point (initialized above)
    
    if turn_around == "on" && i == floor(dir_switch * dir_num) 
        dz_new = -1*dz_new;
        turnt = 1;
    end
    
    z_old = z_new + dz_new';
    dz_old = dz_new;
    
    b_height = b(floor(3*N/4)); % how high is interface currently
    
    subplot(2, 1, 1)
    plot(x, b, 'Color', cm(i, :)) % space x vs biomass b
    xlabel('Space $$x$$', 'FontWeight', 'bold', 'FontSize', 24)
    ylabel('Biomass $$b$$', 'FontWeight', 'bold', 'FontSize', 24)
    title(['Space $$x$$ vs. Biomass $$b$$ with $$b_+$$ = ' ...
             num2str(b_height)], 'FontWeight', 'bold', 'FontSize', 20)
    hold on
    
    subplot(2, 1, 2)
    plot(theta_data, s_data, 'r')
    ylabel('Front speed $$s$$', 'FontWeight', 'bold', 'FontSize', 24)
    xlabel('$$\theta$$', 'FontWeight', 'bold', 'FontSize', 24)
    title(['$$\theta$$ vs. Front speed $$s$$ with $$\theta$$ = ' num2str(theta) ...
             ', $$s$$ = ' num2str(s) ', discrim = ' num2str(discrm)], ...
             'FontWeight', 'bold', 'FontSize', 20)
    axis([1.4 2.6 0 1.5])
    hold on
    
    drawnow
    
    end
    
    %%%%%%%%%%%%%% ANALYSIS OF theta(s) %%%%%%%%%%%%%%%
    
    % change plot titles
    title('$$\theta$$ vs. Front speed $$s$$', ...
            'FontWeight', 'bold', 'FontSize', 24)
    subplot(2, 1, 1)
    title('Space $$x$$ vs. Biomass $$b$$', ...
          'FontWeight', 'bold', 'FontSize', 24)
    
      
      
    % results = [s_data; theta_data; discrm_data];
    
    %%%% TEMPORARY %%%%
    results = [b; v]
    
    dlmwrite(filename, results);
    
%%%%%%%%%%%%%%%%%%% SCALING AND DATA TRANSFORMATIONS %%%%%%%%%%%%%%%%%%%%%%

    case 5
           
    % Load text files containing data and transform continuation data    
        
    % some good files so far are:
    % 's_vs_theta_80317_t3.txt', 's_vs_theta_80417_t1.txt'
    cont_file = 's_vs_theta_80717_t3.txt';
    data_ct = textread(cont_file, '', 'delimiter', ',', 'emptyvalue', NaN);
    s_res = data_ct(1, :); % s data
    t_res = data_ct(2, :); % theta data
    
    % Implement this from data manually transcribed (TO DO)
    % s_file_1 = 'wp_vs_s_corr.txt';
    % data_s1 = textread(s_file_1, '', 'delimiter', ',', 'emptyvalue', NaN);
    % wp_res_1 = data_s1(1, :); % w_plus data
    % ss_res_1 = data_s1(2, :); % s_sim data
    
    %%%%%%%%%%%%%%%%%%%% TRANSFORMATION AND CONSTANTS %%%%%%%%%%%%%%%%%%%%%
    
    c = 0;  % first value will be c + 1
    c_start = c;
    c_len = 4; % actual length will be c_len + 1
    c_stop = c + c_len;
    d_c = 1;
    
    c_new = c;
   
    % initial test
    w_plus = @(c) t_res ./ sqrt(s_res + c);
    w_res = w_plus(c);
    
    c_range = size(c_start : d_c : c_stop);
    dir_num = c_range(2);
    cm = colormap(jet(dir_num));
    c_small = ['c' 'b' 'g' 'r' 'm'];
    
    for i = 1 : dir_num
        
       c_new = c_new + d_c;
       w_res = w_plus(c_new);
       plot(w_res, s_res, c_small(i))
       % plot(w_res, s_res, 'Color', cm(i, :)) % rainbow colors
       axis([0 8 0 8])
       ylabel('Front speed s')
       xlabel('Water w_+')
       title(['Water w_+ vs. Front speed s'])
       
       hold on
       drawnow
       
    end
    
    % Scatter plots of data from direct simulation (veg_model_80317.m)
    
    % c = 1
    scatter([2 3 4 5 6], [1.46 2.52 3.478 4.38 5.33], c_small(1));
    % c = 2
    scatter([1 2 3 4 5 6], [.584 1.82 2.857 3.81 4.776 5.614], c_small(2)); 
    % c = 3
    scatter([1 2 3 4 5 6], [.87 2.065 3.125 4.098 5 5.926], c_small(3)); 
    % c = 4
    scatter([1 2 3 4 5 6], [1.02 2.26 3.33 4.324 5.246 6.25], c_small(4)); 
    % c = 5
    scatter([1 2 3 4 5 6], [1.15 2.449 3.517 4.571 5.517 6.41], c_small(5));  
    
    otherwise
        
    warning(['Unexpected mode: Type 1 for debug, ' ...
              '2 for initial newton, 3 for baby continuation, ' ...
              '4 for tangent continuation, 5 for data transformations']);
        
end