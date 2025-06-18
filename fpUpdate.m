function [x_new,par] = fpUpdate(x,f,par)
% fpUpdate: a simple fixed point solver for x = f(x)
%
% Author: Alex Clymo
% Date: 18/06/2025
%
% Fixed point solver for x = f(x) where x is (N x 1) vector. Problems of
% the form g(x) = 0 can also be solved by transforming to x = g(x) + x and
% so setting f(x) = g(x) + x. 
% fpUpdate performs one step, giving a new guess x_new given current guess
% x and function value f=f(x), using various methods. The idea is that this
% function is useful for several reasons:
%   1) This code can be used in situations where you don't want to
%      explicitly design and build a function f(x). In that case, you are
%      better off sending f(x)-x to fsolve or similar. Instead, this is for
%      when you want to update x as part of a bigger problem that you
%      cannot or do not want to wrap into a function. Key examples are:
%      a) Solving for the equilibrium price sequence following an MIT shock
%      b) Updating certain parameters as part of a calibration routine,
%         where you don't want to send the whole routine to fsolve
%   2) You can quickly switch between updating methods without changing
%      your main code. Basic fixed point updates (x' = f(x)) can be slow or
%      unstable and need dampening. This code allows you to test faster
%      methods that use approximations of the Jacobian with no extra work,
%      and without changing your code at all. 
%   3) The function deals with the creation and storage of past guesses and
%      evaluations in the par structure for methods that use them (e.g.
%      Anderson Acceleration or Broyden's method)
%   4) This is best suited to situations where the computation of f(x)
%      takes a non-trivial amount of time. Then the extra time spent on the
%      updating steps in fpUpdate for more complex methods are more than
%      offset by the savings from converging in fewer iterations. 
%
% Inputs:
%   x          - current x guess (N x 1)
%   f          - current f(x) value (N x 1)
%   par        - structure containing algorithm choice, parameters and past 
%                guesses/evaluations
%   par.method - selects the method:
%                     'fixedPoint' -> basic fixed point update (only works
%                                     for contraction mappings) with
%                                     dampening parameter par.zeta
%                     'anderson'   -> Anderson Acceleration method. See the
%                                     helper function below for required
%                                     par structure and options
%
% Outputs:
%   x_new - updated x guess (N x 1)
%   par   - replicated and possibly updated par structure
%
%
% To do list:
%   1) further testing
%   2) test kkt option for anderson (add conditioning?)
%   3) add more methods: broyden, L-BFGS, ...


switch par.method
    case 'fixedPoint'
        %basic fixed point iteration with dampening:
        % x' = zeta*f(x) + (1-zeta)*x
        zeta = par.zeta;
        x_new = zeta*f + (1-zeta)*x;
    case 'anderson'
        %anderson acceleration
        [x_new,par] = AndersonUpdate(x,f,par);
    otherwise
        error('invalid fixed point method')
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Anderson Acceleration helper function

function [x_new,par] = AndersonUpdate(x,f,par)
% AndersonUpdate performs one Anderson Acceleration step.
% Trying to solve x = f(x)
%
% Inputs:
%   x   - current x guess (N x 1)
%   f   - current f(x) value (N x 1)
%   par - structure containing algorithm pars and past guesses/evaluations
%      .Ma = m+1 - number of past guesses to use (must be >=2)
%                  if Ma < 2 then do simple dampening with zeta0
%      .zeta0    - dampening when not using Anderson because M < Ma
%      .zeta1    - dampening when using Anderson
%                  x_new = zeta*x_new + (1-zeta)*x_hist(:,end)
%      .x_hist   - [N x ?] matrix of past guesses (initialise as  = [])
%      .f_hist   - [N x ?] matrix of f(x) evaluations at those guesses (initialise as  = [])

%
% Output:
%   x_new - [N x 1] updated x guess
%   par   - par structure with updated x_hist, f_hist, and alpha


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up

%unpack structure
Ma = par.Ma;
zeta0 = par.zeta0;
zeta1 = par.zeta1;
m = Ma - 1; %m parameter: past guesses index i=0,1,...,m

% extract x and f histories, or create empty [] if not yet created
if isfield(par,'x_hist') && isfield(par,'f_hist')
    x_hist = par.x_hist;
    f_hist = par.f_hist;
else
    x_hist = [];
    f_hist = [];
end

if Ma < 2, error('Must have Ma >= 2'), end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update: construct x_new

%add current x and f to history
x_hist = [x_hist,x];
f_hist = [f_hist,f];

%extract N and M from x_hist
[N, M] = size(x_hist);
%extract N and M to check x and f have same size
[N2, M2] = size(f_hist);
if (N2 ~= N) || (M2 ~= M), error('x and f history not coherent'), end

if M < Ma
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If less data than anderson columns, just take last f(x) as update
    % with dampening
    x_new = zeta0*f_hist(:,end) + (1-zeta0)*x_hist(:,end);
    alpha = 1; %alpha is just one
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Otherwise do Anderson Acceleration update

    %shrink history to just last Ma values
    x_hist = x_hist(:,end-Ma+1:end);
    f_hist = f_hist(:,end-Ma+1:end);

    % Compute residuals: r_i = f(x_i) - x_i
    R = f_hist - x_hist;  % [N x (m+1)]

    %choose solver method: 'lsqlin' if have optimisation toolbox, if not
    %then use kkt
    method = 'lsqlin';
    switch method
        case 'lsqlin'
            % Compute residuals: r_i = f(x_i) - x_i
            R = f_hist - x_hist;  % [N x (m+1)]
            % lsqlin: min ||R * alpha||^2 s.t. sum(alpha) = 1
            opts = optimoptions('lsqlin','Display','none');
            Aeq = ones(1,m+1);
            beq = 1;
            alpha = lsqlin(R,zeros(N,1),[],[],Aeq,beq,[],[],[],opts);
        case 'kkt'
            % Center residuals by subtracting the last one
            % (optional, but helps numerical stability)
            r_ref = R(:, end);           % r_k
            dR = R(:, 1:end-1) - r_ref;  % [N x m], differences w.r.t. r_k
            % Solve constrained least squares: min ||dR * gamma + r_ref||^2 s.t. sum(gamma) = 1
            % Set up system to solve: [dR, ones(N,1)] * [gamma; lambda] = -r_ref
            % Or use KKT system as:
            G = dR' * dR;              % [m x m]
            c = dR' * r_ref;           % [m x 1]
            % KKT system for constraint sum(alpha) = 1
            A = [G, ones(m,1);
                ones(1,m), 0];
            b = -[c; 1];
            sol = A\b;
            gamma = sol(1:m);
            % Construct alpha vector (length m+1)
            alpha = [gamma; 1 - sum(gamma)];
    end

    % Compute new guess: x_new = sum_i alpha_i * f(x_i)
    x_new = f_hist * alpha;
    %add dampening
    x_new = zeta1*x_new + (1-zeta1)*x_hist(:,end);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update par structure

par.alpha = alpha;
par.x_hist = x_hist;
par.f_hist = f_hist;


end