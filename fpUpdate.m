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
%                                     dampening parameter par.zeta which
%                                     can be scalar or (N x 1) vector
%                     'anderson'   -> Anderson Acceleration method. See the
%                                     helper function below for required
%                                     par structure and options
%   par.xmin   - minimum value for x imposed post updated (scalar or Nx1 vector)
%   par.xmax   - maximum value for x imposed post updated (scalar or Nx1 vector)
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


% check x and f are column vectors
if max(size(x,2),size(f,2)) > 1, error('x and f must be column vectors'), end

% check no NaN, imaginary, or Inf values in x and f 
if max(isnan([x;f])), error('NaN values in x or f'), end
if any(abs([x;f])==Inf), error('Inf or -Inf values in x or f'), end
if max(abs(imag([x;f])))>0, error('Imaginary values in x or f'), end

% select method and update x to x_new
switch par.method
    case 'fixedPoint'
        %basic fixed point iteration with dampening:
        % x' = zeta*f(x) + (1-zeta)*x
        zeta = par.zeta;
        %check that zeta is scalar or column vector of length N
        if ~(size(zeta,2) == 1 && (size(zeta,1) == 1 || size(zeta,1) == size(x,1))), error('zeta must be scalar or N x 1 vector'), end
        %update x
        x_new = zeta.*f + (1-zeta).*x;
    case 'anderson'
        %anderson acceleration
        [x_new,par] = AndersonUpdate(x,f,par);
    otherwise
        error('invalid fixed point method')
end

% impose max and min bounds
x_new = max(x_new,par.xmin);
x_new = min(x_new,par.xmax);

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
%      .zeta0    - dampening when not using Anderson because M < Ma (can be
%                  scalar or (N x 1) vector)
%      .zeta1    - dampening when using Anderson (can be scalar or (N x 1) vector)
%                  x_new = zeta*x_new + (1-zeta)*x_hist(:,end)
%      .maxCondR - maximum condition number of residual maxtrix R before impose
%                  regularisation. if cond(R) > maxCondR then
%                  regularisation term lambda^2 is added to regression so
%                  that the condition number is fixed at maxCondR
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
maxCondR = par.maxCondR;
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
%check that zeta0 and zeta1 are scalar or column vectors of length N
if ~(size(zeta0,2) == 1 && (size(zeta0,1) == 1 || size(zeta0,1) == N)), error('zeta0 must be scalar or N x 1 vector'), end
if ~(size(zeta1,2) == 1 && (size(zeta1,1) == 1 || size(zeta1,1) == N)), error('zeta1 must be scalar or N x 1 vector'), end

if M < Ma
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If less data than anderson columns, just take last f(x) as update
    % with dampening
    x_new = zeta0.*f_hist(:,end) + (1-zeta0).*x_hist(:,end);
    alpha = 1; %alpha is just one
    lambda = 0;
    condR = [];
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Otherwise do Anderson Acceleration update

    %shrink history to just last Ma values
    x_hist = x_hist(:,end-Ma+1:end);
    f_hist = f_hist(:,end-Ma+1:end);

    % Compute residuals: r_i = f(x_i) - x_i
    R = f_hist - x_hist;  % [N x (m+1)]


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up conditioning of R'R matrix


    % Singular values of R from largest to smallest
    s = svd(R);
    % Condition number (ratio svals(1) / svals(end) alternatively, condR = cond(R))
    condR = s(1)/s(end);
    
    % Regularisation: Use parameter lambda to improve conditioning of the least
    % squares problem. In the standard least squares problem, lambda changes
    % the solution from alpha = (R'R)\(R'b) to alpha = (R'R + lambda^2*I)\(R'b)
    % to increasing the conditioning of (R'R + lambda^2*I). We impose a maximum
    % condition number on (R'R + lambda^2*I), where the condition number of R'R
    % is equal to cond(R)^2.
    %condRR = condR^2;
    maxCondRR = maxCondR^2;
    % choose lambda to st cond = maxCondRR, or lambda = 0 if condR < max
    lambda = sqrt( max(s(1)^2 - maxCondRR * s(end)^2,0)/(maxCondRR-1) );
    
    % CHECK: below is condition number of R'R and regularised matrix
    %cond(R'*R)
    %cond(R'*R + lambda^2*eye(size(f_hist,2)))


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run lsqlin to get alpha vector

    if maxCondR > 1 % standard case with some max cond (could be infinite)
        % lsqlin: min ||R * alpha||^2 + lambda^2 ||alpha||^2 s.t. sum(alpha) = 1
        opts = optimoptions('lsqlin','Display','none');
        % constraint sum(alpha) = 1
        Aeq = ones(1,m+1);
        beq = 1;
        % Augment R and RHS for Tikhonov regularization
        R_aug = [R; lambda * eye(m+1)];
        rhs_aug = [zeros(size(R,1),1); zeros(m+1,1)];
        % run lsqlin to get alpha
        alpha = lsqlin(R_aug,rhs_aug,[],[],Aeq,beq,[],[],[],opts);
        % OLD: version without regularisation
        %alpha = lsqlin(R,zeros(N,1),[],[],Aeq,beq,[],[],[],opts);
    else %if maxConR < 1 then tells code to turn of Anderson and just do equal weights alpha = 1/Ma
        alpha = ones(Ma,1)/Ma;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update to x_new with dampening

    % Compute new guess: x_new = sum_i alpha_i * f(x_i)
    x_new = f_hist * alpha;
    %add dampening
    x_new = zeta1.*x_new + (1-zeta1).*x_hist(:,end);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update par structure

par.alpha = alpha;
par.lambda = lambda;
par.x_hist = x_hist;
par.f_hist = f_hist;
par.condR = condR;


end