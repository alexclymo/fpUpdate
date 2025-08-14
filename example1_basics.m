% example1_basics: Script to test fpUdate function for solving x = f(x)
%
% Author: Alex Clymo
% Repository: github->alexclymo->fpUpdate
%
% This code solves the fixed point for a test function f(x) defined in the
% code. Set method to anderson or fixedPoint and see that anderson finds
% the solution in fewer iterations. Code compares solution to that found by
% fsolve to check the fpUpdate method found the same solution. 


clear all
close all
clc

%choose test function from list
testfun = 1

% choose method and set up par structure with default parameters
method = 'broyden'


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define test function and initial guess

%choose test function from list
switch testfun
    case 1 % easy diagonal 2D function
        N = 2;
        f = @(x) [2+cos(x(1)); 1+sin(x(2))];
    case 2 % 2D function with only off-diagonal jacobian (hard, good for testing)
        N = 2;
        f = @(x) [cos(x(2)); sin(x(1))];
    case 3 % another hard 2D function with only off-diagonal jacobian
        N = 2;
        f = @(x) [tanh(5 * x(2)); -tanh(5 * x(1))];
    case 4 % a random N dimension function I wrote
        % x is N x 1 vector
        N = 1000;
        % define f(x) function
        A = magic(N);
        b = (1:N)';
        f = @(x) 1000 + (1:N)' + 0.5*x + 2*x.^0.5 + 1000./max(1,sum((x.^2)'.*x,2));
    otherwise
        error('invalid test function choice')
end

%initial guess
x0 = 0.1*ones(N, 1);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose fixed point solver scheme and set up options

% choose method and set up par structure with default parameters
par = fpSetup(method);

% OPTIONAL: manually make changes to solver parameters by editing par
switch method
    case 'fixedPoint'
        %par.zeta = 0.5;
    case 'anderson'
        %par.Ma = 5; %number of last guesses to use in Anderson
        %par.zeta0 = 0.01; %dampening during pre-Anderson phase
        %par.zeta = 0.5; %dampening during Anderson phase
        %par.maxCondR = 10; %maximum condition number of R before impose regularisation (ridge regression)
    case 'broyden'
        if testfun == 4
            par.H0scale = -1e-3; %flip initial update direction away from 0
        end
        %par.zeta = 0.5;
        %par.H0scale = 0.1; %scale of initial H guess: H0 = par.H0scale * eye(n). Smaller -> smaller initial steps. Can be <0 to flip initial direction
        %par.tolDenom = 1e-10; %if denominator in update below this number then Jacobian update is skipped
        %par.maxCondH = 1e10; %maximum condition number of H before reset Jacobian to identity
    case 'jacob_diag'
        %par.zeta = 0.5;
        %par.zeta0 = 1e-5; %fixed dampening for first iteration
        %par.dmin = 0.0001; %minimum derivative (keeps updates stable)
end

% impose bounds on x if needed, depending on test function choice
switch testfun
    case 4 % a random N dimension function I wrote
        par.xmin = 0;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop to solve fixed point: A while loop repeatedly evalates f(x) and
% constructs a new guess using the fpUdate function, stopping once the
% difference between x and f(x) is sufficiently small

% initialise x at initial guess
x = x0;

% start timer
tic

% start loop over x updates
diff = 1;
tol = 1e-5;
iter = 1;
maxiter = 100000;
while diff > tol && iter <= maxiter

    % evaluate function at new guess
    y = f(x);

    % Perform one x update to get x_new
    [x_new,par] = fpUpdate(x,y,par);

    % compute percentage error (add small constant to denom to handle x=0)
    diff = max( abs( (y-x)./(1e-3+abs(x)) ) )

    % update x by setting x = x_new (only if not converged yet)
    if diff > tol
        x = x_new;
    end

    % increment iteration counter
    iter = iter + 1;

end

% error if failed to converge
if diff > tol
    error('fpUpdate algorithm failed to converge')
end

%number of iterations to reach convergence
iter

%time in fpUpdate algorithm
time_in_fpUpdate = toc


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot diff and dampening across iterations

xPlot = (1:par.iterData.iter);

MM = 1; NN = 2;

subplot(MM,NN,1)
plot(xPlot,par.iterData.rmseList)
xlabel('iteration')
title('RMSE')
grid on

subplot(MM,NN,2)
plot(xPlot,par.iterData.zetaList)
xlabel('iteration')
title('Dampening')
grid on


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare to fsolve solution for validation

%save fpUpdate solution for comparison
x_fpUpdate = x;

tic 
%attempt to solve problem using fsolve: g(x) = f(x) - x = 0
opts = optimoptions('fsolve','Display','off');
%opts = optimoptions('fsolve','Display','off','FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,'StepTolerance',1e-10);
[x_fsolve,fval,exf] = fsolve(@(x) f(x) - x,x0,opts);

%time in fsolve algorithm
time_in_fsolve = toc

%compare fsolve and fpUpdate solutions, if fsolve converged
if exf < 1
    warning('fsolve failed to find solution: not comparing to fsolve')
else
    %compare fpUpdate and fsolve solutions
    diff_vs_fsolve = max(abs((x_fpUpdate - x_fsolve)./(1e-3 + abs(x_fpUpdate))))

    if diff_vs_fsolve > 1e-3
        warning('looks like fpUpdate did not find same solution as fsolve')
    end
end


