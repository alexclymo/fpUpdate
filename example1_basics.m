% example1_basics: Script to test fpUdate function for solving x = f(x)
%
% Author: Alex Clymo
% Date: 18/06/2025
%
% This code solves the fixed point for a test function f(x) defined in the
% code. Set method to anderson or fixedPoint and see that anderson finds
% the solution in fewer iterations. Code compares solution to that found by
% fsolve to check the fpUpdate method found the same solution. 


clear all
close all
clc

%choose test function from list
testfun = 3

% choose method and set up par structure with default parameters
method = 'anderson' %'anderson' or 'fixedPoint'


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define test function and initial guess

%choose test function from list
switch testfun
    case 1 % 2D function with only off-diagonal jacobian (hard, good for testing)
        N = 2;
        f = @(x) [cos(x(2)); sin(x(1))];
    case 2 % another hard 2D function, harder for anderson to converge
        N = 2;
        f = @(x) [tanh(5 * x(2)); tanh(5 * x(1))];
    case 3 % a random N dimension function I wrote
        % x is N x 1 vector
        N = 1000;
        % define f(x) function
        A = magic(N);
        b = (1:N)';
        f = @(x) 100 + (1:N)' + 2*x.^0.5 + 1000./max(1,sum((x.^2)'.*x,2));
    otherwise
        error('invalid test function choice')
end

%initial guess
x0 = 0.1*ones(N, 1);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose fixed point solver scheme and set up options

% choose method and set up par structure with default parameters
parFP = fpSetup(method);

% OPTIONAL: manually make changes to solver parameters by editing parFP
switch method
    case 'anderson'
        parFP.Ma = 5; %number of last guesses to use in Anderson
        parFP.zeta0 = 0.01; %dampening during pre-Anderson phase
        parFP.zeta = 0.5; %dampening during Anderson phase
        parFP.maxCondR = 10; %maximum condition number of R before impose regularisation (ridge regression)
    case 'fixedPoint'
        parFP.zeta = 0.5;
    otherwise
        error('invalid fixed point method')
end

% impose bounds on x if needed, depending on test function choice
switch testfun
    case 3 % a random N dimension function I wrote
        parFP.xmin = 0;
        
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
maxiter = 1000;
while diff > tol && iter <= maxiter

    % evaluate function at new guess
    y = f(x);

    % Perform one x update to get x_new
    [x_new,parFP] = fpUpdate(x,y,parFP);

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

xPlot = (1:parFP.iterData.iter);

MM = 1; NN = 2;

subplot(MM,NN,1)
plot(xPlot,parFP.iterData.rmseList)
xlabel('iteration')
title('RMSE')
grid on

subplot(MM,NN,2)
plot(xPlot,parFP.iterData.zetaList)
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


