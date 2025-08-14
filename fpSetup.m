function par = fpSetup(method)
% fpSolve: initialise options structure for fpUpdate function
%
% Author: Alex Clymo
% Repository: github->alexclymo->fpUpdate
%
% Optional function to set up par structure for fpUpdate code. Can also do
% so manually, but this code will correctly set up the structure for you
% with the default parameter values. These can then be changed later by
% changing the fields of the par structure. 
%
% Input:
%   method          - solver method to use:
%                     'fixedPoint' -> basic fixed point update (only works
%                                     for contraction mappings) with
%                                     dampening parameter par.zeta
%                     'anderson'   -> Anderson Acceleration method. See the
%                                     fpUpdate function for required par
%                                     structure and options
% Output:
%   par - structure with default options for this method


% Set up par structure with default parameter choices
switch method
    case 'fixedPoint'
        par.method = 'fixedPoint';
    case 'anderson'
        par.method = 'anderson';
        par.Ma = 5; %number of last guesses to use in Anderson
        par.zeta0 = 0.01; %fixed dampening during pre-Anderson phase
        par.maxCondR = 10; %maximum condition number of R before impose regularisation (ridge regression)
    case 'broyden'
        par.method = 'broyden';
        par.H0scale = 0.1; %scale of initial H guess: H0 = par.H0scale * eye(n). Smaller -> smaller initial steps. Can be <0 to flip initial direction
        par.tolDenom = 1e-10; %if denominator in update below this number then Jacobian update is skipped
        par.maxCondH = 1e10; %maximum condition number of H before reset Jacobian to identity
    case 'jacob_diag'
        par.method = 'jacob_diag';
        par.zeta0 = 1e-5; %fixed dampening for first iteration
        par.dmin = 0.0001; %minimum derivative (keeps updates stable)
    case 'jacob_rank1'
        par.method = 'jacob_rank1';
    otherwise
        error('invalid fixed point method')
end

% adaptive dampening controls
par.zeta = 0.5; %initial dampening (zeta = 1 is no dampening)
par.adaptiveDampening = 'on'; %if off, then fixed dampening of zeta
par.adSettings.shrinkFactor = 0.5; %must be <= 1 (zeta -> shrinkFactor*zeta)
par.adSettings.growFactor = 1 + 0.8*(1./par.adSettings.shrinkFactor-1); %must be >= 1 (zeta -> growFactor*zeta) and good to be < 1/shrinkFactor
par.adSettings.zetaMin = 1e-3;
par.adSettings.zetaMax = 1;

% no bounds on x by default
par.xmin = -Inf;
par.xmax = Inf;

% initialise structure for carrying past data used in calculations
par.iterData.iter = 0; %initialise counter for current iteration number
par.iterData.rmseList = []; %vector containing root mean squared error at each iteration
par.iterData.zetaList = []; %vector containing zeta at each iteration

par.verbose = 1; %display information / warnings on or off