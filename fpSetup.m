function par = fpSetup(method)
% fpSolve: initialise options structure for fpUpdate function
%
% Author: Alex Clymo
% Date: 18/06/2025
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
        par.zeta = 0.5; %dampening (zeta = 1 is no dampening)
    case 'anderson'
        par.method = 'anderson';
        par.Ma = 5; %number of last guesses to use in Anderson
        par.zeta0 = 0.01; %dampening during pre-Anderson phase
        par.zeta1 = 1; %dampening during Anderson phase
    otherwise
        error('invalid fixed point method')
end