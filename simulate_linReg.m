function [phenotype, X, gTruth] = simulate_linReg(nyVars, nXvars, numObs, varargin)
% Function to simulate phenotype(s) using a linear regression (simple
% additive model) assumption: y = beta * X + e 
%% Mandatory inputs:
% nyVars:       [1 x 1]     number indicating how many phenotypes should
%                           be simulated (v) (default: 1)
%
% nXvars:       [1 x 1]     number indicating how many X variables should
%                           be simualated (p) (default: 5)
% 
% numObs:       [1 x 1]     number indicating how many observations should
%                           be simulated (n) (default: 1000)
% 
%% Optional input(s): should be provided as name-pair values
% intercept:    logical     true/false indicating if the X variables should
%                           include an intercept or not (default: true)
% 
% seed:         [1 x 1]     number used to set the rng seed (default: no
%                           seed is set)
% 
% verify:       logical     true/false indicating if post fit scatter plot
%                           against ground truth should be shown as
%                           verification (default: false)
% 
% beta:         [p x v]     vector or matrix of already known effect sizes;
%                           if this is provided, nXvars has no effect;
%                           additionally, if only a single column is
%                           provided, the same effect sizes are used for
%                           all nyVars
% 
% X:            [n x p]     vector or matrix of X variables to use; if this
%                           is provided, X variables are not synthesied and
%                           nXvars and numObs are not used
% 
%% Outputs:
% phenotype:    [n x v]     synthesised phenotype
% 
% X:            [n x p]     vector or matrix of X variables
% 
% gTruth:       [p x v]     vector of ground truth beta coefficients
% 
%% Examples:
% Example 1:
% Simulate 100 phenotypes with 20 X variables and 10000 observations;
% setting the seed as 123 and requesting post simulation verification
% [phenotype, X, gTruth] = simulate_linReg(100, 20, 10000, 'seed', 123, 'verify', true);
% 
% Example 2:
% Simulate 1 phenotype with 5 X variables and 2000 observations without
% intercept
% [phenotype, X, gTruth] = simulate_linReg(1, 5, 2000, 'intercept', false);
% 
% Example 3:
% Simulate phenotype by explicitly providing beta coefficients
% [phenotype, X, gTruth] = simulate_linReg(1, 2, 2000, 'beta', [0.1; 0.2]);
% 
% Example 4:
% Simulate phenotype by explicitly providing X matrix; this will result in
% two phenotypes with 1000 observations and 5 X variables
% X = randn(1000, 5);
% [phenotype, X, gTruth] = simulate_linReg(2, [], [], 'X', X);
% or equivalently:
% [phenotype, X, gTruth] = simulate_linReg(2, 5, 2000, 'X', X);
% 
%% Notes:
% X variables, ground truth beta coefficients, and noise are drawn from a
% standard normal distributions

%% Check inputs
if ~exist('nyVars', 'var') || isempty(nyVars)
    nyVars = 1;
else
    if nyVars < 1 || ~isscalar(nyVars)
        error('Number of phenotypes should be a scalar greater than or equal to 1');
    end
end

if ~exist('nXvars', 'var') || isempty(nXvars)
    nXvars = 5;
else
    if nXvars < 1 || ~isscalar(nXvars)
        error('Number of X variables should be a scalar greater than or equal to 1');
    end
end

if ~exist('numObs', 'var') || isempty(numObs)
    numObs = 1000;
else
    if numObs < 1 || ~isscalar(numObs)
        error('Number of observations should be a scalar greater than or equal to 1');
    end
end

%% Optional arguments
p = inputParser;
addParameter(p, 'intercept', true,  @islogical);
addParameter(p, 'seed',      [],    @(x) isempty(x) || isscalar(x));
addParameter(p, 'verify',    false, @islogical);
addParameter(p, 'beta',      [],    @isnumeric);
addParameter(p, 'X',         [],    @isnumeric);
parse(p, varargin{:});

intercept = p.Results.intercept;
seed      = p.Results.seed;
verify    = p.Results.verify;
beta      = p.Results.beta;
X         = p.Results.X;

%% Set seed, if required
if ~isempty(seed)
    rng(seed, 'twister');
end

%% Create X variables, if required
if isempty(X)
    if intercept
        X = [ones(numObs, 1), randn(numObs, nXvars-1)];
    else
        X = randn(numObs, nXvars);
    end
else
    [numObs, nXvars] = size(X);
end

%% Generate noise
noise = randn(numObs, nyVars);

%% Generate ground truth effect size(s)
if ~isempty(beta)
    [tmp_nX, tmp_ny] = size(beta);

    % Check if a single column of beta is provided
    if tmp_ny == 1
        gTruth = repmat(beta, 1, nyVars);
    else
        % Make sure dimensions match
        if tmp_nX == nXvars
            if tmp_ny == nyVars
                gTruth = beta;
            else
                error('Mismatch between number of phenotypes and number of columns of beta');
            end
        else
            error('Mismatch between number of X variables and number of rows of beta');
        end
    end
else
    gTruth = randn(nXvars, nyVars);
end

%% Simulate y
phenotype = X * gTruth + noise;

%% Verify
if verify
    tmp_estimates = X \ phenotype;
    scatter(gTruth(:), tmp_estimates(:));
end