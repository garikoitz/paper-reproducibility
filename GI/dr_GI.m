function [GI,meanC,isOverlap,minUpper,maxLower,overlapRange]=...
                                 dr_GI(measVals, upperVals, lowerVals, varargin)
%
% 
%
% Syntax:
% 
%
% Description:
% 
%
% Inputs: (required)
% 
%
% Optional key/val pairs:
% 
%
% Examples in the source code
%
% GLU Vistalab, 2018
% v02: 2019-02: cleaning the code
%
% See also:  


%% Parse inputs
p = inputParser;
addRequired(p, 'measVals' , @isnumeric);
addRequired(p, 'upperVals', @isnumeric);
addRequired(p, 'lowerVals', @isnumeric);
addOptional(p, 'LowerRange', 0, @isnumeric);
addOptional(p, 'UpperRange', 1, @isnumeric);

parse(p,measVals,upperVals,lowerVals,varargin{:});
LowerRange      = p.Results.LowerRange;
UpperRange      = p.Results.UpperRange;

% Do the thing
minUpper         = min(upperVals);
maxLower         = max(lowerVals);
maxUpper         = max(upperVals);
minLower         = min(lowerVals);
expRange         = maxUpper - minLower;
overlapRange     = minUpper - maxLower;
GI               = expRange - overlapRange;
if GI < 0; error('expRange - verlapRange cannot be negative'); end

% Divide by the experimental range to make GI unit free
measRange        = UpperRange - LowerRange;
if measRange < 0; error('Check range of the measurement, it cannot be negative'); end
GI = GI/measRange; 

% Correct for sqrt of number of experiments E
E  = length(upperVals);
GI = GI/sqrt(E);

% Use constant for sensitivity
% See Figure S1. Simulations of GI sensitivity changes to c
    c  = 5;    
GI = c * GI;    
    
% Transformation to make GI between 0 and 1;
GI = exp(-GI);  

% Boolean that expreses if there is an overlap or not
isOverlap = false;
if overlapRange > 0; isOverlap = true; end

% Calculate the mean of all the CI-s, the mean of all the experiments, i.e. the
% mean of the generalization experiment.
allC    = upperVals - lowerVals;
meanC   = mean(allC);
