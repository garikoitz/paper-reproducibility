function [Mean    ,lowerQCI     ,upperQCI, ...
         lowerMean,lowerLowerQCI,lowerUpperQCI, ...
         upperMean,upperLowerQCI,upperUpperQCI, ...
         bootstrapStatMedian]      = dr_bootstrapDistribution(values, varargin)
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
%
% See also:  


%% Parse inputs
p = inputParser;
addRequired(p, 'values', @isnumeric);
addOptional(p, 'nReps' , 1000, @isnumeric);
addOptional(p, 'perc'  , 90  , @isnumeric);
addOptional(p, 'grandMean'  , false  , @islogical);

parse(p,values,varargin{:});
nReps      = p.Results.nReps;
perc       = p.Results.perc;
grandMean  = p.Results.grandMean;

% We need 0-100 based, not 0-1 based, fix it if it comes in 0-1 range. 
% if perc > 1; perc = perc/100; end
if perc < 1; perc = perc*100; end

% Define the required confidence intervals as two percentiles
twoTailedRange = (100 - perc) / 2;
    

if grandMean
    values = values(:);
end

fmedian    = @(X) median(X,'omitnan');
fpctile    = @(X) prctile(X,[twoTailedRange, 100-twoTailedRange]);
bootstrapStatMedian   = bootstrp(nReps, fmedian, values);
bootstrapStatPrctiles = bootstrp(nReps, fpctile, values); % Duplicates the number of columns
lowerBootstrapStatPrctiles = bootstrapStatPrctiles(:,1:2:size(bootstrapStatPrctiles,2)-1);
upperBootstrapStatPrctiles = bootstrapStatPrctiles(:,2:2:size(bootstrapStatPrctiles,2));


% Mean of the median
Mean = mean(bootstrapStatMedian,'omitnan');
rawCIMean = prctile(bootstrapStatMedian,[twoTailedRange, 100-twoTailedRange]);

% Lower quantile, with CI
lowerQCI    = mean(lowerBootstrapStatPrctiles,'omitnan');
rawCIlowerQCI = prctile(lowerBootstrapStatPrctiles,[twoTailedRange, 100-twoTailedRange]);

% Upper quantile, with CI
upperQCI    = mean(upperBootstrapStatPrctiles,'omitnan');
rawCIupperQCI = prctile(upperBootstrapStatPrctiles,[twoTailedRange, 100-twoTailedRange]);



if grandMean
    lowerMean = rawCIMean(1);
    upperMean = rawCIMean(2);
    lowerLowerQCI = rawCIlowerQCI(1);
    upperLowerQCI = rawCIlowerQCI(2);
    lowerUpperQCI = rawCIupperQCI(1);
    upperUpperQCI = rawCIupperQCI(2);
else
    lowerMean = rawCIMean(1,:);
    upperMean = rawCIMean(2,:);
    lowerLowerQCI = rawCIlowerQCI(1,:);
    upperLowerQCI = rawCIlowerQCI(2,:);
    lowerUpperQCI = rawCIupperQCI(1,:);
    upperUpperQCI = rawCIupperQCI(2,:);
end

