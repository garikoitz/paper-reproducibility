function  [longVals, BSVals] = dr_meanBoots(t, varargin)
%
% We can use the distribution of the values and report the confidence intervals,
% or we can provide the confidence intervals of the means. 
% For example for FA, we will be interested in knowing what is the range of
% values across experiments, and the CI should be for the distributions, not for
% the means. The CI of the means after bootstraping will be tight, but the
% distributions depends on population characteristics, and in the noise of the
% measurement intstrumentation. 
% 
% When it is a distribution, the bootstraStats that we are passing is the actual
% values of the distribution that went in. 
%
%
% Description:
%  Input a table with meanFA (or others) values in the columns and
%  subjects in rows. Create a correlation matrix.
%
% Inputs: (required)
%  Matlab table with variable names per each structure. 
%
% Optional key/val pairs:
%  'pMin':        minimum p to visualizing correlations
%
% Examples in the source code
%
% GLU Vistalab, 2018
%
% See also:  


%% Parse inputs
p = inputParser;

addRequired(p, 't'              ,       @istable);

addOptional(p, 'pMin'           , 0.05, @isnumeric);
addOptional(p, 'ExtractNumeric' , true, @islogical);
addOptional(p, 'tvalue'         , true, @islogical);
addOptional(p, 'tractsOrder'    , []  , @iscellstr);
addOptional(p, 'useDistribution', true, @islogical);
addOptional(p, 'newTractNames'  , []  , @iscellstr);
addOptional(p, 'nReps'          , 1000, @isnumeric);
addOptional(p, 'CIrange'        , 90  , @isnumeric);  % alpha <.05 (two-tailed)
parse(p,t,varargin{:});

ExtractNumeric  = p.Results.ExtractNumeric;
pMin            = p.Results.pMin;
tvalue          = p.Results.tvalue;
tractsOrder     = p.Results.tractsOrder;
useDistribution = p.Results.useDistribution;
newTractNames   = p.Results.newTractNames;
nReps           = p.Results.nReps;
CIrange         = p.Results.CIrange;

%% Obtain the tract names from the table:
cat = string(unique(t.SliceCats));
if length(cat) ~= 1
    error('We should only have one category here')
end
if ExtractNumeric
    [structPairs, newStructs, LIndex, RIndex] = dr_obtainPairs(t, 'wide');
    % Test
    if ~isequal(strrep(structPairs(LIndex),'Left','') , strrep(structPairs(RIndex),'Right',''))
        error('There is a problem with the Left and Right structure indexes')
    end
    % Separate the table into two: 
    t    = t(:,  ismember(t.Properties.VariableNames, structPairs));
end
if ~isempty(tractsOrder)
    t = t(:,tractsOrder);
end
if ~isempty(newTractNames)
    t.Properties.VariableNames = strrep(newTractNames,' ','');
end
Structure = t.Properties.VariableNames;
% Create short names for better visualization.
% corrplot only visualizes first 5 chars for example
% Structure = strrep(Structure, 'Left', 'L');
% Structure = strrep(Structure, 'Right', 'R');
bootstrapStat   = [];
lowerCI         = [];
upperCI         = [];

matrixT = t{:,:};



% Set useDistribution if we want to obtain the mean and the CI of the mean
if useDistribution
   [Mean    ,lowerQCI     ,upperQCI, ...
    lowerMean,lowerLowerQCI,lowerUpperQCI, ...
    upperMean,upperLowerQCI,upperUpperQCI, ...
    bootstrapStat]= dr_bootstrapDistribution(matrixT, ...
                                             'nReps', nReps, ...
                                             'perc', CIrange, ...
                                             'grandMean', false);
    meanCI  = Mean; 
    lowerCI = lowerLowerQCI;
    upperCI = upperUpperQCI;
    meanCI  = Mean; 
    lowerCI = lowerQCI;
    upperCI = upperQCI;
else 
    % Define the required confidence intervals
    twoTailedRange = (100 - CIrange) / 2;
    % Create the functions for the correlation with options and the CI
    fmean   = @(X) mean(X,'omitnan');
    % Bootstrap the results
    bootstrapStat  = bootstrp(nReps, fmean, matrixT);
    
    rawCI = prctile(bootstrapStat, [twoTailedRange, 100-twoTailedRange]);
    lowerCI = rawCI(1,:);
    upperCI = rawCI(2,:);
    % Calculate the mean of the distribution or the bootstraps
    meanCI  = mean(bootstrapStat);
end    
    

    
    % Make bootstrapStat a cell array to pass all vectors between functions
    bootstrapStat = mat2cell(matrixT', [ones(1,size(matrixT,2))])';
    
    
    % Change this name: bootstrapStat

    
    % Create a table again for visualization.
    
    % Calculate the mean of means or of the distributions per each tract
    meanTable = array2table(meanCI');
    meanTable.Properties.VariableNames = cat;
    meanTable.CorName = categorical(Structure');
    N = height(meanTable);
    meanTable.Type = categorical(repmat("Corr",[N,1]));
    meanTable = meanTable(:,[2,3,1]);
    
    % Calculate the lower CI bootstrapping of the means
    % or the lower-lower quantile for the distribution, per each tract
    lowerCI = array2table(lowerCI');
    lowerCI.Properties.VariableNames = cat;
    lowerCI.CorName = categorical(Structure');
    N = height(lowerCI);
    lowerCI.Type = categorical(repmat("Lower",[N,1]));
    lowerCI = lowerCI(:,[2,3,1]);
    
    % Calculate the lower CI bootstrapping of the means
    % or the lower-lower quantile for the distribution, per each tract
    upperCI = array2table(upperCI');
    upperCI.Properties.VariableNames = cat;
    upperCI.CorName = categorical(Structure');
    N = height(upperCI);
    upperCI.Type = categorical(repmat("Upper",[N,1]));
    upperCI = upperCI(:,[2,3,1]);
    
    % Now store all values as arrays inside a table
    bootstrapStat = array2table(bootstrapStat');
    bootstrapStat.Properties.VariableNames = cat;
    bootstrapStat.CorName = categorical(Structure');
    bootstrapStat = bootstrapStat(:,[2,1]);
    
    
    % Join the tables and add then to a resul struct
    longVals  = [meanTable;lowerCI;upperCI];
    BSVals    = bootstrapStat;

end


