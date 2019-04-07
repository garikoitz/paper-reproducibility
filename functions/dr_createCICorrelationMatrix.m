function  [corrMatTable, lowerCI, upperCI, bootstrapStat] = dr_createCICorrelationMatrix(t, varargin)
%
% Input a matrix with per-subject mean values and it will create a
% correlation matrix for each pair. It will be lower triangular.
%
% Syntax:
%     corrMatTable = dr_createCorrelationMatrix(t, ...)
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

addRequired(p, 't', @istable);
addOptional(p, 'pMin', 0.05, @isnumeric);
addOptional(p, 'ExtractNumeric', true, @islogical);
addOptional(p, 'tvalue', true, @islogical);
addOptional(p, 'tractsOrder', [], @iscellstr);
addOptional(p, 'newTractNames', [], @iscellstr);
addOptional(p, 'useBootstrap', true, @islogical);
addOptional(p, 'nReps', 1000, @isnumeric);
addOptional(p, 'CIrange', 95, @isnumeric);  % alpha <.05 (two-tailed)
parse(p,t,varargin{:});

ExtractNumeric= p.Results.ExtractNumeric;
pMin          = p.Results.pMin;
tvalue        = p.Results.tvalue;
tractsOrder   = p.Results.tractsOrder;
newTractNames = p.Results.newTractNames;
useBootstrap  = p.Results.useBootstrap;
nReps         = p.Results.nReps;
CIrange       = p.Results.CIrange;

%% Obtain the tract names from the table:
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
bootstrapStat = [];
lowerCI         = [];
upperCI         = [];
if useBootstrap
    % Create the functions for the correlation with options and the CI
    fcorr   = @(X) corr(X,'rows', 'pairwise');
    % Bootstrap the results
    bootstrapStat  = bootstrp(nReps, fcorr, t{:,:});
    % Define the required confidence intervals
    twoTailedRange = (100 - CIrange)/2;
    rawCI = prctile(bootstrapStat, [twoTailedRange, 100-twoTailedRange]);
    
    % Calculate the correlations
    [corrXmeans, Xpval] = corr(t{:,:},'rows', 'pairwise');
    
    % Make 0 the values that include zero in the CI
    % makeTheseZero = logical(reshape(abs((rawCI(1,:)>0) - (rawCI(2,:)>0)), size(corrXmeans)));
    lowerCI = reshape(rawCI(1,:), size(corrXmeans));
    upperCI = reshape(rawCI(2,:), size(corrXmeans));
    
    % Make bootstrapStat a cell array to pass all vectors between functions
    bootstrapStat = reshape(mat2cell(bootstrapStat', ...
                                     [ones(1,size(bootstrapStat,2))])', ...
                            size(corrXmeans));
    
    % Diagonalize the matrix
    corrXmeans(logical(triu(ones(size(corrXmeans)),0))) = 0;
    lowerCI(logical(triu(ones(size(lowerCI)),0))) = 0;
    upperCI(logical(triu(ones(size(upperCI)),0))) = 0;
    bootstrapStat(logical(triu(ones(size(bootstrapStat)),0))) = {0};
    
    % Create a table again for visualization.
    corrMatTable = array2table(corrXmeans);
    corrMatTable.Properties.VariableNames = Structure;
    corrMatTable.Properties.RowNames      = Structure;

    lowerCI = array2table(lowerCI);
    lowerCI.Properties.VariableNames = Structure;
    lowerCI.Properties.RowNames      = Structure;

    upperCI = array2table(upperCI);
    upperCI.Properties.VariableNames = Structure;
    upperCI.Properties.RowNames      = Structure;    
    
    bootstrapStat = array2table(bootstrapStat);
    bootstrapStat.Properties.VariableNames = Structure;
    bootstrapStat.Properties.RowNames      = Structure;    
else
    % Create the correlation matrix
    % [corrXmeans, Xpval]     = corr(t{:,:},'rows', 'pairwise');
    Alpha = (100-CIrange)/100;
    [corrXmeans,Xpval,lowerCI,upperCI] = corrcoef(t{:,:},'rows', 'pairwise','Alpha',Alpha); 
    % Make 0 all correlation that are not > pMin
    % corrXmeans(Xpval >= pMin)      = 0;

    % Diagonalize the matrix
    corrXmeans(logical(triu(ones(size(corrXmeans)),0))) = 0;
    % Create a table again for visualization.
    corrMatTable = array2table(corrXmeans);
    corrMatTable.Properties.VariableNames = Structure;
    corrMatTable.Properties.RowNames      = Structure;
end

end


