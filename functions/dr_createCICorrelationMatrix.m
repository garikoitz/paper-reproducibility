function  [longVals, BSVals] = dr_createCICorrelationMatrix(t, varargin)
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

addRequired(p, 't'             ,       @istable);
addOptional(p, 'pMin'          , 0.05, @isnumeric);
addOptional(p, 'ExtractNumeric', true, @islogical);
addOptional(p, 'tvalue'        , true, @islogical);
addOptional(p, 'tractsOrder'   , []  , @iscellstr);
addOptional(p, 'newTractNames' , []  , @iscellstr);
addOptional(p, 'nReps'         , 1000, @isnumeric);
addOptional(p, 'CIrange'       , 95  , @isnumeric);  % alpha <.05 (two-tailed)
addOptional(p, 'onlyBilateral' , true, @islogical);
parse(p,t,varargin{:});

ExtractNumeric= p.Results.ExtractNumeric;
pMin          = p.Results.pMin;
tvalue        = p.Results.tvalue;
tractsOrder   = p.Results.tractsOrder;
newTractNames = p.Results.newTractNames;
nReps         = p.Results.nReps;
CIrange       = p.Results.CIrange;
onlyBilateral = p.Results.onlyBilateral;

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
    % Recalculate to recover order
    [structPairs, newStructs, LIndex, RIndex] = dr_obtainPairs(t, 'wide');
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

if onlyBilateral
    meanTable = cell2table(newStructs');
    meanTable.Properties.VariableNames = {'CorName'};
    N = height(meanTable);
    meanTable.Type = categorical(repmat("Corr",[N,1]));
    meanTable.(cat) = nan([N,1]);
    for ii=1:N
        corname = meanTable.CorName{ii};
        colcorname = strcat('Left',corname);
        rowcorname = strcat('Right',corname);
        meanTable.(cat)(ii) = corrMatTable{rowcorname,colcorname};
    end
    
    lowerCITable = cell2table(newStructs');
    lowerCITable.Properties.VariableNames = {'CorName'};
    N = height(lowerCITable);
    lowerCITable.Type = categorical(repmat("Lower",[N,1]));
    lowerCITable.(cat) = nan([N,1]);
    for ii=1:N
        corname = lowerCITable.CorName{ii};
        colcorname = strcat('Left',corname);
        rowcorname = strcat('Right',corname);
        lowerCITable.(cat)(ii) = lowerCI{rowcorname,colcorname};
    end
    
    
    upperCITable = cell2table(newStructs');
    upperCITable.Properties.VariableNames = {'CorName'};
    N = height(upperCITable);
    upperCITable.Type = categorical(repmat("Upper",[N,1]));
    upperCITable.(cat) = nan([N,1]);
    for ii=1:N
        corname = upperCITable.CorName{ii};
        colcorname = strcat('Left',corname);
        rowcorname = strcat('Right',corname);
        upperCITable.(cat)(ii) = upperCI{rowcorname,colcorname};
    end    
    
    bstable = cell2table(newStructs');
    bstable.Properties.VariableNames = {'CorName'};
    N = height(bstable);
    bstable.(cat) = repmat({nan(size(bootstrapStat{2,1}{:}))},[N,1]);
    for ii=1:N
        corname = bstable.CorName{ii};
        colcorname = strcat('Left',corname);
        rowcorname = strcat('Right',corname);
        bstable.(cat)(ii) = bootstrapStat{rowcorname,colcorname};
    end    
    
    
    % Join the tables and add then to a resul struct
    longVals  = [meanTable;lowerCITable;upperCITable];
    BSVals    = bstable;
    longVals.CorName = categorical(longVals.CorName);
    BSVals.CorName = categorical(BSVals.CorName);
else
    longVals = struct();
    longVals.corrMatTable = corrMatTable;
    longVals.lowerCI      = lowerCI;
    longVals.upperCI      = upperCI;
    BSVals           = bootstrapStat;
end
end


