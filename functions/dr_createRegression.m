function  [longVals, BSVals] = dr_createRegression(t, varargin)
%
% Input a matrix with per-subject mean values and it will create a
% correlation matrix for each pair. It will be lower triangular.
%
% Syntax:
%     corrMatTable = dr_createRegression(t, ...)
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

addRequired(p, 't'                ,          @istable);
addOptional(p, 'pMin'             , 0.05   , @isnumeric);
addOptional(p, 'ExtractNumeric'   , true   , @islogical);
addOptional(p, 'tvalue'           , true   , @islogical);
addOptional(p, 'tractsOrder'      , []     , @iscellstr);
addOptional(p, 'newTractNames'    , []     , @iscellstr);
addOptional(p, 'nReps'            , 1000   , @isnumeric);
addOptional(p, 'CIrange'          , 95     , @isnumeric);  % alpha <.05 (two-tailed)
addOptional(p, 'onlyBilateral'    , true   , @islogical);
addOptional(p, 'normResidual'     , false  , @islogical);
addOptional(p, 'Group2Individual' , false  , @islogical);
addOptional(p, 'multiExpRegrCoeff', table(), @istable  );


parse(p,t,varargin{:});

ExtractNumeric   = p.Results.ExtractNumeric;
pMin             = p.Results.pMin;
tvalue           = p.Results.tvalue;
tractsOrder      = p.Results.tractsOrder;
newTractNames    = p.Results.newTractNames;
nReps            = p.Results.nReps;
CIrange          = p.Results.CIrange;
onlyBilateral    = p.Results.onlyBilateral;
normResidual     = p.Results.normResidual;
Group2Individual = p.Results.Group2Individual;
multiExpRegrCoeff= p.Results.multiExpRegrCoeff;

%% Obtain the tract names from the table:
cat = string(unique(t.SliceCats));
if length(cat) ~= 1
    error('We should only have one category here')
end

% Extract the tract pair names
[structPairs, newStructs, LIndex, RIndex] = dr_obtainPairs(t, 'wide');

% Test
if ~isequal(strrep(structPairs(LIndex),'Left','') , strrep(structPairs(RIndex),'Right',''))
    error('There is a problem with the Left and Right structure indexes')
end

% Separate the table into two: 
ttracts = t(:,  ismember(t.Properties.VariableNames, structPairs));
trest   = t(:,  ~ismember(t.Properties.VariableNames, structPairs));

if ~isempty(tractsOrder)
    ttracts = ttracts(:,tractsOrder);
    % Recalculate to recover order
    [structPairs, newStructs, LIndex, RIndex] = dr_obtainPairs(ttracts, 'wide');
end
if ~isempty(newTractNames)
    ttracts.Properties.VariableNames = strrep(newTractNames,' ','');
end
Structure = ttracts.Properties.VariableNames;
t = [trest,ttracts];

if onlyBilateral
  N = length(newStructs);
  % Create the table from the beginning
  longVals = table(repmat(newStructs',[3,1]), ...
                   [repmat("Corr",[N,1]);repmat("Lower",[N,1]);repmat("Upper",[N,1])], ...
                   repmat(nan(N,1),[3,1]), ...
                   'VariableNames',{'CorName', 'Type', char(cat)});
  longVals.CorName = categorical(longVals.CorName);
  BSVals   = table(newStructs', ...
                   repmat({nan(1,height(t))},[N,1]), ...
                   'VariableNames',{'CorName', char(cat)});
  BSVals.CorName = categorical(BSVals.CorName);
  
  % Calculate values and fill the tables             
  for biStr=newStructs
       % Calculate the residuals
       % Depending on the options use the provided line parameters or do the fit
       if Group2Individual
            yreal     = t{:,['Right' biStr{:}]};
            x         = t{:,['Left'  biStr{:}]};
            b         = multiExpRegrCoeff{multiExpRegrCoeff.CorName==biStr{:},'b'};
            slope     = multiExpRegrCoeff{multiExpRegrCoeff.CorName==biStr{:},'slope'};
            ypred     = b + (slope * x);
            residuals = yreal - ypred;
            RMSE      = sqrt(mean(residuals.^2));
            if normResidual
                residuals = residuals/RMSE;
            end
       else
            model   = fitlm(t, ['Right' biStr{:} ' ~ ' 'Left' biStr{:}]);
            if normResidual
                residuals = model.Residuals.Pearson;
            else
                residuals = model.Residuals.Raw;
            end
       end
        %        plotResiduals(model)
        %        boxplot(table2array(model.Residuals))
        %        plotResiduals(model,'probability')
        %        plotResiduals(model,'lagged')
        %        plotResiduals(model,'symmetry')
        %        plotResiduals(model,'fitted')
        

        % This was wrong, it was taking the min and max, and we should put the
        % bootstrapped confidence intervals...
        [Mean,lowerQCI,upperQCI,~,lowerLowerQCI,~,~,~,upperUpperQCI]=dr_bootstrapDistribution(...
                                                                residuals, ...
                                                                'nReps', nReps, ...
                                                                'perc', CIrange, ...
                                                                'grandMean', true);
        
        % Store the values in the table
        longVals.(cat)(longVals.CorName==biStr{:} & longVals.Type=="Corr") = Mean;
        % longVals.(cat)(longVals.CorName==biStr{:} & longVals.Type=="Lower") = lowerLowerQCI;
        % longVals.(cat)(longVals.CorName==biStr{:} & longVals.Type=="Upper") = upperUpperQCI;
        longVals.(cat)(longVals.CorName==biStr{:} & longVals.Type=="Lower") = lowerQCI;
        longVals.(cat)(longVals.CorName==biStr{:} & longVals.Type=="Upper") = upperQCI;
        BSVals.(cat){BSVals.CorName==biStr{:}} = residuals';
  end
    
else
    error('Not implemented yet')
end
end


