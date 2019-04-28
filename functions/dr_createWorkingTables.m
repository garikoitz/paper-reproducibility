function  [dt, unstackedProfiles, unstackedMeans, indivTracts, pairedTracts] = dr_createWorkingTables(dt, varargin)
%dr_createWorkingTables Create a datatable with stacked and unstacked values
%                       for posterior analysis. 
%
% Syntax
%     [dt, unstackedProfiles, unstackedMeans] = dr_createWorkingTables(dt, ...)
%
% Description 
%     Input a table  
%
% Inputs (required)
%      
%
% Inputs (optional)
%
% Example: 
%           dr_createWorkingTables(dt, ...
%                                  'sliceBasedOn',           {'Proj', 'SHELL'}, ...
%                                  'SHELL',                  {'1000', '2000', '3000'}, ...
%                                  'removeTractsContaining', {'Callosum'}, ...
%                                  'useThisTractsOnly',       tracts, ...
%                                  'createTractPairs',        true, ...
%                                  'TRT',                    'notRETEST', ...
%                                  'AGE',                     [18,40], ...
%                                  'GENDER',                 {'male', 'female'})
%
% GLU Vistalab, 2018
%
% See also:  


% Wahl tracts
Wahltracts = [{'Cingulate'}, {'Arcuate'},{'IFOF'}, {'ILF'}, {'Uncinate'}, {'Corticospinal'}];

% Create default values and parse the input 
p = inputParser;
addRequired(p, 'dt', @istable);                           

addOptional(p, 'sliceBasedOn',           {'Proj'},     @iscell);      % sliceBasedOn    = {'Proj', 'SHELL'};
addOptional(p, 'SHELL',                  {'1000'},     @iscell);      % SHELL    = {'1000', '2000', '3000'};
addOptional(p, 'removeTractsContaining', {'Callosum'}, @iscell);      % removeTractsContaining = {'Callosum'};
addOptional(p, 'useThisTractsOnly',       Wahltracts,  @iscell);      % useThisTractsOnly         = Wahltracts;
addOptional(p, 'createTractPairs',        true,        @islogical);   % createTractPairs         = true;
addOptional(p, 'TRT',                    'notRETEST',  @ischar);      % TRT    = 'TEST';
addOptional(p, 'AGE',                     [1,99],      @isnumeric);   % AGE = [1,99];
addOptional(p, 'GENDER',                 {'male', 'female'}, @iscell);% GENDER         = {'male', 'female'};
addOptional(p, 'removeOutliers',         false,        @islogical);   % removeOutliers         = true;


parse(p, dt, varargin{:});

SHELL                  = p.Results.SHELL;
removeTractsContaining = p.Results.removeTractsContaining;
tracts                 = p.Results.useThisTractsOnly;
createTractPairs       = p.Results.createTractPairs;
TRT                    = p.Results.TRT;
AGE                    = p.Results.AGE;
GENDER                 = p.Results.GENDER;
sliceBasedOn           = p.Results.sliceBasedOn;
removeOutliers         = p.Results.removeOutliers;

% Extract important variables to first level
try
    dt.SHELL   = dt.AcquMD.Shell;
catch
    dt.SHELL   = categorical(string(dt.AcquMD.scanbValue));
end

dt.AGE     = dt.SubjectMD.AGE;

try
    dt.GENDER  = dt.SubjectMD.GENDER;
catch
    dt.GENDER  = categorical(dt.SubjectMD.sex);
end
dt.Struct     = categorical(dt.Struct);
dt.TRT        = categorical(dt.TRT);
% Filters
% SHELL
if ~ismember({'SHELL'}, p.UsingDefaults)
    dt        = dt(ismember(dt.SHELL, SHELL), :);
    dt.SHELL  = removecats(dt.SHELL);
end

% removeTracts
if ~ismember({'removeTractsContaining'}, p.UsingDefaults)
    dt          = dt(~contains(string(dt.Struct(:)),removeTractsContaining), :);
    dt.Struct   = removecats(dt.Struct);
end

% selectTracts
if ~ismember({'useThisTractsOnly'}, p.UsingDefaults)
    dt          = dt(contains(string(dt.Struct(:)),tracts),:);
    dt.Struct   = removecats(dt.Struct);
end

% create tract pairs
if ~ismember({'createTractPairs'}, p.UsingDefaults)
    dt          = dt(ismember(dt.Struct, dr_obtainPairs(dt, 'long')), :);
    dt.Struct   = removecats(dt.Struct);
end

% Filter by AGE if the parameter was explicitly passed (never uses defaults)
if ~ismember({'AGE'}, p.UsingDefaults)
	dt = dt(dt.AGE >= AGE(1) & dt.AGE <= AGE(2), :);
end

% GENDER
if ~ismember({'GENDER'}, p.UsingDefaults)
    dt         = dt(ismember(dt.GENDER, GENDER), :);
    dt.GENDER  = removecats(dt.GENDER);
end

% TRT
if ~ismember({'TRT'}, p.UsingDefaults)
    switch TRT
        case {'notTRAIN'}
            dt = dt(dt.TRT ~= 'TRAIN', :);
        case {'notTEST'}
            dt = dt(dt.TRT ~= 'TEST', :);
        case {'notRETEST'}
            dt = dt(dt.TRT ~= 'RETEST', :);
        case {'TRAIN'}
            dt = dt(dt.TRT == 'TRAIN', :);
        case {'TEST'}
            dt = dt(dt.TRT == 'TEST', :);
        case {'RETEST'}
            dt = dt(dt.TRT == 'RETEST', :);
        otherwise
            error('This case does not make sense')
    end
    
    dt.TRT  = removecats(dt.TRT);
end

% Create the new category combining variables (Proj+Shell covers all data we have for now)
if ~ismember({'sliceBasedOn'}, p.UsingDefaults)
    switch length(sliceBasedOn)
        case {1}
            dt.SliceCats = categorical(string(dt.(sliceBasedOn{1})), 'Ordinal', true);
        case {2}
            dt.SliceCats = categorical(strcat(string(dt.(sliceBasedOn{1})), ...
                                              string(dt.(sliceBasedOn{2}))), ...
                                              'Ordinal', true);
        case {3}
            dt.SliceCats = categorical(strcat(string(dt.(sliceBasedOn{1})), ...
                                              string(dt.(sliceBasedOn{2})), ...
                                              string(dt.(sliceBasedOn{3}))), ...
                                              'Ordinal', true);        
        case {4}
            dt.SliceCats = categorical(strcat(string(dt.(sliceBasedOn{1})), ...
                                              string(dt.(sliceBasedOn{2})), ...
                                              string(dt.(sliceBasedOn{3})), ...
                                              string(dt.(sliceBasedOn{4}))), ...
                                              'Ordinal', true);
        otherwise 
            error('Too many variables to slice. Add it to the switch command.')
    end
end

dt    = dr_addRGBcolumn(dt);
% Make the test-retest colors the same
dt = sortrows(dt, 'SubjID');
listOfSubjIDs = dt.SubjID(dt.TRT=="RETEST");

dt((ismember(dt.SubjID,listOfSubjIDs) & dt.TRT=="RETEST"),'SliceCatsRGB') = ...
         dt((ismember(dt.SubjID,listOfSubjIDs) & dt.TRT=="TEST"),'SliceCatsRGB');
     
     
     
% if size(unique(dt.SliceCats),1) ~= size(unique(dt.SliceCatsRGB),1)
%     error('Define more colors in dr_addRGBcolumn(), there are more experiments in project than colors defined')
% end

% Unstack the tract profiles: one column per Structure
unsProf    = unstack(dt (:,{'SliceCats', 'TRT', 'Proj','SubjID','SHELL','Val','Struct','AGE','GENDER'}),'Val','Struct'); 
% Add an RGB color column to separate values in scatterplots or t-sne
unsProf    = dr_addRGBcolumn(unsProf);
unsProf((ismember(unsProf.SubjID,listOfSubjIDs) & unsProf.TRT=="RETEST"),'SliceCatsRGB') = ...
         unsProf((ismember(unsProf.SubjID,listOfSubjIDs) & unsProf.TRT=="TEST"),'SliceCatsRGB');



% Create the same table but with mean profile values
dt.meanVal= mean(dt.Val, 2);

% Unstack it as well
unsMeans   = unstack(dt  (:,{'SliceCats', 'TRT','Proj','SubjID','SHELL','meanVal','Struct','AGE','GENDER'}),'meanVal','Struct'); 

% Check it worked: isequal(mean(unsProf{1,'LeftArcuate'}), unsMeans{1,'LeftArcuate'})
unsMeans   = dr_addRGBcolumn(unsMeans);
unsMeans((ismember(unsMeans.SubjID,listOfSubjIDs) & unsMeans.TRT=="RETEST"),'SliceCatsRGB') = ...
         unsMeans((ismember(unsMeans.SubjID,listOfSubjIDs) & unsMeans.TRT=="TEST"),'SliceCatsRGB');




% Remove Outliers?
if removeOutliers
    unsMeans   = dr_removeOutliers(unsMeans, 3);
end


% Rename for returning
unstackedProfiles           = unsProf;
unstackedMeans              = unsMeans;
[indivTracts, pairedTracts] = dr_obtainPairs(unsProf, 'wide');
      

end