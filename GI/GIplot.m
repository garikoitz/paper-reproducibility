function GIplot(unsMeans,includeExp,varargin)
%CREATEPROFILE Summary of this function goes here
%   Detailed explanation goes here
%
% 
% Syntax:
%     createBars(dt)
%
% Description:
%  Input a table with profiles and it will plot them
%
% Inputs: (required)
%  dt: datatable
% 
% Optionals: 
% fnameRoot : string
% saveItHere: string
% savePng   : boolean
% saveSvg   : boolean
% WahlOrder : boolean
% cmapname  : string 
% nRep      : float
% CIrange   : array of floats [1:100], can be just 90, for example
% useDistribution = boolean
% includeExp:  cell string with experiment names to include in the calculation
%
% Examples:
%{
%}
% 
% GLU Vistalab, 2018


%% PARSE INPUTS
p = inputParser;

addRequired(p, 'unsMeans');
addRequired(p, 'includeExp');

addOptional(p, 'fnameRoot'      , "changeThisName" , @isstring);
addOptional(p, 'saveItHere'     , "~/tmp"          , @isstring);
addOptional(p, 'savePng'        , false            , @islogical);
addOptional(p, 'saveSvg'        , false            , @islogical);
addOptional(p, 'WahlOrder'      , false            , @islogical);
addOptional(p, 'cmapname'       , "copper"         , @isstring);
addOptional(p, 'nRep'           , 500              , @isfloat);
addOptional(p, 'CIrange'        , 90               , @isfloat);
addOptional(p, 'useDistribution', true             , @islogical);
addOptional(p, 'ResultType'     , 'corr'           , @ischar);

parse(p,unsMeans,includeExp,varargin{:});

fnameRoot       = p.Results.fnameRoot;
saveItHere      = p.Results.saveItHere;
savePng         = p.Results.savePng;
saveSvg         = p.Results.saveSvg;
WahlOrder       = p.Results.WahlOrder;
cmapname        = p.Results.cmapname;
nRep            = p.Results.nRep;
CIrangeOrVals   = p.Results.CIrange;
useDistribution = p.Results.useDistribution;
ResultType      = p.Results.ResultType;

%% PREPARE THE DATA

if ~istable(unsMeans)
    error("uneMeans needs to be a table")
end

if ~iscellstr(includeExp) | ~(length(includeExp)>0)
    error("includeExp need to be a cellstr with more than one experiment")
end

tractsOrder = { 'LeftCingulumCingulate'  , 'RightCingulumCingulate'  , ...
                'LeftArcuate'            , 'RightArcuate'            , ...
                'LeftIFOF'               , 'RightIFOF'               , ...
                'LeftILF'                , 'RightILF'                , ...
                'LeftUncinate'           , 'RightUncinate'           , ...
                'LeftCorticospinal'      , 'RightCorticospinal'      };
WahlTractNames = {  'CBleft'  , 'CBright'  , ...
                    'AFleft'  , 'AFright'  , ...
                    'IFOleft' , 'IFOright' , ...
                    'ILFleft' , 'ILFright' , ...
                    'UFleft'  , 'UFright'  , ...
                    'CSTleft' , 'CSTright' };

allValues = struct();
for CIrange=CIrangeOrVals
    CI = strcat("CI",num2str(CIrange));
    allValues.(CI) = struct();
% Per each CI in, OBTAIN mean value, lower and upper (and pass all the values)
    for ns = 1:length(includeExp)
        cat = string(includeExp{ns});
        [longVals, BSVals]   = dr_meanBoots(unsMeans(unsMeans.SliceCats==cat,:), ... 
                               'tractsOrder',tractsOrder,'newTractNames',WahlTractNames, ...
                               'useDistribution',useDistribution,'nRep',nRep,'CIrange',CIrange);
        if ns==1
            allValues.(CI).longVals = longVals;
            allValues.(CI).BSVals   = BSVals;
        else % for the rest of the cases we will join the table using CorName and Type
            allValues.(CI).longVals = join(allValues.(CI).longVals,longVals);
            allValues.(CI).BSVals   = join(allValues.(CI).BSVals,BSVals);
        end
    end
end


%% CREATE FIGURE AND PLOT
bigfig = figure('Name',fnameRoot, ...
                'NumberTitle','off', ...
                'visible',   'on', ...
                'color','w', ...
                'WindowStyle','normal', ...
                'Units','inches', ...
                'OuterPosition',[0 0 10 14]);
% Instead of the bar plot use the generalization index
referenceColumn = 'NoReference';
versus          = 'NoRetest';
plotIt          = true;
plotIndex       = false;
showLegend = false; showXnames = false; refLine=true;
GIcoloringMethod = 'bars'; % 'circle', 'GIband', 'GRband', 'background','none','bars'

ncol=3; nrow=4;
for nc = 1:length(WahlTractNames)
    long  = allValues.(CI).longVals;
    longB = allValues.(CI).BSVals;
    % Calculate the position
    tn  = string(long.CorName(nc));
    sp = subplot(ncol,nrow,nc);
    % DO the calculation/plots
    [WahlInsideCI,isOverlap,minUpper,maxLower,midCI,GI]=dr_compareCI(long, longB, tn,...
                                           'referenceColumn',referenceColumn,...
                                           'refLine' , refLine,'plotIndex',plotIndex, ...
                                           'plotIt',plotIt,'showLegend',showLegend, ...
                                           'includeExp',includeExp, 'showXnames',showXnames, ...
                                           'GIcoloringMethod', GIcoloringMethod, ...
                                           'ResultType', ResultType);
    % xlabel('FA (TEST)','FontWeight','bold');
    ylabel('FA','FontWeight','bold');
    set(gca,'FontSize',18)
    title(sprintf('%s',tn))
    
    if ismember(nc, [1:4:12])
        ylabel('FA','FontWeight','bold');
    else
        ylabel('');
        set(gca,'YTickLabel',[]);
    end
    
    
end
colormap(cmapname)
h=colorbar;
set(h, 'Position', [.93 0.115 0.01 .83])
suptitle({strrep(fnameRoot,'_','\_'), 'GI plots'})

%% SAVE
if ~exist(saveItHere,'dir')
    mkdir(saveItHere)
end
if saveSvg
    saveas(gcf,fullfile(saveItHere, strcat(fnameRoot,'.svg')),'svg');
end
if savePng
    saveas(gcf,fullfile(saveItHere, strcat(fnameRoot,'.png')),'png');
end
