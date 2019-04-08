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
% fnameRoot      : string
% saveItHere     : string
% savePng        : boolean
% saveSvg        : boolean
% WahlOrder      : boolean
% cmapname       : string 
% nRep           : float
% CIrange        : array of floats [1:100], can be just 90, for example
% useDistribution: boolean
% includeExp     :  cell string with experiment names to include in the calculation
% nrowcol        :  [nrow, ncol] for plotting
% calcType       : string we need to know if this is a correlation, a distribution...
% onlyBilateral  : boolean
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
addOptional(p, 'nrowcol'        , [4,3]            , @isfloat);
addOptional(p, 'calcType'       , 'distribution'   , @ischar);
addOptional(p, 'onlyBilateral'  , true             , @islogical);
addOptional(p, 'winSizeInch'    , [0,0,10,14]      , @isfloat);
addOptional(p, 'ylab'           , 'nothing'        , @ischar);


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
nrowcol         = p.Results.nrowcol;
calcType        = p.Results.calcType;
onlyBilateral   = p.Results.onlyBilateral;
winSizeInch     = p.Results.winSizeInch;
ylab            = p.Results.ylab;

%% PREPARE THE DATA

if ~istable(unsMeans)
    error("uneMeans needs to be a table")
end

if ~iscellstr(includeExp) || isempty(includeExp)
    error("includeExp need to be a cellstr with more than one experiment")
end

tractsOrder = { 'LeftCingulumCingulate'  , 'RightCingulumCingulate'  , ...
                'LeftArcuate'            , 'RightArcuate'            , ...
                'LeftIFOF'               , 'RightIFOF'               , ...
                'LeftILF'                , 'RightILF'                , ...
                'LeftUncinate'           , 'RightUncinate'           , ...
                'LeftCorticospinal'      , 'RightCorticospinal'      };
% To make tract names shorter or to look the same to Wahl', we can substitute
% them with the followgin names. 
% The functions take this argument to make the change: 
%      'newTractNames',WahlTractNames, ...
% As in the other sections we are using the AFQ names, and now we don' have
% space problems because we are using only the 6 bilateral correlations, there
% is no need of shortening the tract names. 
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
        switch lower(calcType)
        case {'distribution','distr'}             
            [longVals,BSVals] = dr_meanBoots(unsMeans(unsMeans.SliceCats==cat,:), ... 
                           'tractsOrder',tractsOrder, ... 
                           'useDistribution',useDistribution,'nRep',nRep,'CIrange',CIrange);
        case {'correlation','corr'}
            % The following functions it is programmed to return a diagonal matrix
            % Select option 'onlyBilateral',true to return only bilateral tracts
            [longVals,BSVals] = dr_createCICorrelationMatrix(unsMeans(unsMeans.SliceCats==cat,:), ... 
                           'tractsOrder',tractsOrder, ...
                           'nRep',nRep,'CIrange',CIrange, 'onlyBilateral',onlyBilateral);
            
            
        otherwise
            error('caltType %s not implemented yet.', calcType)
        end
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

if onlyBilateral
    bigfig = figure('Name',fnameRoot, ...
                    'NumberTitle','off', ...
                    'visible',   'on', ...
                    'color','w', ...
                    'WindowStyle','normal', ...
                    'Units','inches', ...
                    'OuterPosition',winSizeInch);
    % Instead of the bar plot use the generalization index
    referenceColumn = 'NoReference';
    versus          = 'NoRetest';
    plotIt          = true;
    plotIndex       = false;
    showLegend = false; showXnames = false; refLine=true;
    GIcoloringMethod = 'bars'; % 'circle', 'GIband', 'GRband', 'background','none','bars'

    nrow=nrowcol(1); ncol=nrowcol(2);
    for nc = 1:size(allValues.(CI).BSVals,1)
        long  = allValues.(CI).longVals;
        longB = allValues.(CI).BSVals;
        % Calculate the position
        tn  = string(long.CorName(nc));
        sp = subplot(nrow,ncol,nc);
        % DO the calculation/plots
        dr_compareCI(long, longB, tn,'referenceColumn',referenceColumn,...
                     'refLine' , refLine,'plotIndex',plotIndex, ...
                     'plotIt',plotIt,'showLegend',showLegend, ...
                     'includeExp',includeExp, 'showXnames',showXnames, ...
                     'GIcoloringMethod', GIcoloringMethod, ...
                     'ResultType', ResultType);
        % xlabel('FA (TEST)','FontWeight','bold');
        ylabel('FA','FontWeight','bold');
        set(gca,'FontSize',18)
        title(sprintf('%s',tn))

        if nrow==1
            if ismember(nc, [1])
                ylabel(ylab,'FontWeight','bold');
            else
                ylabel('');
                set(gca,'YTickLabel',[]);
            end
        else
            if ismember(nc, [1:4:12])
                ylabel(ylab,'FontWeight','bold');
            else
                ylabel('');
                set(gca,'YTickLabel',[]);
            end
        end

    end
    colormap(cmapname)
    h=colorbar;
    set(h, 'Position', [.93 0.115 0.01 .83])
    suptitle({strrep(fnameRoot,'_','\_'), 'GI plots'})
else  % I copied this from the original file for the case of lower diagonal but it needs to be refactored
    warning('Lower diagonal needs to be refactored')
    % Now the data preparation part have been done above, and it returns
    % structs. Adecuate this. Here is all the code, even the call to
    % dr_createCICorrelationMatrix().
    %{
    
% WAHL 2010: correlation colors and CI tables: DATA PREP
% Create tables comparing the WHL test/retest, and the HCP test/retest. The
% first one is not test/retest per se, as it is the same data, but different
% tool
% It is the same code with different data, and it will produce 4 figures. 

% Code to create:
% - 80% CI per correlation, per project. 
% - Automatically obtain the range of common correlations. 
% - Plot the result for some cases
% - Create the resume plot, the same as the previous ones. 

recalculate  = true;
useBootstrap = true;
nRep = 10000;
CIrangeOrVals = 90;
% CIrangeOrVals = [1:100];


tractsOrder = { 'LeftCingulumCingulate'  , 'RightCingulumCingulate'  , ...
                'LeftArcuate'            , 'RightArcuate'            , ...
                'LeftIFOF'               , 'RightIFOF'               , ...
                'LeftILF'                , 'RightILF'                , ...
                'LeftUncinate'           , 'RightUncinate'           , ...
                'LeftCorticospinal'      , 'RightCorticospinal'      };
WahlTractNames = {  'CB left'  , 'CB right'  , ...
                    'AF left'  , 'AF right'  , ...
                    'IFO left' , 'IFO right' , ...
                    'ILF left' , 'ILF right' , ...
                    'UF left'  , 'UF right'  , ...
                    'CST left' , 'CST right' };




cats    = unique(unsMeans.SliceCats);
catsRGB = unique(unsMeans.SliceCatsRGB);
reorder = [4,5,6,1,2,3];
cats    = cats(reorder);
catsRGB = catsRGB(reorder,:);
WahlorigCorr = [
0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0
0.58, 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0
0.57, 0.56, 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0
0.44, 0.55, 0.57, 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0
0.48, 0.72, 0.63, 0.65, 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0
0.42, 0.58, 0.50, 0.60, 0.81, 0   , 0   , 0   , 0   , 0   , 0   , 0
0.44, 0.71, 0.76, 0.59, 0.77, 0.73, 0   , 0   , 0   , 0   , 0   , 0
0.41, 0.46, 0.62, 0.65, 0.59, 0.66, 0.69, 0   , 0   , 0   , 0   , 0
0.60, 0.58, 0.66, 0.49, 0.62, 0.53, 0.66, 0.52, 0   , 0   , 0   , 0
0.43, 0.47, 0.42, 0.38, 0.52, 0.64, 0.57, 0.66, 0.59, 0   , 0   , 0
0.49, 0.45, 0.49, 0.26, 0.44, 0.27, 0.36, 0.22, 0.48, 0.28, 0   , 0
0.33, 0.37, 0.47, 0.42, 0.45, 0.31, 0.32, 0.47, 0.36, 0.25, 0.62, 0];   
WahlorigCorrTable = dr_createCICorrelationMatrix(...
                               unsMeans(unsMeans.SliceCats=="WHL1000",:), ... 
                               'tractsOrder',tractsOrder,'newTractNames',WahlTractNames, ...
                               'useBootstrap',useBootstrap, 'nRep',nRep, 'CIrange',CIrangeOrVals);
WahlorigCorrTable{:,:} = WahlorigCorr;



if recalculate
    for CIrange=CIrangeOrVals
    % OBTAIN CORR AND CIs
        for ns = 1:3
            [corrMatTable,lowerCI,upperCI,bootstrapStat] = dr_createCICorrelationMatrix(...
                                   unsMeans(unsMeans.SliceCats==string(cats(ns)),:), ... 
                                   'tractsOrder',tractsOrder,'newTractNames',WahlTractNames, ...
                                   'useBootstrap',useBootstrap,'nRep',nRep,'CIrange',CIrange);
            corrMatrices{ns}      = corrMatTable;
            corrLowers{ns}        = lowerCI;
            corrUppers{ns}        = upperCI;
            corrBootstrapStat{ns} = bootstrapStat;   
            % Do HCP-T
            [corrMatTable,lowerCI,upperCI,bootstrapStat] = dr_createCICorrelationMatrix(...
                                   unsMeans(unsMeans.SliceCats==string(cats(ns+3)),:), ... 
                                   'tractsOrder',tractsOrder,'newTractNames',WahlTractNames, ...
                                   'useBootstrap',useBootstrap,'nRep',nRep,'CIrange',CIrange);
            corrMatrices{ns+3} = corrMatTable;
            corrLowers{ns+3}   = lowerCI;
            corrUppers{ns+3}   = upperCI;
            corrBootstrapStat{ns+3} = bootstrapStat;
            % Do HCP-RT
            [corrMatTable,lowerCI,upperCI,bootstrapStat] = dr_createCICorrelationMatrix(...
                                   unsMeansrt(unsMeansrt.SliceCats==string(cats(ns+3)),:), ... 
                                   'tractsOrder',tractsOrder,'newTractNames',WahlTractNames, ...
                                   'useBootstrap',useBootstrap,'nRep',nRep,'CIrange',CIrange);
            corrMatrices{ns+6} = corrMatTable;
            corrLowers{ns+6}   = lowerCI;
            corrUppers{ns+6}   = upperCI;
            corrBootstrapStat{ns+6} = bootstrapStat;
        end

        % Add the results from the paper
        AllCorrMatrices = [{WahlorigCorrTable}, corrMatrices]';
        AllLower        = [{WahlorigCorrTable}, corrLowers]';
        AllUpper        = [{WahlorigCorrTable}, corrUppers]';
        AllBStats       = [{WahlorigCorrTable}, corrBootstrapStat]';

        Experiments     = [{'WHLorig'};cellstr(cats(1:3));...
                           strcat(cellstr(cats(4:6)),'TEST');strcat(cellstr(cats(4:6)),'RETEST')];
        sumPlotTable    = table(Experiments,AllCorrMatrices,AllLower,AllUpper,AllBStats);

        UnstackedCorr       = {};
        UnstackedLower      = {};
        UnstackedUpper      = {};
        UnstackedBStats     = {};
        for nt=1:height(sumPlotTable)
            sumPlotTable.AllCorrMatrices{nt} = sumPlotTable.AllCorrMatrices{nt}(2:end,1:end-1);
            sumPlotTable.AllLower{nt} = sumPlotTable.AllLower{nt}(2:end,1:end-1);
            sumPlotTable.AllUpper{nt} = sumPlotTable.AllUpper{nt}(2:end,1:end-1);
            sumPlotTable.AllBStats{nt}= sumPlotTable.AllBStats{nt}(2:end,1:end-1);

            UnstackedCorr  = [UnstackedCorr;  {dr_stackTable(sumPlotTable.AllCorrMatrices{nt})}];
            UnstackedLower = [UnstackedLower; {dr_stackTable(sumPlotTable.AllLower{nt})}];
            UnstackedUpper = [UnstackedUpper; {dr_stackTable(sumPlotTable.AllUpper{nt})}];
            UnstackedBStats= [UnstackedBStats; {dr_stackTable(sumPlotTable.AllBStats{nt})}];
        end
        sumPlotTable.UnstackedCorr  = UnstackedCorr;
        sumPlotTable.UnstackedLower = UnstackedLower;
        sumPlotTable.UnstackedUpper = UnstackedUpper;
        sumPlotTable.UnstackedBStats= UnstackedBStats;

        noTable = sumPlotTable(:,~contains(sumPlotTable.Properties.VariableNames,'All'));
        % Now make the unstacked col numeric
        long = noTable{1,'UnstackedCorr'}{1};
        long.Properties.VariableNames = {'CorName', noTable{1,1}{1}};
        for na = 2:height(noTable)
            tmp = noTable{na,'UnstackedCorr'}{1};
            long.(noTable{na,1}{1}) = tmp.VAL;
        end
        long.Type = cellstr(repmat('Corr', [height(long),1]));
        % Now make the unstacked col numeric
        longL = noTable{1,'UnstackedLower'}{1};
        longL.Properties.VariableNames = {'CorName', noTable{1,1}{1}};
        for na = 2:height(noTable)
            tmp = noTable{na,'UnstackedLower'}{1};
            longL.(noTable{na,1}{1}) = tmp.VAL;
        end
        longL.Type = cellstr(repmat('Lower', [height(long),1]));

        % Now make the unstacked col numeric
        longU = noTable{1,'UnstackedUpper'}{1};
        longU.Properties.VariableNames = {'CorName', noTable{1,1}{1}};
        for na = 2:height(noTable)
            tmp = noTable{na,'UnstackedUpper'}{1};
            longU.(noTable{na,1}{1}) = tmp.VAL;
        end
        longU.Type = cellstr(repmat('Upper', [height(long),1]));

        % Now make the unstacked col numeric
        longB = noTable{1,'UnstackedBStats'}{1};
        longB.Properties.VariableNames = {'CorName', noTable{1,1}{1}};
        for na = 2:height(noTable)
            tmp = noTable{na,'UnstackedBStats'}{1};
            longB.(noTable{na,1}{1}) = tmp.VAL;
        end
        longB.Type = cellstr(repmat('BStats', [height(long),1]));
        

        % Aggregate
        long = [long; longL; longU]; 
        long.Type = categorical(long.Type);
        long.CorName = categorical(long.CorName);
        longcell{CIrange} = long;
    end
    % I did not create the longB for all the CIRange
    % save('allCI1to100correlationsv2.mat','longcell')
else
    load(fullfile(paperReprPath,'local','cache','allCI1to100correlationsv2.mat'))
end

% Create the plots of interests
% Create empty matrix that can go to display_matrix
plotIt     = true; 
CIrange = 90;
long    = longcell{CIrange};
tmp = WahlorigCorrTable; tmp{:,:} = NaN(size(tmp{:,:})); tmp{:,:}(ismissing(tmp))= 0;
TWahlInsideCI = tmp;
TisOverlap    = tmp;
TminUpper     = tmp;
TmaxLower     = tmp;
TmidCI        = tmp;

% PLOT LOWER DIAGONAL (substitutes both figures)
if plotIt
    mrvNewGraphWin('Lower Diag. FA Correlation CI Plots');
    % figHdl = figure('Name','FA profiles', ...
    %                 'NumberTitle','off', ...
    %                 'visible',   'on', ...
    %                 'color','w', ...
    %                 'Units','pixel', ...
    %                 'Position',[0 0 1900 1100]);
    set(gcf,'WindowStyle','normal');
    set(gcf, 'Units', 'inches', 'OuterPosition', [0, 0, 19, 12]);
    set(gcf,'color','w');
end
ncol=11;nrow = 11;
cPos         = tmp(2:end,1:end-1);
cPos{:,:}    = reshape([1:size(cPos{:,:},1) .* size(cPos{:,:},2)],[size(cPos{:,:},1),size(cPos{:,:},2)])';
% Figure paper plots 4 and 9 

referenceColumn = 'WHLorig'; % 'WHLorig' | 'NoReference'
versus          = 'WHL1000';
plotIndex       = false;
showLegend      = false; 
showXnames      = false; 
refLine         = true; 
GIcoloringMethod = 'GRbars'; % 'circle', 'GIband', 'GRband', 'background','none','bars', 'GRbars'
% takeout =  {'WHLorig','WHL1000', ...
%             'YWM1000','YWM2000', ...
%             'HCP3000TEST','HCP3000RETEST', ...
%             'HCP1000TEST','HCP1000RETEST'};
takeout =  {'YWM1000','YWM2000', ...
            'HCP2000TEST','HCP2000RETEST', ...
            'HCP3000TEST','HCP3000RETEST', ...
            'HCP1000TEST','HCP1000RETEST'};
% takeout =  {'WHLorig'};
% takeout =  {'WHLorig','WHL1000', ...
%             'YWM1000','YWM2000'};
% takeout =  {'WHLorig','YWM2000', ...
%             'HCP2000TEST','HCP2000RETEST', ...
%             'HCP3000TEST','HCP3000RETEST'};
GIlist = [];
for nc=1:length(unique(long.CorName))
    % Calculate the position
    corname  = string(long.CorName(nc));
    cornames = split(corname,'_');
    bat = cornames{1};
    bi  = cornames{2};
    
    if plotIt; subplot(ncol,nrow,cPos{bat,bi}); end    
    % DO the calculation/plots
    [WahlInsideCI, isOverlap, minUpper, maxLower, midCI, GI]=dr_compareCI(long, longB, corname,...
                                                    'referenceColumn',referenceColumn,...
                                                    'refLine' , refLine,'plotIndex',plotIndex, ...
                                                    'plotIt',plotIt,'showLegend',showLegend, ...
                                                    'takeout',takeout, 'showXnames',showXnames, ...
                                                    'GIcoloringMethod',GIcoloringMethod, ...
                                                    'ResultType','correlation');
    if ismember(cPos{bat,bi}, cPos{size(cPos,1),:});xlabel(bi,'fontsize',14,'color','k','fontweight','bold');end
    if ismember(cPos{bat,bi}, cPos{:,1});ylabel(bat,'fontsize',14,'color','k','fontweight','bold');end                  

    % Now populate the tables

    % fprintf('%s_%s :: %f\n',bat, bi, mc.mean(ii))
    TWahlInsideCI{bat,bi} = WahlInsideCI;
    TisOverlap{bat,bi}    = isOverlap;
    TminUpper{bat,bi}     = minUpper;
    TmaxLower{bat,bi}     = maxLower;
    TmidCI{bat,bi}        = midCI;
    GIlist = [GIlist,GI];
end
if plotIndex
    markInd = '_withIndex';
else
    markInd = '_noIndex';
end
mainTitle = sprintf('GI%s_%svs%s_CI%d_%s', markInd, referenceColumn, versus, CIrange,GIcoloringMethod);
csvwrite(fullfile(saveItHere,[mainTitle '.csv']), GIlist);
if plotIt
    % Add a colorbar
    colormap(cmapname)
    h=colorbar;
    % set()
    set(h, 'Position', [.93 0.115 0.01 .83])
%     for i=1:4
%           pos=get(ax(i), 'Position');
%           set(ax(i), 'Position', [pos(1) pos(2) 0.85*pos(3) pos(4)]);
%     end

    
    suptitle(strrep(mainTitle,'_','\_'))
    set(gcf,'color','w');

    saveas(gcf,fullfile(saveItHere, [mainTitle '.png']),'png');
    saveas(gcf,fullfile(saveItHere, [mainTitle '.svg']),'svg');
    close(gcf)
end

    %}
end








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
