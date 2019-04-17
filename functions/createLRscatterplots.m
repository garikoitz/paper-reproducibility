function [TRTtests,allCOVS] = createLRscatterplots(unsMeans,includeExp,varargin)
%CREATEPROFILE Summary of this function goes here
%   Detailed explanation goes here
%
% 
% Syntax:
%     createLRscatterplotsGIplot(dt)
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
% winSizePix     : array of 4 floats [0,0,1900,1100]
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
addOptional(p, 'winSizePix'     , []               , @isfloat);  % Almost full screen in laptop, [0 0 1900 1100]
addOptional(p, 'winSizeInch'    , []               , @isfloat);  % [0 0 10 14]
addOptional(p, 'ylab'           , 'nothing'        , @ischar);
addOptional(p, 'xlab'           , 'nothing'        , @ischar);
addOptional(p, 'ylimits'        , [0,1]            , @isfloat);
addOptional(p, 'xlimits'        , [0,1]            , @isfloat);

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
winSizePix      = p.Results.winSizePix;
winSizeInch     = p.Results.winSizeInch;
ylab            = p.Results.ylab;
xlab            = p.Results.xlab;
ylimits         = p.Results.ylimits;
xlimits         = p.Results.xlimits;

%% PREPARE THE DATA
if ~iscellstr(includeExp) || isempty(includeExp)
    error("includeExp need to be a cellstr with at least one experiment")
end

for ncat = 1:length(includeExp)

    t = unsMeans(unsMeans.SliceCats==includeExp{ncat},:);
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
    newTractNames = {  'CBleft'  , 'CBright'  , ...
                        'AFleft'  , 'AFright'  , ...
                        'IFOleft' , 'IFOright' , ...
                        'ILFleft' , 'ILFright' , ...
                        'UFleft'  , 'UFright'  , ...
                        'CSTleft' , 'CSTright' };
    newTractNames = {};

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

    %% CREATE FIGURE AND PLOT

    if ~isempty(winSizePix) && isempty(winSizeInch) 
        units      = 'pixel';
        winSize    = winSizePix;
        winPosType = 'Position';
    elseif isempty(winSizePix) && ~isempty(winSizeInch) 
        units = 'inches';
        winSize = winSizeInch;
        winPosType = 'OuterPosition';
    else
        units      = 'pixel';
        winSizePix = [0 0 1900 1100];
        winPosType = 'Position';
        warning('Using default window size, specify winSizePix or winSizePix (and not both)')
    end


        bigfig = figure('Name',strcat(fnameRoot,'_',cat), ...
                        'NumberTitle','off', ...
                        'visible',   'on', ...
                        'color','w', ...
                        'WindowStyle','normal', ...
                        'Units',units, ...
                        winPosType,winSize);


        nrow=nrowcol(1); ncol=nrowcol(2);
        for nc = 1:length(newStructs)
            subplot(nrow,ncol,nc);
            tn = newStructs{nc};
            Left     = t{:,strcat("Left",tn)};
            Right    = t{:,strcat("Right",tn)};
            
            % Esto es para una prueba que hicimos sin el valor constante
            % slope = inv(Left' * Left) * (Left'* Right);
            
            catcolor = t.SliceCatsRGB{1,:};
            title(sprintf('%s',tn))
            identityLine(gca);
            scatter(Left,Right,20,catcolor,'filled');hold on;
            % scatter(b2000T,b2000RT,18,catcolors{3,:});
            % scatter(b3000T,b3000RT,20,catcolors{5,:},'filled');
            % xticks([0.1:0.1:0.6]);yticks([0.1:0.1:0.6]);
            set(gca,'FontSize',18)
            title(sprintf('%s',tn))
            axis equal;

            % Here's the line
            p = line(xlimits,ylimits, 'color',[.5 .5 .5], 'linestyle','--', 'linewidth',2);
            % Set line properties.  These probably want to come in as an argument
            grid on
            xlim(xlimits);
            ylim(ylimits);
            lsl = lsline;
            set(lsl, 'color', catcolor, 'linewidth',2)
            

            if nrow==1
                if ismember(nc, [1])
                    ylabel(ylab,'FontWeight','bold');
                    yticks([0.4, 0.6])
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
            xlabel(xlab,'FontWeight','bold');
        end
        suptitle({strrep(fnameRoot,'_','\_'), cat})

    %% SAVE
    if ~exist(saveItHere,'dir')
        mkdir(saveItHere)
    end
    if saveSvg
        saveas(gcf,fullfile(saveItHere, strcat(fnameRoot,'_',cat,'.svg')),'svg');
    end
    if savePng
        saveas(gcf,fullfile(saveItHere, strcat(fnameRoot,'_',cat,'.png')),'png');
    end
end