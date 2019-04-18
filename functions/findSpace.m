function [band,bigfig] = findSpace(groupStats,varargin)
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
% winSizePix     : array of 4 floats [0,0,1900,1100]
%
% Examples:
%{
%}
% 
% GLU Vistalab, 2018


%% PARSE INPUTS
p = inputParser;

addRequired(p, 'groupStats');

addOptional(p, 'fnameRoot'       , "changeThisName" , @isstring);
addOptional(p, 'saveItHere'      , "~/tmp"          , @isstring);
addOptional(p, 'savePng'         , false            , @islogical);
addOptional(p, 'saveSvg'         , false            , @islogical);
addOptional(p, 'nrowcol'         , [4,3]            , @isfloat);
addOptional(p, 'winSizePix'      , []               , @isfloat);  % Almost full screen in laptop, [0 0 1900 1100]
addOptional(p, 'winSizeInch'     , []               , @isfloat);  % [0 0 10 14]
addOptional(p, 'ylab'            , 'nothing'        , @ischar);
addOptional(p, 'xlab'            , 'nothing'        , @ischar);
addOptional(p, 'ylims'           , [0,1]            , @isfloat);
addOptional(p, 'xlims'           , [0,1]            , @isfloat);
addOptional(p, 'showXnames'      , false            , @islogical);
addOptional(p, 'LineColor'       , [0 0 0]          , @isfloat);
addOptional(p, 'LineStyle'       , '-.'             , @ischar);
addOptional(p, 'showLineForm'    , false            , @islogical);
addOptional(p, 'useStdAllBS'     , 'Std'            , @ischar);
addOptional(p, 'howManySD'       , [2,1]            , @isfloat);
addOptional(p, 'CIperc'          , [95,68]          , @isfloat);
addOptional(p, 'nReps'           , 1000             , @isfloat);
addOptional(p, 'plotSubj'        , true             , @islogical);
addOptional(p, 'SubjValues'      , table()          , @istable);
addOptional(p, 'findBand'        , true             , @islogical);
parse(p,groupStats,varargin{:});

fnameRoot       = p.Results.fnameRoot;
saveItHere      = p.Results.saveItHere;
savePng         = p.Results.savePng;
saveSvg         = p.Results.saveSvg;
nrowcol         = p.Results.nrowcol;
winSizePix      = p.Results.winSizePix;
winSizeInch     = p.Results.winSizeInch;
ylab            = p.Results.ylab;
xlab            = p.Results.xlab;
ylims           = p.Results.ylims;
xlims           = p.Results.xlims;
showXnames      = p.Results.showXnames;
LineColor       = p.Results.LineColor;
LineStyle       = p.Results.LineStyle;
showLineForm    = p.Results.showLineForm;
useStdAllBS     = p.Results.useStdAllBS;
howManySD       = p.Results.howManySD;
CIperc          = p.Results.CIperc;
nReps           = p.Results.nReps;
plotSubj        = p.Results.plotSubj;
SubjValues      = p.Results.SubjValues;
findBand        = p.Results.findBand;

%% PREPARE THE DATA

if ~istable(groupStats)
    error("groupStats needs to be a table")
end
GS       = groupStats;

band = [0,1];
if findBand
   unsMeans = SubjValues;
end

% Depending on the method, calculate the lines to be in or out
switch lower(useStdAllBS)
case {'std','sd'}
    fmean   = @(X) mean(X,'omitnan');
    fstd    = @(X) std(X,'omitnan');
    Mean_Cell= cellfun(fmean, GS{:,'residuals'}, 'Uni',0);
    Std_Cell = cellfun(fstd, GS{:,'residuals'}, 'Uni',0);
    GS.Mean     = cell2mat(Mean_Cell); 
    GS.Std      = cell2mat(Std_Cell); 
case {'bs','bootstrap'}
    
    for nc = 1:height(GS)
        tn     = string(GS.CorName(nc));
        values = GS.residuals{GS.CorName==tn};
        if length(CIperc)==1
            [Mean1,lowerQCI1,upperQCI1,~,lowerLowerQCI1,~,~,~,upperUpperQCI1]=dr_bootstrapDistribution(...
                                                                values, ...
                                                                'nReps', nReps, ...
                                                                'perc', CIperc, ...
                                                                'grandMean', true);
        else
            [Mean1,lowerQCI1,upperQCI1,~,lowerLowerQCI1,~,~,~,upperUpperQCI1]=dr_bootstrapDistribution(...
                                                                values, ...
                                                                'nReps', nReps, ...
                                                                'perc', CIperc(1), ...
                                                                'grandMean', true);
            [Mean2,lowerQCI1,upperQCI1,~,lowerLowerQCI2,~,~,~,upperUpperQCI2]=dr_bootstrapDistribution(...
                                                                values, ...
                                                                'nReps', nReps, ...
                                                                'perc', CIperc(2), ...
                                                                'grandMean', true);
            GS.Mean2(GS.CorName==tn)     = Mean2;
            GS.low2(GS.CorName==tn)   = lowerQCI2;
            GS.up2(GS.CorName==tn)     = upperQCI2;
        end
        
                                                            
        GS.Mean1(GS.CorName==tn)     = Mean1;
        GS.low1(GS.CorName==tn)   = lowerQCI1;
        GS.up1(GS.CorName==tn)     = upperQCI1;

    end
case {'all','allvalues'}
    Max_Cell = cellfun(@max, GS{:,'residuals'}, 'Uni',0);
    Min_Cell = cellfun(@min, GS{:,'residuals'}, 'Uni',0);
    GS.UpperLim = cell2mat(Max_Cell); 
    GS.LowerLim = cell2mat(Min_Cell); 
otherwise
    error('This case has not been implemented yet, select: std, bs or all')
end


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




bigfig = figure('Name',fnameRoot, ...
                'NumberTitle','off', ...
                'visible',   'on', ...
                'color','w', ...
                'WindowStyle','normal', ...
                'Units',units, ...
                winPosType,winSize);
nrow=nrowcol(1); ncol=nrowcol(2);
for nc = 1:height(GS)
    % Calculate the position
    tn  = string(GS.CorName(nc));
    sp = subplot(nrow,ncol,nc);
    
    % Do the calculations/plots
    b         = GS.b(GS.CorName==tn);
    slope     = GS.slope(GS.CorName==tn);
    % Depending on the method, calculate the lines to be in or out
    plotContour2 = false;
    switch lower(useStdAllBS)
        case {'std','sd'}
            if length(howManySD)==1
                alto1 = GS.Mean(GS.CorName==tn) + howManySD * GS.Std(GS.CorName==tn);
                bajo1 = GS.Mean(GS.CorName==tn) - howManySD * GS.Std(GS.CorName==tn);
            else
                alto1 = GS.Mean(GS.CorName==tn) + howManySD(1) * GS.Std(GS.CorName==tn);
                bajo1 = GS.Mean(GS.CorName==tn) - howManySD(1) * GS.Std(GS.CorName==tn);
                alto2 = GS.Mean(GS.CorName==tn) + howManySD(2) * GS.Std(GS.CorName==tn);
                bajo2 = GS.Mean(GS.CorName==tn) - howManySD(2) * GS.Std(GS.CorName==tn);
                plotContour2 = true;
                contours2 = [alto2,bajo2];
            end
        case {'bs','bootstrap'}
            alto1 = GS.up1(GS.CorName==tn);
            bajo1 = GS.low1(GS.CorName==tn);
            if length(CIperc)==2
                alto2 = GS.up2(GS.CorName==tn);
                bajo2 = GS.low2(GS.CorName==tn);
                plotContour2 = true; 
                contours2 = [alto2,bajo2];
            end
        case {'all','allvalues'}
            alto1     = GS.UpperLim(GS.CorName==tn);
            bajo1     = GS.LowerLim(GS.CorName==tn);
        otherwise
            error('This case has not been implemented yet, select: std, bs or all')
    end
    L         = linspace(0.01,1);
    R         = linspace(0.01,1);
    [X,Y]     = meshgrid(L,R);
    Z         = Y - (b + slope * X);
    [c,h]     = contour(X,Y,Z, [alto1,bajo1], ...
                        'LineColor', LineColor, ...
                        'LineStyle', LineStyle, ...
                        'LineWidth', 2); hold on;
                    
     % Plot valid band values for right at position 0.5
     Z05 = Y - (b + slope * 0.5);
     topright = max(unique(Y(Z05 < alto1 & Z05 > bajo1)));
     botright = min(unique(Y(Z05 < alto1 & Z05 > bajo1)));
     
     
    %      h1=plot(0.5*[1,1],[0 topright],'Color',[.3 .3 .3],'LineStyle','-.','LineWidth',2);hold on;
    %      jbfill([0.01:0.01:1], ...
    %              topright * ones(size([0.01:0.01:1])), ...
    %              botright * ones(size([0.01:0.01:1])),...
    %              'k','k',0,0.25);
     
    %       text(xlims(1)+.01,ylims(2)-0.05, ...
    %           sprintf('RightFA'), ...
    %           'fontsize',14,'color',[.3, .3, .3],'FontWeight','bold');
      text(xlims(1)+.01,ylims(2)-0.03, ...
          sprintf('Right CI(95%%)=%.2f FA',topright-botright),...
          'fontsize',14,'color',[.3, .3, .3],'FontWeight','bold');
                    
                    
                    
                    
                    
    if plotContour2
        hold on;
        [c,h]     = contour(X,Y,Z, contours2, ...
                        'LineColor', LineColor, ...
                        'LineStyle', ':', ...
                        'LineWidth', 1); 
    end
    % Plot the slope formula
    if showLineForm
        text(xlims(1)+.05,ylims(1)-0.05,sprintf('Right=%.2f+%.2fxLeft',b,slope),'fontsize',12,'color',[.3, .3, .3],'FontWeight','bold');
    end

    if plotSubj
        for ns = 1:height(SubjValues)
            if findBand
                color = table2array(SubjValues{ns,'SliceCatsRGB'});
            else
                if ns==1
                    color=[80 220 100]/255;
                else
                    color=[237 41 57]/255;
                end
            end
            Left  = SubjValues{ns,strcat('Left',tn)};
            Right = SubjValues{ns,strcat('Right',tn)};
            plot(Left, Right, '.', 'color',color, 'MarkerSize',8)
        end
        
    end
    
    
    
    set(gca,'FontSize',18)
    title(sprintf('%s',tn))
    grid on; axis equal;
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
    xlim(xlims);ylim(ylims);
    xlabel(xlab,'FontWeight','bold');
    % xticks([xlims(1):0.1:xlims(2)]);
    
end
suptitle({strrep(fnameRoot,'_','\_')})


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
