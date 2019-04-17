function  [RefInsideCI,isOverlap,minUpper,maxLower,meanCI,GI]=dr_compareCI(...
                                                        t, B, corname, varargin)
%
% 
%
% Syntax:
%     corrMatTable = dr_compareCI(t, ...)
%
% Description:
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

addRequired(p, 't'               , @istable);
addRequired(p, 'B'               , @istable);
addRequired(p, "corname"         , @isstring);

addOptional(p, 'CIrange'         , 90        , @isnumeric);
addOptional(p, 'referenceColumn' , 'WAHLorig', @ischar);
addOptional(p, 'plotIt'          , false     , @islogical);
addOptional(p, 'plotIndex'       , false     , @islogical);
addOptional(p, 'refLine'         , false     , @islogical);
addOptional(p, 'showLegend'      , false     , @islogical);
addOptional(p, 'showXnames'      , false     , @islogical);
addOptional(p, 'includeExp'      , {}        , @iscellstr);
addOptional(p, 'GIcoloringMethod', 'none'    , @ischar);
addOptional(p, 'cmapname'        , 'copper'  , @ischar);
addOptional(p, 'ResultType'      , 'corr'    , @ischar);
addOptional(p, 'normResidual'    , false     , @islogical);
addOptional(p, 'groupStats'      , table()   , @istable);
parse(p,t,B,corname, varargin{:});

CIrange         = p.Results.CIrange;
referenceColumn = p.Results.referenceColumn;
plotIt          = p.Results.plotIt;
plotIndex       = p.Results.plotIndex;
refLine         = p.Results.refLine;
showLegend      = p.Results.showLegend;
includeExp      = p.Results.includeExp;
showXnames      = p.Results.showXnames;
GIcoloringMethod= p.Results.GIcoloringMethod;
cmapname        = p.Results.cmapname;
ResultType      = p.Results.ResultType;
normResidual    = p.Results.normResidual;
groupStats      = p.Results.groupStats;

switch lower(ResultType)
    case {'fa',"FA"}  % add other values with 0-1 ranges
        LowerRange = 0;
        UpperRange = 1;
    case {'corr','correlation',"corr","correlation",'regressionslope'}
        LowerRange = -1;
        UpperRange = 1;
    case {'faregressionresiduals'}
        LowerRange = -0.2;
        UpperRange = 0.2;
    case {'faregressionnormresiduals'}
        LowerRange = -5;
        UpperRange = 5;        
    otherwise
        error('Please add the range of this result type to dr_compareCI.m')
end


% Do the thing
% Reference (in WAHLorig the three values will be the same)
RefInsideCI = false;
noreference = strcmp(lower(referenceColumn),'noreference');

if strcmp(t.Properties.VariableNames{1},'TractNames')
    names = t.Properties.VariableNames;
    names{1} = 'CorName';
    t.Properties.VariableNames = names;
    B.Properties.VariableNames = names;
end

% Obtain Values for the CIs
restCols        = ~ismember(t.Properties.VariableNames, {'CorName','Type',referenceColumn});
restColumnNames = t.Properties.VariableNames(restCols);

% Obtain Values for the bootstraps
restColsB        = ~ismember(B.Properties.VariableNames, {'CorName','Type',referenceColumn});
restColumnNamesB = B.Properties.VariableNames(restColsB);

if ~isequal(restColumnNames, restColumnNamesB); error('Not using the same Experiments');end


if isempty(includeExp)
    error('Include some experiments for the GI calculation')
else
    restColumnNames  = restColumnNames(ismember(restColumnNames, includeExp));
    restColumnNamesB = restColumnNamesB(ismember(restColumnNamesB, includeExp));
end
corrVals     = t{t.CorName==corname & t.Type=="Corr" ,restColumnNames};
upperVals    = t{t.CorName==corname & t.Type=="Upper",restColumnNames};
lowerVals    = t{t.CorName==corname & t.Type=="Lower",restColumnNames};
if ismember({'ALL'},includeExp)
    Bvalues      = B{B.CorName==corname, ...
                     restColumnNamesB(~ismember(restColumnNamesB,{'ALL'})) };
    BvaluesALL   = B{B.CorName==corname, ...
                     restColumnNamesB(ismember(restColumnNamesB,{'ALL'})) };
    valuesALL    = cell2mat(BvaluesALL)';
else
    Bvalues      = B{B.CorName==corname, restColumnNamesB};
end


values = cell2mat(Bvalues)';

[Mean    ,lowerQCI     ,upperQCI, ...
         lowerMean,lowerLowerQCI,lowerUpperQCI, ...
         upperMean,upperLowerQCI,upperUpperQCI]= dr_bootstrapDistribution(values, ...
                                                            'nReps', 500, ...
                                                            'perc', 90, ...
                                                            'grandMean', true);

if ~(length(restColumnNames)==length(upperVals))
    error('Missmatch in experiment numbers'); end
[GI,meanCI,isOverlap,minUpper,maxLower,OR] = dr_GI(corrVals, upperVals, lowerVals, ...
                                                        'LowerRange', LowerRange, ...
                                                        'UpperRange', UpperRange);

if ~noreference
    refCorr  = t{t.CorName==corname & t.Type=="Corr" ,referenceColumn};
    refUpper = t{t.CorName==corname & t.Type=="Upper",referenceColumn};
    refLower = t{t.CorName==corname & t.Type=="Lower",referenceColumn}; 

    if refCorr >= maxLower && refCorr <= minUpper; RefInsideCI = true; end
end

if plotIt
    % Prepare data for plotting
    X     = [1:length(upperVals); 1:length(upperVals)];
    Y     = [upperVals; lowerVals];
    
    % Mark the band where the residuals shouldn't be
    if normResidual
        jbfill([0.8, length(upperVals)+.2],[-2,-2],[2,2],[.1,.1,.1],[.2,.2,.2],0,0.1);
    end
    
    
    if refLine && strcmpi(GIcoloringMethod, 'grbars')
        l3 = line([0.8, length(upperVals)+.2],[refCorr,refCorr],'linewidth',2,'color','k','linestyle','-.'); hold on;
    else
        l1 = line([0.8, length(upperVals)+.2],[upperUpperQCI,upperUpperQCI],'linewidth',2,'color',[0.3, 0.3, 0.3],'linestyle',':');hold on;
        l2 = line([0.8, length(upperVals)+.2],[lowerLowerQCI,lowerLowerQCI],'linewidth',2,'color',[0.3, 0.3, 0.3],'linestyle',':');
    end
        
    xlim([0.8, length(upperVals)+.2]); 
    ylim([LowerRange UpperRange]);
    cmap = colormap(cmapname);
    cind = ceil(GI*length(cmap));
    if cind==0; cind=1; end
    % ylabel('Correlation'); 
    % title([ char(strrep(corname,'_','\_'))])
    % HERE COMES THE COLORING PART
    yPos = 0.95;if minUpper > 0.7; yPos = 0.1;end
    switch lower(GIcoloringMethod)
        case {'grband'}
            if isOverlap
                if length(upperVals)==1 && ~RefInsideCI 
                    jbfill([0.8, length(upperVals)+.2],[minUpper,minUpper],[maxLower,maxLower],'r','k',0,0.3);
                else
                    jbfill([0.8, length(upperVals)+.2],[minUpper,minUpper],[maxLower,maxLower],[0.4660, 0.6740, 0.1880],'k',0,0.3);
                end
            else
                jbfill([0.8, length(upperVals)+.2],[minUpper,minUpper],[maxLower,maxLower],'r','k',0,0.3);
            end
            plot(X,Y,'linewidth',3,'color','k'); 
        case {'giband'}
            jbfill([0.8, length(upperVals)+.2],[max(upperVals),max(upperVals)],[min(lowerVals),min(lowerVals)],cmap(cind,:),'k',0,1);
            l3 = line([0.8, length(upperVals)+.2],[grandMean,grandMean],'linewidth',1,'color','k','linestyle','-');
            plot(X,Y,'linewidth',3,'color','k'); 
        case {'circle'}
            plot(X,Y,'linewidth',3,'color','k'); 
        case {'bars'}
            l3 = line([0.8, length(upperVals)+.2],[Mean,Mean],'linewidth',1,'color','k','linestyle','-');
            plot(X,Y,'linewidth',12,'color',cmap(cind,:)); 
            if plotIndex
                text(1,yPos,sprintf('G_%i(%i%%)=%.2f%',length(upperVals),CIrange,GI),'fontsize',12,'color',[.3, .3, .3]);
            end
        case {'barsdots'}
            if ismember({'ALL'},includeExp)
                ALLpos = find(ismember(includeExp,{'ALL'}),1);
                XB    = repmat([1:(length(upperVals)-1)], [length(values)/(length(upperVals)-1),1]);
                % Add some jittering for plotting
                XBjitt    = XB + 0.25*(rand(size(XB))-.5);
                YB    = reshape(values, [length(values)/(length(upperVals)-1),length(upperVals)-1]);
                
                YBALL = valuesALL;
                XBALL = ALLpos * ones(size(YBALL));
                XBALLjitt = XBALL + 0.3*(rand(size(XBALL))-.5);
            else
                XB    = repmat([1:length(upperVals)], [length(values)/length(upperVals),1]);
                % Add some jittering for plotting
                XB    = XB + 0.3*(rand(size(XB))-.5);
                YB    = reshape(values, [length(values)/length(upperVals),length(upperVals)]);
            end
            if ~isclose(YB(:,1),Bvalues{1}')
                error('There was a problem converting the data')
            end

            l3 = line([0.8, length(upperVals)+.2],[Mean,Mean],'linewidth',1,'color','k','linestyle','-');
            plot(X,Y,'linewidth',12,'color',cmap(cind,:)); hold on;
            plot(XB,YB,'k.','MarkerSize',8);
            plot(XBALL,YBALL,'w.','MarkerSize',6);
            if plotIndex
                text(1,yPos,sprintf('G_%i(%i%%)=%.2f%',length(upperVals),CIrange,GI),'fontsize',12,'color',[.3, .3, .3]);
            end
        case {'grbars'}
            green = [0,204,102]/255; red = [249,73,73]/255;
            if RefInsideCI
                plot(X,Y,'linewidth',12,'color',green); 
            else
                plot(X,Y,'linewidth',12,'color',red); 
            end         
        case {'background'}
            set(gca,'color',cmap(cind,:))
            plot(X,Y,'linewidth',3,'color','k'); 
        case {'none'}
            plot(X,Y,'linewidth',3,'color','k'); 
    end
    
    if plotIndex
        text(1,UpperRange,sprintf('G_%i(%i%%)=%.2f%',length(upperVals),CIrange,GI),'fontsize',12,'color',[.3, .3, .3]);
    end   
    % Plot the slope formula
    if ~isempty(groupStats)
        b     = groupStats{groupStats.CorName==corname,'b'};
        slope = groupStats{groupStats.CorName==corname,'slope'};
        text(1,.9*UpperRange,sprintf('Right=%.2f+%.2fLeft',b,slope),'fontsize',14,'color',[.3, .3, .3],'FontWeight','bold');
    end
    
    
    if showXnames
        % xticks(1:length(upperVals)); xtickangle(45); 
        % xticklabels(restColumnNames);
        tt = text(1:length(upperVals),LowerRange * ones([1,length(upperVals)]), restColumnNames);
        set(tt,'Rotation',90,'fontsize',14);  % ,'FontWeight','bold'
        set(gca,'xtick',[])
        set(gca,'xcolor',[1 1 1]);
        box off;
    else
        set(gca,'xtick',[]);
        % set(gca,'ytick',[]);
        box off;
        set(gca,'xcolor',[1 1 1]);
        set(gca,'linewidth',2)
        % ylabel('FA', 'FontWeight','bold'); % xlabel('Divisions')
        set(gca,'xtick',[])
        
    end
    if refLine
        if showLegend; legend([l1,l2,l3], {'minimum upper C.I.','maximum lower C.I.',['Reference ' referenceColumn]});end;
    else
        if showLegend; legend([l1,l2], {'minimum upper C.I.','maximum lower C.I.'});end;
    end
    
    
end

    
    
    
    
end
