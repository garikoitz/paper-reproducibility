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
parse(p,t,B,corname, varargin{:});

CIrange         = p.Results.CIrange;
referenceColumn = p.Results.referenceColumn;
plotIt          = p.Results.plotIt;
plotIndex       = p.Results.plotIndex;
refLine         = p.Results.refLine;
showLegend      = p.Results.showLegend;
includeExp         = p.Results.includeExp;
showXnames      = p.Results.showXnames;
GIcoloringMethod= p.Results.GIcoloringMethod;
cmapname        = p.Results.cmapname;
ResultType      = p.Results.ResultType;


switch lower(ResultType)
    case {'fa',"FA"}  % add other values with 0-1 ranges
        LowerRange = 0;
        UpperRange = 1;
    case {'corr','correlation',"corr","correlation"}
        LowerRange = -1;
        UpperRange = 1;
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
Bvalues      = B{B.CorName==corname, restColumnNamesB};

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
    X     = [1:length(upperVals); 1:length(upperVals)];
    Y     = [upperVals; lowerVals];
    
    if refLine && strcmpi(GIcoloringMethod, 'grbars')
        l3 = line([0.8, length(upperVals)+.2],[refCorr,refCorr],'linewidth',2,'color','k','linestyle','-.'); hold on;
    else
        l1 = line([0.8, length(upperVals)+.2],[upperUpperQCI,upperUpperQCI],'linewidth',2,'color',[0.3, 0.3, 0.3],'linestyle',':');hold on;
        l2 = line([0.8, length(upperVals)+.2],[lowerLowerQCI,lowerLowerQCI],'linewidth',2,'color',[0.3, 0.3, 0.3],'linestyle',':');
    end
        
    xlim([0.8, length(upperVals)+.2]); ylim([0 1]);
    cmap = colormap(cmapname);
    cind = ceil(GI*length(cmap));
    if cind==0; cind=1; end
    if showXnames
        xticks(1:length(upperVals)); xtickangle(45); 
        xticklabels(restColumnNames);
    else
        set(gca,'xtick',[]);
        % set(gca,'ytick',[]);
        box off;
        set(gca,'xcolor',[1 1 1]);
        set(gca,'linewidth',2)
        % ylabel('FA', 'FontWeight','bold'); % xlabel('Divisions')
        set(gca,'xtick',[])
        yPos = 0.95;if minUpper > 0.7; yPos = 0.1;end
    end
    % ylabel('Correlation'); 
    % title([ char(strrep(corname,'_','\_'))])
    % HERE COMES THE COLORING PART
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
        case {'grbars'}
            green = [0,204,102]/255; red = [249,73,73]/255;
            if RefInsideCI
                plot(X,Y,'linewidth',12,'color',green); 
            else
                plot(X,Y,'linewidth',12,'color',red); 
            end
            if plotIndex
                text(1,yPos,sprintf('G_%i(%i%%)=%.2f%',length(upperVals),CIrange,GI),'fontsize',12,'color',[.3, .3, .3]);
            end            
        case {'background'}
            set(gca,'color',cmap(cind,:))
            plot(X,Y,'linewidth',3,'color','k'); 
        case {'none'}
            plot(X,Y,'linewidth',3,'color','k'); 
    end
    
    
    if refLine
        if showLegend; legend([l1,l2,l3], {'minimum upper C.I.','maximum lower C.I.',['Reference ' referenceColumn]});end;
    else
        if showLegend; legend([l1,l2], {'minimum upper C.I.','maximum lower C.I.'});end;
    end
    
    
end

    
    
    
    
end
