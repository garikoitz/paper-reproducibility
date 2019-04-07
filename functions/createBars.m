function createBars(dt,varargin)
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
% meanOrSd  : string
% savePng   : boolean
% saveSvg   : boolean
% WahlOrder : boolean
% HCPTRT    : boolean
%
% Examples:
%{
%}
% 
% GLU Vistalab, 2018


%% 0.- Parse inputs
p = inputParser;

addRequired(p, 'dt');

addOptional(p, 'fnameRoot'   , "changeThisName" , @isstring);
addOptional(p, 'saveItHere'  , "~/tmp"          , @isstring);
addOptional(p, 'meanOrSd'    , "mean"           , @isstring);
addOptional(p, 'savePng'     , false            , @islogical);
addOptional(p, 'saveSvg'     , false            , @islogical);
addOptional(p, 'WahlOrder'   , false            , @islogical);
addOptional(p, 'HCPTRT'      , false            , @islogical);
parse(p,dt,varargin{:});

fnameRoot   = p.Results.fnameRoot;
saveItHere  = p.Results.saveItHere;
meanOrSd    = p.Results.meanOrSd;
savePng     = p.Results.savePng;
saveSvg     = p.Results.saveSvg;
WahlOrder   = p.Results.WahlOrder;
HCPTRT      = p.Results.HCPTRT;

%% Prepare the data
cats      = categories(dt.SliceCats);
catcolors = unique(dt.SliceCatsRGB);
tracts    = unique(dt.Struct);
% Wahl Order
if WahlOrder
    newOrdering = [3,4,  11,12,  5,6,   7,8,   9,10,   1,2];
    tracts    = tracts(newOrdering);
end
% HCPTRT
if HCPTRT
    catcolors = catcolors([1,1,2,2,3,3],:);
    linestyles= {':','-',':','-',':','-'};
end



omean   = @(x) mean(x,'omitnan');
ostd    = @(x) std(x,'omitnan');
left_color = [0 0 .55]; right_color = [0.3 0.3 0.3];


tmp     = varfun(ostd,dt,'GroupingVariables',{'Struct','SliceCats'},'InputVariables',{'meanVal'});
replRes = varfun(omean,dt,'GroupingVariables',{'Struct','SliceCats'},'InputVariables',{'meanVal'});
replRes.Properties.VariableNames{'Fun_meanVal'} = 'Mean';
replRes.SD  = tmp.Fun_meanVal;
replRes     = replRes(:,[1,2,4,5]);
% Unstack it to show it as consecutive tables. Reorder to replicate
replRes     = unstack(replRes,{'Mean','SD'},{'SliceCats'});
if WahlOrder
    replRes     = replRes(newOrdering,:);
end
% Table copied from Wahl 2010 paper
replRes.Mean_WHLorig = [0.565,0.522,0.517,0.490,0.549,0.534,0.510,0.497,0.489,0.470,0.587,0.575]';
replRes.SD_WHLorig   = [0.028,0.028,0.023,0.03,0.026,0.024,0.028,0.026,0.022,0.024,0.023,0.021]';

if HCPTRT
    replResSorted = replRes(:,[1, 3,9,   2,8,   5,11,   4,10,  7,13,   6,12]);
    cats          = categories(dt.SliceCats);
    catcolors     = unique(dt.SliceCatsRGB);
    catcolors     = catcolors([1,1,2,2,3,3],:);
    catcolors.cat = categorical(cats);
    inBarColor = {'w','w','k','k','k','k'};
else
    replResSorted = replRes(:,[1, 14,15, 5,11, 6,12, 7,13, 2,8, 3,9, 4,10]);
    cats          = categories(dt.SliceCats);
    catcolors     = unique(dt.SliceCatsRGB);
    catcolors.cat = categorical(cats);
    catcolors     = [table(0,0,0,{'WHLorig'},'VariableNames',{'R','G','B','cat'}); ...
                     catcolors];
    catcolors     = catcolors([1,5,6,7,2,3,4],:);
    inBarColor = {'w','k','w','k','w','k','k'};
end
% CREATE FIGURE AND PLOT
bigfig = figure('Name',fnameRoot, ...
                'NumberTitle','off', ...
                'visible',   'on', ...
                'color','w', ...
                'Units','pixel', ...
                'Position',[0 0 1900 1100]);
for nt = 1: height(replResSorted)
    set(bigfig,'defaultAxesColorOrder',[left_color; right_color]);
    subplot(3,4,nt)
    tn = string(replResSorted{nt,'Struct'});
    Names = replResSorted.Properties.VariableNames(contains(replResSorted.Properties.VariableNames, 'Mean_'));
    Names = strrep(Names,'Mean_','');
    means = replResSorted{nt,contains(replResSorted.Properties.VariableNames, 'Mean_')};
    stds  = replResSorted{nt,contains(replResSorted.Properties.VariableNames, 'SD_')};
    CV    = 100*(stds ./ means);
    if meanOrSd=="mean"
        yyaxis right; ylim([0 20]); ylabel('CoV (%)');
        plot([1:length(CV)],CV,'--','Color',right_color,'linewidth',2);
    end
    for nc=1:length(means)
        yyaxis left; 
        if meanOrSd=="SD"
            ylabel('SD of FA'); 
            bar(nc, stds(nc), 'FaceColor', catcolors{nc,{'R','G','B'}},'EdgeColor','none','BarWidth',.9); hold on;
            ylim([0, .05])
        else
            ylabel('Mean FA'); 
            bar(nc, means(nc), 'FaceColor', catcolors{nc,{'R','G','B'}},'EdgeColor','none','BarWidth',.9); hold on;
            ylim([0, .8])   
            d  = errorbar(nc,means(nc),stds(nc),'.', 'color',[.5 .5 .5 ]); d.LineWidth = 2;
            yyaxis right; ylim([0 20]);ylabel('CoV (%)');
            plot(nc,CV(nc),'o', 'MarkerEdgeColor',right_color,'MarkerFaceColor','w','MarkerSize',8);hold on;
        end
        set(gca,'xtick',[]); set(gca,'xtickLabel',[]);
        T1 = text(nc, -0.05, string(catcolors.cat(nc)),'Color',[.0 .0 .0 ], 'Rotation', 45, 'HorizontalAlignment', 'Right');
        end
    set(gca,'FontSize',14)
    title(sprintf('%s',tn))
end
if HCPTRT
    suptitle({strrep(fnameRoot,'_','\_'), 'HCP TEST-RETEST'})
else
    suptitle({strrep(fnameRoot,'_','\_')})
end

%% Save it
if ~exist(saveItHere,'dir')
    mkdir(saveItHere)
end
if saveSvg
    saveas(gcf,fullfile(saveItHere, strcat(fnameRoot,'.svg')),'svg');
end
if savePng
    saveas(gcf,fullfile(saveItHere, strcat(fnameRoot,'.png')),'png');
end
