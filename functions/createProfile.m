function createProfile(dt,varargin)
%CREATEPROFILE Summary of this function goes here
%   Detailed explanation goes here
%
% 
% Syntax:
%     dr_fwDownloadFileFromZip(dt)
%
% Description:
%  Input a table with profiles and it will plot them
%
% Inputs: (required)
%  dt: datatable
% 
% Optionals: 
% fnameRoot : char
% saveItHere: char
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
addOptional(p, 'savePng'     , false            , @islogical);
addOptional(p, 'saveSvg'     , false            , @islogical);
addOptional(p, 'WahlOrder'   , false            , @islogical);
addOptional(p, 'HCPTRT'      , false            , @islogical);
parse(p,dt,varargin{:});

fnameRoot   = p.Results.fnameRoot;
saveItHere  = p.Results.saveItHere;
savePng     = p.Results.savePng;
saveSvg     = p.Results.saveSvg;
WahlOrder   = p.Results.WahlOrder;
HCPTRT      = p.Results.HCPTRT;

%% 1.- Create the profile
cats      = categories(dt.SliceCats);
catcolors = unique(dt.SliceCatsRGB);
tracts    = unique(dt.Struct);
% Wahl Order
if WahlOrder
    tracts    = tracts([3,4,  11,12,  5,6,   7,8,   9,10,   1,2]);
end
% HCPTRT
if HCPTRT
    catcolors = catcolors([1,1,2,2,3,3],:);
    linestyles= {':','-',':','-',':','-'};
end


% Create full screen window conrolling size in pixels 
figHdl = figure('Name',fnameRoot, ...
                'NumberTitle','off', ...
                'visible',   'on', ...
                'color','w', ...
                'Units','pixel', ...
                'Position',[0 0 1900 1100]);
for nt = 1: length(tracts)
    subplot(3,4,nt)
    tn = tracts(nt);
    a = [];
    for nc=1:length(cats)
        cat      = cats{nc}; 
        cat      = cats{nc}; 
        if HCPTRT
            catColor  = catcolors{nc,:}; 
            lineStyle = linestyles{nc}; 
        else
            catColors = dt{dt.SliceCats==cat,'SliceCatsRGB'}; 
            catColor  = catColors{1,:};
            lineStyle = '-';
        end
        values   = mean(dt{dt.Struct==tn & dt.SliceCats==cat,'Val'},1,'omitnan');
        sdvalues = std(dt{dt.Struct==tn & dt.SliceCats==cat,'Val'},1,'omitnan');
        upper    = values + sdvalues; lower = values - sdvalues;
        a = [a; plot([1:100],values,'color',catColor,'linewidth',2,'LineStyle',lineStyle)]; 
        hold on;
        jbfill([1:100],upper,lower,catColor,catColor,0,0.1);
    end
    if WahlOrder
        if (nt > 8 && nt < 11)
            legend(a, strrep(cats,'_','\_'), 'location','northwest'); 
        else
            legend(a, strrep(cats,'_','\_'), 'location','southwest'); 
        end
    end
    ylabel('FA', 'FontWeight','bold'); % xlabel('Divisions')
    set(gca,'xtick',[])
    ylim([0.1, 0.75]); yticks([0.2,.4,.6])
    set(gca,'FontSize',18)
    title(sprintf('%s',tn))
end
suptitle({strrep(fnameRoot,'_','\_'), 'Group average (1 SD)'})

if ~exist(saveItHere,'dir')
    mkdir(saveItHere)
end

if saveSvg
    saveas(gcf,fullfile(saveItHere, strcat(fnameRoot,'.svg')),'svg');
end
if savePng
    saveas(gcf,fullfile(saveItHere, strcat(fnameRoot,'.png')),'png');
end


