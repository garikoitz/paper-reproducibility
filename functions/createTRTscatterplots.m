function [TRTtests,allCOVS] = createTRTscatterplots(dt,unsMeans,varargin)
%CREATEPROFILE Summary of this function goes here
%   Detailed explanation goes here
%
% 
% Syntax:
%     createTRTscatterplots(dt)
%
% Description:
%  Input a table with profiles and it will plot them
%
% Inputs: (required)
%  dt: datatable
%  unsMeans: datatable
% 
% Optionals: 
% fnameRoot : string
% saveItHere: string
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
addRequired(p, 'unsMeans');

addOptional(p, 'fnameRoot'   , "changeThisName" , @isstring);
addOptional(p, 'saveItHere'  , "~/tmp"          , @isstring);
addOptional(p, 'savePng'     , false            , @islogical);
addOptional(p, 'saveSvg'     , false            , @islogical);
addOptional(p, 'WahlOrder'   , false            , @islogical);
addOptional(p, 'HCPTRT'      , false            , @islogical);
parse(p,dt,unsMeans,varargin{:});

fnameRoot   = p.Results.fnameRoot;
saveItHere  = p.Results.saveItHere;
savePng     = p.Results.savePng;
saveSvg     = p.Results.saveSvg;
WahlOrder   = p.Results.WahlOrder;
HCPTRT      = p.Results.HCPTRT;

%% Prepare the data
tracts    = unique(dt.Struct);
% Wahl Order
if WahlOrder
    newOrdering = [3,4,  11,12,  5,6,   7,8,   9,10,   1,2];
    tracts    = tracts(newOrdering);
end


omean   = @(x) mean(x,'omitnan');
ostd    = @(x) std(x,'omitnan');
left_color = [0 0 .55]; right_color = [0.3 0.3 0.3];


tmp     = varfun(ostd,dt,'GroupingVariables',{'Struct','SliceCats'},'InputVariables',{'meanVal'});
replRes = varfun(omean,dt,'GroupingVariables',{'Struct','SliceCats'},'InputVariables',{'meanVal'});
replRes.Properties.VariableNames{'Fun_meanVal'} = 'Mean';
replRes.SD  = tmp.Fun_meanVal;
replRes     = replRes(:,[1,2,4,5]);
replRes     = unstack(replRes,{'Mean','SD'},{'SliceCats'});
replRes     = replRes([3,4,  11,12,   5,6,   7,8,   9,10,   1,2],:);
replResSortedTRT      = replRes(:,[1, 3,9,   2,8,   5,11,   4,10,  7,13,   6,12]);


figure('Name',fnameRoot, ...
                'NumberTitle','off', 'visible',   'on', ...
                'color','w', 'Units','pixel', ...
                'Position',[0 0 1900 1100]);
allrmse1000  = []; allrmse2000  = []; allrmse3000  = [];
allrrmse1000 = []; allrrmse2000 = []; allrrmse3000 = [];
allicc1000   = []; allicc2000   = []; allicc3000   = [];
allcov1000   = []; allcov2000   = []; allcov3000   = [];
Tstd1000     = []; Tstd2000     = []; Tstd3000     = [];
RTstd1000    = []; RTstd2000    = []; RTstd3000    = [];
for nt = 1: height(replResSortedTRT)
    subplot(3,4,nt)
    tn = string(replResSortedTRT{nt,'Struct'});
    b1000T    = unsMeans(unsMeans.TRT=='TEST' & unsMeans.SHELL=='1000',{char(tn),'SubjID'});
    b1000RT   = unsMeans(unsMeans.TRT=='RETEST' & unsMeans.SHELL=='1000', {char(tn),'SubjID'});
    b1000T    = sortrows(b1000T,'SubjID'); b1000RT=sortrows(b1000RT,'SubjID');
    color1000 = unique(unsMeans{unsMeans.TRT=='TEST' & unsMeans.SHELL=='1000','SliceCatsRGB'});
    if isequal(b1000T.SubjID,b1000RT.SubjID)
        b1000T = b1000T.(tn); b1000RT = b1000RT.(tn);
    else
        error('The test-retest is not comparing the same subjects')
    end

    b2000T    = unsMeans(unsMeans.TRT=='TEST' & unsMeans.SHELL=='2000', {char(tn),'SubjID'});
    b2000RT   = unsMeans(unsMeans.TRT=='RETEST' & unsMeans.SHELL=='2000', {char(tn),'SubjID'});
    b2000T=sortrows(b2000T,'SubjID'); b2000RT=sortrows(b2000RT,'SubjID');
    color2000 = unique(unsMeans{unsMeans.TRT=='TEST' & unsMeans.SHELL=='2000','SliceCatsRGB'});
    if isequal(b2000T.SubjID,b2000RT.SubjID)
        b2000T = b2000T.(tn); b2000RT = b2000RT.(tn);
    else
        error('The test-retest is not comparing the same subjects')
    end
    
    b3000T    = unsMeans(unsMeans.TRT=='TEST' & unsMeans.SHELL=='3000', {char(tn),'SubjID'});
    b3000RT   = unsMeans(unsMeans.TRT=='RETEST' & unsMeans.SHELL=='3000', {char(tn),'SubjID'});
    b3000T=sortrows(b3000T,'SubjID'); b3000RT=sortrows(b3000RT,'SubjID');
    color3000 = unique(unsMeans{unsMeans.TRT=='TEST' & unsMeans.SHELL=='3000','SliceCatsRGB'});
    if isequal(b3000T.SubjID,b3000RT.SubjID)
        b3000T = b3000T.(tn); b3000RT = b3000RT.(tn);
    else
        error('The test-retest is not comparing the same subjects')
    end
    
    scatter(b1000T,b1000RT,22,color1000{:,:}, 'filled');hold on;
    scatter(b2000T,b2000RT,22,color2000{:,:}, 'filled');
    scatter(b3000T,b3000RT,22,color3000{:,:}, 'filled');
    % xlim([0.3, 0.66]);ylim([0.3, 0.66]);
    % xticks([0.3:0.1:0.6]);yticks([0.3:0.1:0.6]);
    
    
    identityLine(gca);
    axis equal
    
    xlim([0.3, 0.66]);ylim([0.3, 0.66]);
    xticks([0.4:0.1:0.6]);yticks([0.4:0.1:0.6]);
    location='southeast';if nt >=11;location='southeast';end
    [N1, Nold1, Xm1,Ym1,rho1,pval1,rhom1,pvalm1,rmse1,rmsem1,rrmse1,rrmsem1,sdX1,sdY1,sdXm1,sdYm1,icc1,iccm1,CoV1,CoVm1] = dr_corrrmse(b1000T,b1000RT);
    [N2, Nold2, Xm2,Ym2,rho2,pval2,rhom2,pvalm2,rmse2,rmsem2,rrmse2,rrmsem2,sdX2,sdY2,sdXm2,sdYm2,icc2,iccm2,CoV2,CoVm2] = dr_corrrmse(b2000T,b2000RT);
    [N3, Nold3, Xm3,Ym3,rho3,pval3,rhom3,pvalm3,rmse3,rmsem3,rrmse3,rrmsem3,sdX3,sdY3,sdXm3,sdYm3,icc3,iccm3,CoV3,CoVm3] = dr_corrrmse(b3000T,b3000RT);
%     legend({sprintf('b=1000s/mm^2 (rmse=%.3f)',rmse1),...
%             sprintf('b=2000s/mm^2 (rmse=%.3f)',rmse2),...
%             sprintf('b=3000s/mm^2 (rmse=%.3f)',rmse3)}, ...
%             'location',location, 'Box','off')
    allrmse1000  = [allrmse1000 , rmse1 ];allrmse2000  = [allrmse2000 , rmse2 ];allrmse3000  = [allrmse3000 , rmse3 ];
    allrrmse1000 = [allrrmse1000, rrmse1];allrrmse2000 = [allrrmse2000, rrmse2];allrrmse3000 = [allrrmse3000, rrmse3];
    allicc1000   = [allicc1000  , icc1  ];allicc2000   = [allicc2000  , icc2  ];allicc3000   = [allicc3000  , icc3  ];
    allcov1000   = [allcov1000  , CoV1  ];allcov2000   = [allcov2000  , CoV2  ];allcov3000   = [allcov3000  , CoV3  ];
    Tstd1000     = [Tstd1000    , sdXm1 ];Tstd2000     = [Tstd2000  , sdXm2   ];Tstd3000     = [Tstd3000    , sdXm3 ];
    RTstd1000    = [RTstd1000   , sdYm1 ];RTstd2000    = [RTstd2000  , sdYm2  ];RTstd3000    = [RTstd3000   , sdYm3 ];
    title(sprintf('%s',tn));
%     T1 = text(.25, 0.625, sprintf('CoV (b=1000)=%1.2f%%',CoVm1),'Color',color1000{:,:}, 'Rotation', 0, 'HorizontalAlignment', 'Left','fontweight','normal','fontsize',16);
%     T2 = text(.20, 0.575, sprintf('CoV (b=2000)=%1.2f%%',CoVm2),'Color',color2000{:,:}, 'Rotation', 0, 'HorizontalAlignment', 'Left','fontweight','normal','fontsize',16);
%     T3 = text(.15, 0.525, sprintf('CoV (b=3000)=%1.2f%%',CoVm3),'Color',color3000{:,:}, 'Rotation', 0, 'HorizontalAlignment', 'Left','fontweight','normal','fontsize',16);
    T1 = text(.55, 0.625, sprintf('%1.2f',rmsem1),'Color',color1000{:,:}, 'Rotation', 0, 'HorizontalAlignment', 'center','fontweight','normal','fontsize',16);
    T2 = text(.40, 0.575, sprintf('%1.2f',rmsem2),'Color',color2000{:,:}, 'Rotation', 0, 'HorizontalAlignment', 'center','fontweight','normal','fontsize',16);
    T3 = text(.35, 0.525, sprintf('%1.2f',rmsem3),'Color',color3000{:,:}, 'Rotation', 0, 'HorizontalAlignment', 'center','fontweight','normal','fontsize',16);
    xlabel('MEAN FA (TEST)','FontWeight','bold');
    ylabel('MEAN FA (RETEST)','FontWeight','bold');
    set(gca,'FontSize',18)
    grid off

end

if HCPTRT
    suptitle({strrep(fnameRoot,'_','\_'), 'HCP TEST-RETEST'})
else
    suptitle({strrep(fnameRoot,'_','\_')})
end




%% Create output tables
RMSE     = [mean(allrmse1000) , mean(allrmse2000) ,mean(allrmse3000,'omitnan')]';
rRMSE    = [mean(allrrmse1000), mean(allrrmse2000),mean(allrrmse3000,'omitnan')]';
ICC      = [mean(allicc1000)  , mean(allicc2000)  ,mean(allicc3000,'omitnan')]';
CoV      = [mean(allcov1000)  , mean(allcov2000)  ,mean(allcov3000,'omitnan')]';
minCOV   = [min(allcov1000)   , min(allcov2000)   ,min(allcov3000)]';
maxCOV   = [max(allcov1000)   , max(allcov2000)   ,max(allcov3000)]';
TRTtests = table(RMSE,rRMSE,ICC,CoV,minCOV,maxCOV);
TRTtests.Properties.RowNames = {'b1000','b2000','b3000'};

allCOVS = replResSortedTRT(:,'Struct');
allCOVS.TRT_b1000_CoV = allcov1000';
allCOVS.TRT_b2000_CoV = allcov2000';
allCOVS.TRT_b3000_CoV = allcov3000';


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
