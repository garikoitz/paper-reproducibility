function createDistributionsPlots(dt,unsMeans,varargin)
%CREATEPROFILE Summary of this function goes here
%   Detailed explanation goes here
%
% 
% Syntax:
%     createDistributionsPlots(dt)
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


nBoot     = 5000;
normalKsdensity = 'normal';
ylims = [0, 20];
xlims = [0.1, .8];
omean   = @(x) mean(x,'omitnan');
ostd    = @(x) std(x,'omitnan');

tmp     = varfun(ostd,dt,'GroupingVariables',{'Struct','SliceCats'},'InputVariables',{'meanVal'});
replRes = varfun(omean,dt,'GroupingVariables',{'Struct','SliceCats'},'InputVariables',{'meanVal'});
replRes.Properties.VariableNames{'Fun_meanVal'} = 'Mean';
replRes.SD  = tmp.Fun_meanVal;
replRes     = replRes(:,[1,2,4,5]);
% Unstack it to show it as consecutive tables. Reorder to replicate
replRes     = unstack(replRes,{'Mean','SD'},{'SliceCats'});
replRes     = replRes([3,4,  11,12,   5,6,   7,8,   9,10,   1,2],:);
replRes.Mean_WHLorig = [0.565,0.522,0.517,0.490,0.549,0.534,0.510,0.497,0.489,0.470,0.587,0.575]';
replRes.SD_WHLorig   = [0.028,0.028,0.023,0.03,0.026,0.024,0.028,0.026,0.022,0.024,0.023,0.021]';
% replResSorted        = replRes(:,[1, 14,15,  5,11,   6,12,   7,13,   2,8,  3,9,   4,10]);
replResSorted = replRes;


if ~HCPTRT

    figure('Name',fnameRoot, ...
                    'NumberTitle','off', ...
                    'visible',   'on', ...
                    'color','w', ...
                    'Units','pixel', ...
                    'Position',[0 0 1900 1100]);
    for nt = 1: height(replResSorted)
        subplot(3,4,nt)
        tn = char(replResSorted{nt,'Struct'});
        X       = unsMeans.(tn);
        [X_values, MU, SIGMA, MN, MX] = dr_distPlottingVals(X);
        if strcmp(normalKsdensity, 'normal')
            PD = fitdist(X,'normal');
            group_pdf = pdf(PD, X_values);
            mySupTitle = 'Normal distribution of mean FA-s per project/tract.';
        else
            [group_pdf, X_values] = ksdensity(X);
            mySupTitle = 'Density distribution of mean FA-s per project/tract.';
        end
        a = plot(X_values,group_pdf,'Color','k','LineStyle',':','LineWidth',3); hold on;
        title(tn); ylim(ylims);xlim(xlims);
        [ph,msg]=jbfill(X_values,group_pdf,zeros(size(group_pdf)),'k','k',0,0.4);
        h1=plot(MU*[1,1],[0 max(group_pdf)],'Color','k','LineStyle',':','LineWidth',1);
        xlabel('FA', 'FontWeight','bold'); 
        % ylim([0.1, 0.75]); yticks([0.2,.4,.6])
        set(gca,'FontSize',18)
        title(sprintf('%s',tn))
        cats      = categories(unsMeans.SliceCats);
        catcolors = unique(unsMeans.SliceCatsRGB);
        for sc=1:length(cats)
            x = unsMeans.(tn)(unsMeans.SliceCats==cats{sc});
            [x_values, mu, sigma, mn, mx] = dr_distPlottingVals(x);
            if strcmp(normalKsdensity, 'normal')
                pd = fitdist(x,'normal');
                indiv_pdf = pdf(pd, x_values);
            else
                [indiv_pdf, x_values] = ksdensity(x);
            end
            a = [a;plot(x_values, indiv_pdf/length(cats), 'LineWidth',2, 'color', catcolors{sc,:})]; 
            % Now calculate the effect size with bootstrapping and mark it in plot
            stats = mes( X, x, ...
                         'hedgesg','isDep',0,'nBoot',nBoot, ...
                         'doPlot',0,'missVal','listwise','confLevel',.95,'ROCtBoot',false);
            if stats.t.p > 0.05; texto = 'n.s.';
            else;                texto = sprintf('d=%6.3f',stats.hedgesg); end
            hArrw = (length(cats)-sc+1)*(max(group_pdf)/length(cats));
            h1=plot(mu*[1,1],[0 hArrw],'Color',catcolors{sc,:},'LineStyle','-','LineWidth',.7);
            drawArrow([mu MU],[hArrw,hArrw],{'Color',catcolors{sc,:},'LineWidth',.7, ...
                      'string',texto, 'FontSize', 14});
        end
        % legend(a,[{'All Samples'};strrep(cats,'_','\_')],'Location','NorthEast'); hold off;                             

    end
else

    % {
    % Same plot for TEST-RETEST
    HCPunsMeans=unsMeans(unsMeans.Proj=='HCP' & unsMeans.TRT=='TEST',:);
    HCPunsMeans.SliceCats = removecats(HCPunsMeans.SliceCats);
    HCPunsMeansrt=unsMeans(unsMeans.Proj=='HCP' & unsMeans.TRT=='RETEST',:);
    HCPunsMeansrt.SliceCats = removecats(HCPunsMeansrt.SliceCats);

    cats      = categories(HCPunsMeans.SliceCats);
    catcolors = unique(HCPunsMeansrt.SliceCatsRGB);
    catsrt      = categories(HCPunsMeansrt.SliceCats);
    catcolorsrt = catcolors;
    
    
    
    nBoot     = 5000;
    normalKsdensity = 'normal';
    figure('Name',fnameRoot, ...
                    'NumberTitle','off', ...
                    'visible',   'on', ...
                    'color','w', ...
                    'Units','pixel', ...
                    'Position',[0 0 1900 1100]);
    for nt = 1: height(replResSorted)
        subplot(3,4,nt)
        tn = char(replResSorted{nt,'Struct'});
        
        a = [];
        l = {};
        for sc=1:length(cats)
            x   = HCPunsMeans.(tn)(HCPunsMeans.SliceCats==cats{sc});
            xrt = HCPunsMeansrt.(tn)(HCPunsMeansrt.SliceCats==catsrt{sc});
            [x_values, mu, sigma, mn, mx] = dr_distPlottingVals(x);
            [x_valuesrt, murt, sigmart, mnrt, mxrt] = dr_distPlottingVals(xrt);
            if strcmp(normalKsdensity, 'normal')
                pd = fitdist(x,'normal');
                indiv_pdf = pdf(pd, x_values);
                pdrt = fitdist(xrt,'normal');
                indiv_pdfrt = pdf(pdrt, x_valuesrt);
                mySupTitle = 'Normal distribution of mean FA-s per project/tract.';
            else
                [indiv_pdf, x_values] = ksdensity(x);
                [indiv_pdfrt, x_valuesrt] = ksdensity(xrt);
                mySupTitle = 'Density distribution of mean FA-s per project/tract.';
            end
            at  = plot(x_values, indiv_pdf, 'LineWidth',2, 'color',catcolors{sc,:});hold on;
            art = plot(x_valuesrt, indiv_pdfrt, 'LineWidth',2, 'color', catcolors{sc,:}, 'LineStyle',':');
            xlim([0, .8]);ylim([0, 20]);
            xticks([0:0.2:.8]);yticks([0:5:20]);
            title(tn)
            a = [a;at;art];  % Concatenate plots to be used with legend
            l = [l; strcat(strrep(cats{sc},'_','\_'),'TEST'); ...
                    strcat(strrep(cats{sc},'_','\_'),'RETEST')];
            % Now calculate the effect size with bootstrapping and mark it in plot
    %         stats = mes( x, xrt, ...
    %                      'hedgesg','isDep',0,'nBoot',nBoot, ...
    %                      'doPlot',0,'missVal','listwise','confLevel',.95,'ROCtBoot',false);
    %         if stats.t.p > 0.05; sig = 'n.s.';
    %         else;                sig = sprintf('d=%6.2f',stats.t.p); end
    %         texto = sprintf('d=%6.3f (%s)',stats.hedgesg,sig);
    %         hArrw = max(indiv_pdf);
    %         hArrwrt = max(indiv_pdfrt);
    %         h1=plot(mu*[1,1],[0 hArrw],'Color',catcolors{sc,:},'LineStyle','-','LineWidth',1);
    %            plot(murt*[1,1],[0 hArrwrt],'Color',catcolors{sc,:},'LineStyle',':','LineWidth',1);
    %         drawArrow([mu murt],[hArrw,hArrw],{'Color',catcolors{sc,:},'LineWidth',1,'string',texto});
        end
        xlabel('FA', 'FontWeight','bold'); 
        set(gca,'FontSize',18)
        title(sprintf('%s',tn))

        % legend(a, l,'Location','SouthWest'); hold off;
    end
end
%}

if HCPTRT
    suptitle({strrep(mySupTitle,'_','\_'), 'HCP TEST-RETEST'})
else
    suptitle({strrep(mySupTitle,'_','\_')})
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
