%% INSTRUCTIONS


% Download or clone the repository
%       !git clone https://github.com/garikoitz/paper-reproducibility.git
% 
clear all; close all; clc;
% Add the root of the repository to the Matlab path (or run this code):

    % cd('<path-to-your-code>/paper-reproducibility')
    cd('~/soft/paper-reproducibility')
    rootDir = pwd;
    addpath(genpath(rootDir));

% Specify a path to save the output figures. 
paperPath = '~/gDrive/STANFORD/PROJECTS/2018 Computational reproducibility (Gari-Pratik-Zhimei-Brian)';
saveItHere = fullfile(paperPath, '01_REVIEW01/Figures/REV01_v01/sources');


% Read the data 
% Check if there is a local cache, otherwise download it from FW
DataVersion    = '03';
collectionName = 'ComputationalReproducibility';
measure        = 'fa';

fname          = sprintf('AllV%s_%s_%s.mat',DataVersion, collectionName, measure);
localfname = fullfile(paperReprPath,'local',fname);
if exist(localfname,'file')
    data = load(localfname);
else  % Download it from the Flywheel collection attachment
    serverName     = 'stanfordlabs';
    st  = scitran(serverName);
    cc  = st.search('collection','collection label contains',collectionName);
    data=load(st.fw.downloadFileFromCollection(cc{1}.collection.id,fname,localfname));
end

% This script can be run directly using the Run button or step by step. 

%% Checks and conversions
DT = data.dt;
% Rename project names coming from FW
DT.Proj = renamecats(DT.Proj,'HCP_preproc','HCP');
DT.Proj = renamecats(DT.Proj,'PRATIK','WHL');
DT.Proj = renamecats(DT.Proj,'Weston Havens','YWM');
% Change the names of the TRT fields
DT.TRT(DT.Proj=="WHL") = 'TEST';
DT.TRT(DT.Proj=="YWM")   = 'TEST';
DT.TRT = removecats(DT.TRT);

okk = @(x)  x;
unique(varfun(okk,DT,'GroupingVariables',{'Proj','SubjID','TRT'},'InputVariables',{'TRT'}))


summary(DT)
cmapname = 'copper';
if ~exist(saveItHere); mkdir(saveItHere); end

%% Create working tables
% Select only the Wahl tracts from the AFQ output
useTracts = [{'CingulumCingulate'}, {'Arcuate'}, {'IFOF'}, {'ILF'}, {'Uncinate'}, ...
          {'Corticospinal'}];

% Do the filtering and obtain the unstacked tract profiles
[dtt, unsProf, unsMeans, indivTracts, pairedTracts] = dr_createWorkingTables(DT, ...
                                 'sliceBasedOn',           {'Proj', 'SHELL'}, ...
                                 'SHELL',                  {'1000', '2000', '3000'}, ...
                                 'removeTractsContaining', {'Callosum'}, ... 
                                 'useThisTractsOnly',       useTracts, ...
                                 'createTractPairs',        true, ... 
                                 'TRT',                    'TEST', ...
                                 'AGE',                     [14,58], ...
                                 'GENDER',                 {'male', 'female'});                             
                             
[dtTRT, unsProfTRT, unsMeansTRT] = dr_createWorkingTables(DT(DT.Proj=='HCP',:), ...
                                 'sliceBasedOn',           {'Proj', 'SHELL','TRT'}, ...
                                 'SHELL',                  {'1000', '2000', '3000'}, ...
                                 'removeTractsContaining', {'Callosum'}, ... 
                                 'useThisTractsOnly',       useTracts, ...
                                 'createTractPairs',        true, ...
                                 'TRT',                    'notTRAIN', ...
                                 'AGE',                     [14,58], ...
                                 'GENDER',                 {'male', 'female'});
                             
[dtt1000, unsProft1000, unsMeanst1000] = dr_createWorkingTables(DT, ...
                                 'sliceBasedOn',           {'Proj', 'SHELL'}, ...
                                 'SHELL',                  {'1000'}, ...
                                 'removeTractsContaining', {'Callosum'}, ... 
                                 'useThisTractsOnly',       useTracts, ...
                                 'createTractPairs',        true, ...
                                 'TRT',                    'TEST', ...
                                 'AGE',                     [14,58], ...
                                 'GENDER',                 {'male', 'female'});
                             
[dtt2000, unsProft2000, unsMeanst2000] = dr_createWorkingTables(DT, ...
                                 'sliceBasedOn',           {'Proj', 'SHELL'}, ...
                                 'SHELL',                  {'2000'}, ...
                                 'removeTractsContaining', {'Callosum'}, ... 
                                 'useThisTractsOnly',       useTracts, ...
                                 'createTractPairs',        true, ...
                                 'TRT',                    'TEST', ...
                                 'AGE',                     [14,58], ...
                                 'GENDER',                 {'male', 'female'});
                             
                             
% See what's going on with the nans, if the pipeline worked there shuold be
% almost no NaNs (NaNs mean that the tract was not found)
fprintf('Total number of NaN: %d, total values: %d, %% of NaNs: %2.2f%%\n', sum(isnan(dtt.meanVal)), length(dtt.meanVal), 100*(sum(isnan(dtt.meanVal))/length(dtt.meanVal)))
fprintf('\n\nSame thing but per project: \n')
onans = @(x) nnz(isnan(x));
nansByProject = varfun(onans,dtt,'GroupingVariables',{'Proj','SHELL'},'InputVariables',{'meanVal'});
nansByProject.Percentage = 100*(nansByProject.Fun_meanVal./nansByProject.GroupCount)
fprintf('\n\nSame thing but per tract and project: \n')
nansByTractProject = varfun(onans,dtt,'GroupingVariables',{'Proj','SHELL','Struct'},'InputVariables',{'meanVal'});
nansByTractProject.Percentage = 100*(nansByTractProject.Fun_meanVal./nansByTractProject.GroupCount)
                       
%% Look at the populations for the three projects
AGE_GENDER_TABLE = unsProf(unsProf.SHELL=='1000', {'Proj','AGE','GENDER'});
AGE_means        = varfun(@mean ,AGE_GENDER_TABLE,'GroupingVariables','Proj','InputVariables',{'AGE'});
AGE_std          = varfun(@std  ,AGE_GENDER_TABLE,'GroupingVariables','Proj','InputVariables',{'AGE'});
AGE_means.std_AGE = AGE_std.std_AGE;
GENDER_mean      = varfun(@mean,AGE_GENDER_TABLE,'GroupingVariables',{'Proj','GENDER'},'InputVariables',{'AGE'});
GENDER_std       = varfun(@std,AGE_GENDER_TABLE,'GroupingVariables',{'Proj','GENDER'},'InputVariables',{'AGE'});
GENDER_mean.std_GENDER = GENDER_std.std_AGE;
% ANOVA
[P,ANOVATAB,STATS] = anova1(AGE_GENDER_TABLE.AGE, AGE_GENDER_TABLE.Proj);
saveas(gcf,fullfile(saveItHere, 'AGE_ANOVA_a.png'),'png');close(gcf);close(gcf);
% ANOVATAB
% Source     SS      df      MS       F    Prob>F
% -----------------------------------------------
% Groups      36.8     2   18.3939   0.2   0.8209
% Error    12007.9   129   93.0847               
% Total    12044.7   131

%% CREATE FA PROFILES 

% REPLICATION EXPERIMENT: HCP TEST-RETEST (FIGURE 4A)
fnameRoot = "FA_Profiles_HCP_TEST_RETEST";
createProfile(dtTRT,'fnameRoot',fnameRoot,'saveItHere',string(saveItHere), ...
                    'saveSvg'  ,true     ,'WahlOrder' ,true, 'HCPTRT',true)

% GENERALIZATION EXPERIMENT: b1000 and b2000 (FIGURE 5A)
bvals     = [1000   , 2000];
dts       = {dtt1000, dtt2000};
% Calculate for b1000 and b2000
for nb=1:length(bvals)
    bval      = bvals(nb);
    fnameRoot = strcat("FA_Profiles_b", num2str(bval));
    createProfile(dts{nb},'saveSvg',true,'saveItHere',string(saveItHere),'fnameRoot',fnameRoot,'WahlOrder',true)
end

%% TEST: MEAN FA BARS/GENERALIZATION PLOTS
% Plot bars with covar
fnameRoot = "FA_Means_withCoVdots";
createBars(dtt,'fnameRoot',fnameRoot,'saveItHere',string(saveItHere), ...
                    'saveSvg'  ,true     ,'WahlOrder' ,true, 'HCPTRT',false)




% ALL: MEAN FA GI
% Prepare the data to be able to plot later on. 
% This has bee copied from below, were the correlations were calculated
recalculate  = true;
nRep = 500;
CIrangeOrVals = 90;  % [1:100];
useDistribution = true;


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

for CIrange=CIrangeOrVals
% OBTAIN CORR AND CIs
    for ns = 1:3
        [meanTable, lowerCI, upperCI, bootstrapStat] = dr_meanBoots(...
                               unsMeans(unsMeans.SliceCats==string(cats(ns)),:), ... 
                               'tractsOrder',tractsOrder,'newTractNames',WahlTractNames, ...
                               'useDistribution',useDistribution,'nRep',nRep,'CIrange',CIrange);
        meanTables{ns}        = meanTable;
        corrLowers{ns}        = lowerCI;
        corrUppers{ns}        = upperCI;
        corrBootstrapStat{ns} = bootstrapStat;   
        % Do HCP-T
        [meanTable, lowerCI, upperCI, bootstrapStat] = dr_meanBoots(...
                               unsMeans(unsMeans.SliceCats==string(cats(ns+3)),:), ... 
                               'tractsOrder',tractsOrder,'newTractNames',WahlTractNames, ...
                               'useDistribution',useDistribution,'nRep',nRep,'CIrange',CIrange);
        meanTables{ns+3}   = meanTable;
        corrLowers{ns+3}   = lowerCI;
        corrUppers{ns+3}   = upperCI;
        corrBootstrapStat{ns+3} = bootstrapStat;
        % Do HCP-RT
        [meanTable, lowerCI, upperCI, bootstrapStat] = dr_meanBoots(...
                               unsMeansrt(unsMeansrt.SliceCats==string(cats(ns+3)),:), ... 
                               'tractsOrder',tractsOrder,'newTractNames',WahlTractNames, ...
                               'useDistribution',useDistribution,'nRep',nRep,'CIrange',CIrange);
        meanTables{ns+6} = meanTable;
        corrLowers{ns+6}   = lowerCI;
        corrUppers{ns+6}   = upperCI;
        corrBootstrapStat{ns+6} = bootstrapStat;
    end

    % Add the results from the paper (do not do it here but maintain structure)
    AllMeans        = [meanTables]';
    AllLower        = [corrLowers]';
    AllUpper        = [corrUppers]';
    AllBStats       = [corrBootstrapStat]';

    Experiments     = [cellstr(cats(1:3));strcat(cellstr(cats(4:6)),'TEST');strcat(cellstr(cats(4:6)),'RETEST')];
    sumPlotTable    = table(Experiments,AllMeans,AllLower,AllUpper,AllBStats);

%     UnstackedMean       = {};
%     UnstackedLower      = {};
%     UnstackedUpper      = {};
%     UnstackedBStats     = {};
%     for nt=1:height(sumPlotTable)
%         UnstackedMean  = [UnstackedMean;  {sumPlotTable.AllMeans{nt}}];
%         UnstackedLower = [UnstackedLower; {sumPlotTable.AllLower{nt}}];
%         UnstackedUpper = [UnstackedUpper; {sumPlotTable.AllUpper{nt}}];
%         UnstackedBStats= [UnstackedBStats; {sumPlotTable.AllBStats{nt}}];
%     end
%     sumPlotTable.UnstackedMean  = UnstackedMean;
%     sumPlotTable.UnstackedLower = UnstackedLower;
%     sumPlotTable.UnstackedUpper = UnstackedUpper;
%     sumPlotTable.UnstackedBStats= UnstackedBStats;
%     noTable = sumPlotTable(:,~contains(sumPlotTable.Properties.VariableNames,'All'));
    noTable = sumPlotTable;
    
    
    % Now make the unstacked col numeric
    long       = noTable{1,'AllMeans'}{1};
    AllMeans   = long{:,:}';
    TractNames = long.Properties.VariableNames';
    long = table(TractNames,AllMeans);
    long.Properties.VariableNames = {'TractNames', noTable{1,1}{1}};
    for na = 2:height(noTable)
        tmp = noTable{na,'AllMeans'}{1};
        long.(noTable{na,1}{1}) = tmp{:,:}';
    end
    long.Type = cellstr(repmat('Mean', [height(long),1]));
    
    % Now make the unstacked col numeric
    longL = noTable{1,'AllLower'}{1};
    AllLower   = longL{:,:}';
    TractNames = longL.Properties.VariableNames';
    longL = table(TractNames,AllLower);
    longL.Properties.VariableNames = {'TractNames', noTable{1,1}{1}};
    for na = 2:height(noTable)
        tmp = noTable{na,'AllLower'}{1};
        longL.(noTable{na,1}{1}) = tmp{:,:}';
    end
    longL.Type = cellstr(repmat('Lower', [height(long),1]));

    % Now make the unstacked col numeric
    longU = noTable{1,'AllUpper'}{1};
    AllUpper   = longU{:,:}';
    TractNames = longU.Properties.VariableNames';
    longU = table(TractNames,AllUpper);
    longU.Properties.VariableNames = {'TractNames', noTable{1,1}{1}};
    for na = 2:height(noTable)
        tmp = noTable{na,'AllUpper'}{1};
        longU.(noTable{na,1}{1}) = tmp{:,:}';
    end
    longU.Type = cellstr(repmat('Upper', [height(long),1]));

    % Now make the unstacked col numeric
    longB = noTable{1,'AllBStats'}{1};
    AllBstats   = longB{:,:}';
    TractNames = longB.Properties.VariableNames';
    longB = table(TractNames,AllBstats);
    longB.Properties.VariableNames = {'TractNames', noTable{1,1}{1}};
    for na = 2:height(noTable)
        tmp = noTable{na,'AllBStats'}{1};
        longB.(noTable{na,1}{1}) = tmp{:,:}';
    end
    longB.Type = cellstr(repmat('BStats', [height(long),1]));

    % Aggregate
    long = [long; longL; longU]; 
    long.Type = categorical(long.Type);
    long.TractNames = categorical(long.TractNames);
    longcell{CIrange} = long;
end
% I did not create the longB for all the CIRange
% save('allCI1to100correlationsv2.mat','longcell')




% Instead of the bar plot use the generalization index
referenceColumn = 'NoReference';
versus          = 'NoRetest';
plotIt          = true;
plotIndex       = false;
showLegend = false; showXnames = false; refLine=true;
GIcoloringMethod = 'bars'; % 'circle', 'GIband', 'GRband', 'background','none','bars'
% takeout =  {'WHLorig','YWM2000', ...
%             'HCP2000TEST','HCP2000RETEST', ...
%             'HCP3000TEST','HCP3000RETEST'};
% takeout =  {'WHLorig'};        
takeout =  {'WHLorig','WHL1000', ...
            'YWM1000','YWM2000'};
takeout =  {'WHLorig', ...
            'HCP1000RETEST','HCP2000RETEST', ...
            'HCP3000RETEST'};
% takeout =  {'WHLorig','YWM2000', ...
%             'WHL1000','YWM1000', ...
%             'HCP1000TEST','HCP1000RETEST', ...
%             'HCP2000TEST','HCP2000RETEST'};
% takeout =  {'YWM1000','YWM2000', ...
%             'HCP1000TEST','HCP1000RETEST', ...
%             'HCP2000TEST','HCP2000RETEST', ...
%             'HCP3000TEST','HCP3000RETEST'};

bigfig = mrvNewGraphWin('Mean FA per project'); 
% figHdl = figure('Name','FA profiles', ...
%                 'NumberTitle','off', ...
%                 'visible',   'on', ...
%                 'color','w', ...
%                 'Units','pixel', ...
%                 'Position',[0 0 1900 1100]);
set(bigfig,'WindowStyle','normal')
set(bigfig, 'Units', 'inches', 'OuterPosition', [0, 0, 10, 14]);
ncol=3; nrow=4;
for nc = 1:length(TractNames)
    % Calculate the position
    tn  = string(long.TractNames(nc));
    sp = subplot(ncol,nrow,nc);
    % DO the calculation/plots
    [WahlInsideCI, isOverlap, minUpper, maxLower, midCI, GI]=dr_compareCI(long, longB, tn,...
                                                    'referenceColumn',referenceColumn,...
                                                    'refLine' , refLine,'plotIndex',plotIndex, ...
                                                    'plotIt',plotIt,'showLegend',showLegend, ...
                                                    'takeout',takeout, 'showXnames',showXnames, ...
                                                    'GIcoloringMethod', GIcoloringMethod, ...
                                                    'ResultType', 'FA');
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
set(bigfig,'color','w');
suptitle({sprintf('GI_FA_Meansvs%s.png',versus), 'FA distributions'})

saveas(bigfig,fullfile(saveItHere, sprintf('GI_withIndex_FA_Meansvs%s.png',versus)),'png');
saveas(bigfig,fullfile(saveItHere, sprintf('GI_withIndex_FA_Meansvs%s.svg',versus)),'svg');
close(bigfig)

%% HCP TEST-RETEST: MEAN FA BARS/TRT SCATTERPLOTS/GENERALIZATION PLOTS
% Plot bars with covar with Mean FA
fnameRoot = "FA_Means_withCoVdots_HCP_TRT";
createBars(dtTRT,'fnameRoot',fnameRoot,'saveItHere',string(saveItHere), ...
                    'saveSvg',true     ,'WahlOrder' ,true, 'HCPTRT',true)
% Plot bars with FA SD
fnameRoot = "FA_SD_HCP_TRT";
createBars(dtTRT, 'fnameRoot',fnameRoot, 'saveItHere',string(saveItHere), ...
           'meanOrSd',"SD", 'saveSvg',true, 'WahlOrder',true, 'HCPTRT',true)
% PLOT HCP TRT scatterplots and obtain ICC and rmse values
fnameRoot = "FA_Scatterplots_rmse_CoV_HCP_TRT";
[TRTtests,allCOVS] = createTRTscatterplots(dtTRT,unsMeansTRT,'fnameRoot',fnameRoot, ...
                                           'saveItHere',string(saveItHere), ...
                                           'meanOrSd',"SD", 'saveSvg',true,...
                                           'WahlOrder',true, 'HCPTRT',true)

%% MEAN FA DISTRIBUTIONS (normality, Cohens'd)

% Create distributions for b1000
fnameRoot = "FA_Distributions_b1000";
createDistributionsPlots(dtt1000,unsMeanst1000,'fnameRoot',fnameRoot,'saveItHere',string(saveItHere), ...
                         'saveSvg',true, 'WahlOrder',true, 'HCPTRT',false)

% Create distributions for b2000
fnameRoot = "FA_Distributions_b2000";
createDistributionsPlots(dtt2000,unsMeanst2000,'fnameRoot',fnameRoot,'saveItHere',string(saveItHere), ...
                         'saveSvg',true, 'WahlOrder',true, 'HCPTRT',false)

                     
                     
fnameRoot = "FA_Distributions_HCP_TRT";
createDistributionsPlots(dtTRT,unsMeansTRT,'fnameRoot',fnameRoot,'saveItHere',string(saveItHere), ...
                         'saveSvg',true, 'WahlOrder',true, 'HCPTRT',true)
       
%% WAHL 2010: correlation colors and CI tables: DATA PREP
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

%% Create the plots of interests
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

%% CALCULATE PERCENTAGE GIs
cd(saveItHere)
A = dir('GI_noIndex_noreference*.csv');A(:).name
GIs   = {};
Names = {};
for na=1:length(A)
    a = A(na);
    aName = a.name;
    % b = strrep(aName,'GI_noIndex_WHLorigvs','');
    % Names{na} = strrep(b,'_CI90_GRbars.csv','');
    b = strrep(aName,'GI_noIndex_noreferencevs','');
    Names{na} = strrep(b,'_CI90_bars.csv','');
    GIs{na} = dlmread(a.name)';
end
gitable = cell2table(GIs);
gitable.Properties.VariableNames = Names;

gitableS = stack(gitable, gitable.Properties.VariableNames, 'NewDataVariableName', 'GI');
gitableS.perc = NaN*ones(height(gitableS),1);
GIVALUES = .5;
for nu = 1:length(unique(gitableS.GI_Indicator))
    howMany = length(gitableS.GI{nu});
    higherThanGivalues = sum(gitableS.GI{nu} > GIVALUES);
    perc = 100*higherThanGivalues/howMany;
    gitableS.perc(nu) = perc;
end

%% TRT of the correlations
% TRT scatterplots
% HCP TEST-RETEST
tmp     = varfun(ostd,dtTRT,'GroupingVariables',{'Struct','SliceCats'},'InputVariables',{'meanVal'});
replRes = varfun(omean,dtTRT,'GroupingVariables',{'Struct','SliceCats'},'InputVariables',{'meanVal'});
replRes.Properties.VariableNames{'Fun_meanVal'} = 'Mean';
replRes.SD  = tmp.Fun_meanVal;
replRes     = replRes(:,[1,2,4,5]);
replRes     = unstack(replRes,{'Mean','SD'},{'SliceCats'});
replRes     = replRes([3,4,  11,12,   5,6,   7,8,   9,10,   1,2],:);
replResSortedTRT      = replRes(:,[1, 3,9,   2,8,   5,11,   4,10,  7,13,   6,12]);
dt = dtTRT;
cats      = categories(dt.SliceCats);
catcolors = unique(dt.SliceCatsRGB);
catcolors     = catcolors([1,1,2,2,3,3],:);

mrvNewGraphWin('Correlation TEST-RETEST Scatterplots'); 
% figHdl = figure('Name','FA profiles', ...
%                 'NumberTitle','off', ...
%                 'visible',   'on', ...
%                 'color','w', ...
%                 'Units','pixel', ...
%                 'Position',[0 0 1900 1100]);
long = longcell{CIrange};
for nt = 1: height(replResSortedTRT)
    scatter(long.HCP1000TEST, long.HCP1000RETEST,20,catcolors{1,:},'filled');hold on;
    scatter(long.HCP2000TEST, long.HCP2000RETEST,18,catcolors{3,:});
    scatter(long.HCP3000TEST, long.HCP3000RETEST,20,catcolors{5,:},'filled');
    xlim([-0.4, 1]);ylim([-0.4, 1]);
    xticks([-0.4:0.1:1]);yticks([-0.4:0.1:1]);
    % axis equal; % identityLine(gca);
    location='southeast';if nt >=11;location='southeast';end
    legend({sprintf('b=1000s/mm^2 (rmse=%.3f)',rmse1),...
            sprintf('b=2000s/mm^2 (rmse=%.3f)',rmse2),...
            sprintf('b=3000s/mm^2 (rmse=%.3f)',rmse3)}, ...
            'location',location)
    xlabel('TEST: Mean FA correlations');ylabel('RETEST: Mean FA correlations')
end    
set(gcf,'color','w');
suptitle({'Mean FA correlation scatterplots', 'HCP TEST-RETEST'})
saveas(gcf,fullfile(saveItHere, 'FACorrelation_Scatterplots_rmse_HCP_TRT.png'),'png');
close(gcf)   
    
% Create the table with the averages
[N1, Nold1, Xm1,Ym1,rho1,pval1,rhom1,pvalm1,rmse1,rmsem1,rrmse1,rrmsem1,sdX1,sdY1,sdXm1,sdYm1,icc1,iccm1] = dr_corrrmse(long.HCP1000TEST,long.HCP1000RETEST);
[N2, Nold2, Xm2,Ym2,rho2,pval2,rhom2,pvalm2,rmse2,rmsem2,rrmse2,rrmsem2,sdX2,sdY2,sdXm2,sdYm2,icc2,iccm2] = dr_corrrmse(long.HCP2000TEST,long.HCP2000RETEST);
[N3, Nold3, Xm3,Ym3,rho3,pval3,rhom3,pvalm3,rmse3,rmsem3,rrmse3,rrmsem3,sdX3,sdY3,sdXm3,sdYm3,icc3,iccm3] = dr_corrrmse(long.HCP3000TEST,long.HCP3000RETEST);

RMSE     = [rmse1,rmse2,rmse3]';
rRMSE    = [rrmse1,rrmse2,rrmse3]';
ICC      = [icc1,icc2,icc3]';
TRTtests = table(RMSE,rRMSE,ICC);
TRTtests.Properties.RowNames = {'b1000','b2000','b3000'};

%% To visualize all the data, plot some scatterplots as well
tractsOrder = { 'LeftCingulumCingulate'  , 'RightCingulumCingulate'  , ...
                'LeftArcuate'            , 'RightArcuate'            , ...
                'LeftIFOF'               , 'RightIFOF'               , ...
                'LeftILF'                , 'RightILF'                , ...
                'LeftUncinate'           , 'RightUncinate'           , ...
                'LeftCorticospinal'      , 'RightCorticospinal'      };
cats      = categories(unsMeans.SliceCats);
catcolors = unique(unsMeans.SliceCatsRGB);
for nc = 1:length(cats)
    corrMatTable = dr_createCorrelationMatrix(...
                               unsMeans(unsMeans.SliceCats==string(cats(nc)),:), ... 
                               'tractsOrder',tractsOrder, ...
                               'useCorrplot',true, 'plotIt', true, ...
                               'plotTitle', string(cats(nc)),  ...
                               'useBootstrap',false, 'nRep',10000, 'CIrange',95);
    set(gcf,'color','w');
    saveas(gcf,fullfile(saveItHere, sprintf('%s_FA_corrplots.png',string(cats(nc)))),'png');
    close(gcf)

end
cats      = categories(unsMeansrt.SliceCats);
catcolors = unique(unsMeansrt.SliceCatsRGB);
for nc = 1:length(cats)
    corrMatTable = dr_createCorrelationMatrix(...
                               unsMeansrt(unsMeansrt.SliceCats==string(cats(nc)),:), ... 
                               'tractsOrder',tractsOrder, ...
                               'useCorrplot',true, 'plotIt', true, ...
                               'plotTitle', strcat(string(cats(nc)),"RETEST"), ...
                               'useBootstrap',false, 'nRep',10000, 'CIrange',95);
    set(gcf,'color','w');
    saveas(gcf,fullfile(saveItHere, sprintf('%sRETEST_FA_corrplots.png',string(cats(nc)))),'png');
    close(gcf)

end

%% Generalizability (Consistency...) Index Line
plotIt = true; showLegend=false; showXnames=false;
mrvNewGraphWin('FA Correlation GI plots');
% figHdl = figure('Name','FA profiles', ...
%                 'NumberTitle','off', ...
%                 'visible',   'on', ...
%                 'color','w', ...
%                 'Units','pixel', ...
%                 'Position',[0 0 1900 1100]);
ncol=6;nrow=11;
for nc=1:length(unique(long.CorName))
    subplot(ncol,nrow,nc)
    corname  = string(long.CorName(nc));
    x = [1:100];
    y = [];
    for CIrange=1:100
        referenceColumn = 'NoReference';
        versus          = 'AllOptions';
        takeout =  {'WHLorig'};
        [WahlInsideCI, isOverlap, minUpper, maxLower, midCI, GI]=dr_compareCI(longcell{CIrange}, corname,...
                                                        'referenceColumn',referenceColumn,...
                                                        'plotIt',plotIt,'showLegend',showLegend, ...
                                                        'takeout',takeout, 'showXnames',showXnames);
        y = [y , GI];
    end
    plot(x,y); xlabel('Confidence Interval [0~100]'); ylabel('Consistency Index [0~1]');
    title(strrep(corname,'_','\_')); hold on;
    first = find(y>0);
    if ~isempty(first)
        first = first(1);
        if first==1;first=2;end;
        [maxGI, maxInd] = max(y);
        ylim([0,1.1]);
        l1 = line([x(first-1),x(first-1)],[0,1],'linewidth',2,'color','g','linestyle',':');
        l2 = line([x(1),x(100)],[maxGI,maxGI],'linewidth',2,'color','b','linestyle',':');
        l3 = line([x(maxInd),x(maxInd)],[maxGI-.05,maxGI+.05],'linewidth',2,'color','k','linestyle','-');
    end
end

%% Simulations with meanVar, distrNoise and OI
%{
% meanVar:    variance between the means of the distributions
% distrNoise: how pointy/flat is the probability distribution of the results of
%             the bootstrapping
% OI:         Overlap Index: the one called GI above. Relation between overlap
%             and the mean(CI) 

% Create the matrix of normal distributions
rng('default'); % Reset the random num generator every time to reproduce plots
bNum     = 10000;
plotHist = false;
A = zeros(bNum,100,100);
myMean = 0;
if plotHist; mrvNewGraphWin('Simulated Result Distributions');end
% figHdl = figure('Name','FA profiles', ...
%                 'NumberTitle','off', ...
%                 'visible',   'on', ...
%                 'color','w', ...
%                 'Units','pixel', ...
%                 'Position',[0 0 1900 1100]);
nh=0;
for a=1:size(A,2) 
    myMean = myMean + 0.01;
    mySD   = 0.1;
    for b=1:size(A,3) 
        mySD = mySD + 0.01;
        r = myMean + mySD .* randn(bNum,1);
        A(:,a,b) = r;
        if plotHist
            nh=nh+1;subplot(size(A,2),size(A,3),nh);hist(r);xlim([-.2,1.2]);ylim([0,4000])
        end
    end
end

simulateAgain = true;
% I created pretty big simulations, read them from here, otherwise just run it
% everytime
% readCacheFile = fullfile(paperReprPath,'local','cache','BmatrixWithSimulationsSD005.mat');
% readCacheFile = fullfile(paperReprPath,'local','cache','BmatrixWithSimulationsSD01.mat');
if simulateAgain

    % Now that we have all the possible normals, we need to combine them in the
    % plot. 
    % Define the required confidence intervals
    % B = zeros(4,size(A,2)^2,size(A,3)^2);
    B = zeros(4,size(A,2),size(A,3));
    % First 2 loops, navigate through all 10*10 values. Select each of them and then
    % loop again to do everything against everything
    
    ii=0;
    % Select one with the lowest mean and a mean SD to compare with the rest of
    % the values. 
    ae=1; % for ae=1:size(A,2);
        be=round(size(A,3)/2);  % for be=1:size(A,3) 
            ii=ii+1;
            jj=0;
            for ai=1:size(A,2)
                for bi=1:size(A,3) 
                    jj=jj+1;
                    % Obtain distributions
                    d1 = A(:,ae,be);
                    d2 = A(:,ai,bi);
                    % Do the calculation here
                    vOI=[];vmeanResult=[];vsdResult=[];vmeanCI=[];
                    for CIrange=1:100
                        twoTailedRange = (100 - CIrange)/2;
                        rawCI    = prctile([d1,d2], [twoTailedRange, 100-twoTailedRange]);
                        corrVals = [mean(d1), mean(d2)]; 
                        upperVals=rawCI(2,:); 
                        lowerVals=rawCI(1,:);
                        [OI, meanResult, sdResult, meanCI] = dr_GI(corrVals,upperVals,lowerVals);
                        vOI        =[vOI        ,OI        ];
                        vmeanResult=[vmeanResult,meanResult];
                        vsdResult  =[vsdResult  ,sdResult  ];
                        vmeanCI    =[vmeanCI    ,meanCI    ];
                    end
                    first = find(vOI>0);
                    if ~isempty(first)
                        first = first(1);
                        % if first==1;first=2;end;
                        % [maxGI, maxInd] = max(y);
                        % Write result
                        % B(:,ii,jj) = [first           ,vmeanResult(first), ...
                        %               vsdResult(first),vmeanCI(first)];
                        B(:,ai,bi) = [first           ,vmeanResult(first), ...
                                      vsdResult(first),vmeanCI(first)];
                    end

                end
            end
        % end
    % end
else
    load(readCacheFile)    
end

sdOfMean     = squeeze(B(3,:,:));
meanRes      = squeeze(B(2,:,:));
meanOfCIs    = squeeze(B(4,:,:));
CI1stOverlap = squeeze(B(1,:,:));

mrvNewGraphWin('Simulations with SD of mean, mean CIs and firstOverlapCI'); pointsize = 5;
% figHdl = figure('Name','FA profiles', ...
%                 'NumberTitle','off', ...
%                 'visible',   'on', ...
%                 'color','w', ...
%                 'Units','pixel', ...
%                 'Position',[0 0 1900 1100]);    
subplot(2,2,1);
    plot3(sdOfMean,meanOfCIs,CI1stOverlap,'b.'); hold on; grid
    xlabel('SD of Result Means'); ylabel('Mean CI at first overlap');zlabel('CI of first Overlap');

    subplot(2,2,2);
    scatter(sdOfMean(:),CI1stOverlap(:), pointsize, meanOfCIs(:));
    hold on;
    xlabel('SD of Result Means'); ylabel('CI of first Overlap')
    title('Color: Mean CI at min overlap'); colorbar

    subplot(2,2,3);
    scatter(meanOfCIs(:), sdOfMean(:), pointsize, CI1stOverlap(:));
    hold on;identityLine(gca); ylim([0,.75]);xlim([0,1.1])
    xlabel('Mean CI at min overlap'); ylabel('SD of Result Means'); 
    title('Color: CI of first overlap'); colorbar

    subplot(2,2,4);
    scatter(meanOfCIs(:),CI1stOverlap(:), pointsize, sdOfMean(:));
    xlabel('Mean CI at first overlap');ylabel('CI at first overlap'); 
    title('Color: SD of Means'); colorbar



   

% Remove some of the results, the distance between means between 0.5 and 0.6 is
% the same as the distance between 0.4 and 0.5 (and then all the variations
% inside)
sdOfMean     = squeeze(B(3,:,1:10));
meanOfCIs    = squeeze(B(4,:,1:10));
CI1stOverlap = squeeze(B(1,:,1:10));   
pointsize = 10;
scatter(sdOfMean(:),meanOfCIs(:), pointsize, CI1stOverlap(:));
hold on;% identityLine(gca); 
xlabel('SD of Result Means'); ylabel('Mean CI at min overlap');
title('Color: CI of first overlap'); colorbar
xlim([0.27,0.30]);ylim([0.39,0.43])








% Create some normal distributions for figures
rng('default'); % Reset the random num generator every time to reproduce plots
bNum     = 10000;
plotHist = false;
A = zeros(bNum,100,100);
myMean = 0;
nh=0;
for a=1:size(A,2) 
    myMean = myMean + 0.01;
    mySD   = 0.1;
    for b=1:size(A,3) 
        mySD = mySD + 0.01;
        r = myMean + mySD .* randn(bNum,1);
        A(:,a,b) = r;
    end
end
% GOOD Replicability - BAD Generalizability
X1 = A(:,1,10);X2 = A(:,45,10);
vals{1} = {X1,X2};
normalNames{1} = "GOODrep-BADgen";
% BAD Replicability - BAD Generalizability
X1 = A(:,1,40); X2 = A(:,45,40);
vals{2} = {X1,X2};
normalNames{2} = "BADrep-BADgen";
% GOOD Replicability - GOOD Generalizability
X1 = A(:,1,10);X2 = A(:,5,10);
vals{3} = {X1,X2};
normalNames{3} = "GOODrep-GOODgen";
% BAD Replicability - GOOD Generalizability
X1 = A(:,1,40);X2 = A(:,5,40);
vals{4} = {X1,X2};
normalNames{4} = "BADrep-GOODgen";

for nv=1:length(vals)
    X1 = vals{nv}{1,1};
    X2 = vals{nv}{1,2};
    [X1_values, MU1, SIGMA1, MN1, MX1] = dr_distPlottingVals(X1);
    [X2_values, MU2, SIGMA2, MN2, MX2] = dr_distPlottingVals(X2);
    PD1 = fitdist(X1,'normal');
    PD2 = fitdist(X2,'normal');
    group_pdf1 = pdf(PD1, X1_values);
    group_pdf2 = pdf(PD2, X2_values);
    darkerBlue  = [0, 0.1875, 1];
    lighterBlue = [0, 0.8125, 1];
    plot(X1_values,group_pdf1,'Color','k','LineStyle','-','LineWidth',4)
    hold on; axis off;set(gcf,'color','w');
    plot(MU1*[1,1],[0 max(group_pdf1)],'Color','k','LineStyle','-','LineWidth',2);
    xlim([-1.5,1.7]);ylim([0,2])
    plot(X2_values,group_pdf2,'Color','k','LineStyle','--','LineWidth',4);
    plot(MU2*[1,1],[0 max(group_pdf2)],'Color','k','LineStyle','--','LineWidth',2);
    % xkcdify(gca)
    set(gcf,'color','w');
    saveas(gcf,fullfile(saveItHere, sprintf('NormalComparisons_BW_%s.svg',normalNames{nv})),'svg');
    close(gcf)
end
%}

%% SIMULATION GRAPH
% close all; clear all; clc; 

E = 9;  %Number of experiments
mu0 = 2;  %ground truth of the parameter to be estimated
sigma_mu = 0.1;  %standard deviation of the means
sigma_e = 5;  %standard deviation of each experients
n = 44;  %number of samples used in each experiment; same for all experiment
B = 10000;  %Number of boostrap samples; in the simulation we use the theretical quantile instead because of the speed. The results are the same
perc = 0.9;  %The level of the CI
n_realization = 100;  %For each setting, GI is calculated for n_realiazation times and we record the mean;
ratio = 50;  % ratio of sigma_e/sigma_mu


% Generating heatmap: 
%Specify the candidate values for (sigma_mu,sigma_e,n,E) respectively
sigma_mu_cand = 0.025:0.0125:1;
sigma_e_cand  = 0.025:0.025:2;
n_cand        = [100];
E_cand        = [5];
%Initialization: GI_ne is a kxk cell, with each cell represent a setting of
%(n,e)
GI_ne = cell(size(n_cand,2), size(E_cand,2));
% The loop of different values of (n,e)
for nind = 1:size(n_cand,2)
    for Eind = 1:size(E_cand,2)
        n = n_cand(nind);
        E = E_cand(Eind); 
        
        
        % From now on, the individual plots
        GI_rec = zeros(size(sigma_mu_cand,2),size(sigma_e_cand,2));
        for i = 1:size(sigma_mu_cand,2)
            sigma_mu  = sigma_mu_cand(i);
            for j=1:size(sigma_e_cand,2)
                sigma_e = sigma_e_cand(j);
                % Calculate the index
                GI_temp = 0;
                %Repeat n_realization times and record the mean of GI
                for i_realization = 1:n_realization
                    mu = normrnd(mu0,sigma_mu,E,1);
                    mu_est = [];
                    CI = zeros(E,2);    
                    CI_theory = zeros(E,2);
                    result    = zeros(E,1);
                    %Generate E CIs
                    for e = 1:E
                        X = normrnd(mu(e),sigma_e,n,1);
                        % using boorstrap to compute the CI
                        % boots_rec = [];
                        % mu_est = [mu_est mean(X)];
                        %for b= 1:B
                        %    bootsind = randsample(1:n,n,true);
                        %    boots_rec = [boots_rec mean(X(bootsind)) - mean(X)];
                        %end
                        %CI(e,1) = mean(X) + quantile(boots_rec,(1-perc)/2);
                        %CI(e,2) = mean(X) + quantile(boots_rec,(1+perc)/2);

                        %check the bootstrap CI with theoretical CI
                        CI_theory(e,1) = mean(X) + norminv((1-perc)/2)*std(X)/sqrt(n);
                        CI_theory(e,2) = mean(X) + norminv((1+perc)/2)*std(X)/sqrt(n);
                        result         = mean(X);
                    end
                    % GI_temp = GI_temp + dr_GI(E,CI_theory,n);
                    upper = CI_theory(:,2);
                    lower = CI_theory(:,1);
                    GI_temp = GI_temp + dr_GI(result,upper,lower, ...
                                                        'LowerRange', -1,...
                                                        'UpperRange', 1);
                end
                % Record the mean for (i,j)
                GI_rec(i,j) = GI_temp/n_realization;
            end
        end
        
        % Put it on the plots
        GI_ne{nind, Eind} = GI_rec;
    end
end

% Plot the heatmap
mrvNewGraphWin('Simulations for the GI'); 
% figHdl = figure('Name','FA profiles', ...
%                 'NumberTitle','off', ...
%                 'visible',   'on', ...
%                 'color','w', ...
%                 'Units','pixel', ...
%                 'Position',[0 0 1900 1100]);
np = 0;
for nn = 1:size(n_cand,2)
    for nE = 1:size(E_cand,2)
        np = np + 1;
        subplot(size(n_cand,2),size(E_cand,2),np)
        imagesc(sigma_e_cand, sigma_mu_cand, GI_ne{nn,nE});
        colormap(cmapname)
        colorbar;
        set(gca,'YDir','normal');
        axis on;
        % title(sprintf('Heatmap for GI^{n:%i}_{E:%i}(CI:90%%)', n_cand(nn), E_cand(nE)))
        xlabel('Mean variance of the result distributions (\sigma_{E})')
        ylabel('Variance of the means (\sigma_{\mu})')
        xticks([.5:.5:4]);yticks([.5:.5:2])
        set(gca,'fontsize',16,'fontweight','bold')
    end
end
set(gcf,'color','w');
set(gcf,'WindowStyle','normal')
set(gcf, 'Units', 'inches', 'OuterPosition', [10, 10, 10, 10]);
% suptitle({'GI Simulations'})
saveas(gcf,fullfile(saveItHere, 'GI_Simulations_corr.png'),'png');
saveas(gcf,fullfile(saveItHere, 'GI_Simulations_corr.svg'),'svg');
close(gcf)


% Plot the normals, results of the upper corner
%{

%% Parameters
sigma_mu_cand = [0.1 2];
sigma_e_cand = [0.1 2];
mu0 = 2; 
E = 2;
n = 20;
perc = 0.9;
n_realization = 10000;
%% Generating data for plots
plot_all = cell(2,2);
for i = 1:2
    for j=1:2
        sigma_mu = sigma_mu_cand(i);
        sigma_e = sigma_e_cand(j);
        plot_mu = normrnd(mu0,sigma_mu,n_realization,1);
        plot_x = normrnd(plot_mu(1),sigma_e,n_realization,1);
        plot_all{i,j} = [plot_mu,plot_x];
        
    end
end



%%
subplot(2,2,3)
        h_mu = histfit(plot_all{1,1}(:,1));
        h_mu(1).FaceColor = [0.8 0.8 0];
        h_mu(1).EdgeColor = [1,1,1];
        h_mu(2).Color = [0,0,0];
        alpha(0)
        xlim([-8,12]);
        hold on;
        h_e = histfit(plot_all{1,1}(:,2));
        h_e(1).FaceColor = [0.8 0.8 0];
        h_e(1).EdgeColor = [1,1,1];
        h_e(2).LineStyle = '--';
        h_e(2).Color = [0,0,0];
        alpha(0)
        xlim([-8,12]);
        ylim([0,400]);
subplot(2,2,4)
        h_mu = histfit(plot_all{1,2}(:,1));
        h_mu(1).FaceColor = [0.8 0.8 0];
        h_mu(1).EdgeColor = [1,1,1];
        h_mu(2).Color = [0,0,0];
        alpha(0)
        xlim([-8,12]);
        ylim([0,400]);
        
        hold on;
        h_e = histfit(plot_all{1,2}(:,2));
        h_e(1).FaceColor = [0.8 0.8 0];
        h_e(1).EdgeColor = [1,1,1];
        h_e(2).LineStyle = '--';
        h_e(2).Color = [0,0,0];
        alpha(0)
        xlim([-8,12]);
        ylim([0,400]);

subplot(2,2,1)
        h_mu = histfit(plot_all{2,1}(:,1));
        h_mu(1).FaceColor = [0.8 0.8 0];
        h_mu(1).EdgeColor = [1,1,1];
        h_mu(2).Color = [0,0,0];
        alpha(0)
        xlim([-8,12]);
        ylim([0,400]);
        hold on;
        h_e = histfit(plot_all{2,1}(:,2));
        h_e(1).FaceColor = [0.8 0.8 0];
        h_e(1).EdgeColor = [1,1,1];
        h_e(2).LineStyle = '--';
        h_e(2).Color = [0,0,0];
        alpha(0)
        xlim([-8,12]);
        ylim([0,400]);
        
subplot(2,2,2)
        h_mu = histfit(plot_all{2,2}(:,1));
        h_mu(1).FaceColor = [0.8 0.8 0];
        h_mu(1).EdgeColor = [1,1,1];
        h_mu(2).Color = [0,0,0];
        alpha(0)
        xlim([-8,12]);
        ylim([0,400]);
        hold on;
        h_e = histfit(plot_all{2,2}(:,2));
        h_e(1).FaceColor = [0.8 0.8 0];
        h_e(1).EdgeColor = [1,1,1];
        h_e(2).LineStyle = '--';
        h_e(2).Color = [0,0,0];
        alpha(0)
        xlim([-8,12]);
        ylim([0,400]);
%}        


% Sensitivity analysis of c over the GI
% Generating plot data
c_cand = 1:10;
x = 0:0.001:2;
y = [];
group = [];
for c = c_cand
    y = [y;exp(-c*x)];
end

% Plotting
leg = cell(10,1);
for i =1:10
    plot(x,y(i,:));
    xlim([0,2]);
    leg{i} = ['c=' num2str(i)];
    hold on;
end
legend(leg);
ylabel('GI')
xlabel('$\displaystyle\frac{(L-O)}{(R_R-R_L)\sqrt{E}}$','Interpreter','latex')

% add shaded areas
% Add lines
h1 = line([0.05 0.05],[0 1]);
h2 = line([0.8 0.8],[0 1]);
% Set properties of lines
set([h1 h2],'Color','k','LineWidth',1,'LineStyle','-.')

%% Check the correlation example
LCST1000 =  unsMeans.LeftArcuate(unsMeans.SliceCats=='HCP1000');
RCST1000 =  unsMeans.RightArcuate(unsMeans.SliceCats=='HCP1000');
LCST2000 =  unsMeans.LeftArcuate(unsMeans.SliceCats=='HCP2000');
RCST2000 =  unsMeans.RightArcuate(unsMeans.SliceCats=='HCP2000');
LCST3000 =  unsMeans.LeftArcuate(unsMeans.SliceCats=='HCP3000');
RCST3000 =  unsMeans.RightArcuate(unsMeans.SliceCats=='HCP3000');
scatter(LCST1000,RCST1000,10,'b');hold on;
scatter(LCST2000,RCST2000,10,'r')
scatter(LCST3000,RCST3000,10,'g')
identityLine(gca)

std(LCST1000)
std(RCST1000)
std(LCST2000)
std(RCST2000)
std(LCST3000)
std(RCST3000)

%% Volume-FA (profile) comparisons for the 042 subject
dv = load(fullfile(paperReprPath,'local','cache','WestonHavens_Subj_042_volume.mat')); DV=dv.dt;
dt = load(fullfile(paperReprPath,'local','cache','WestonHavens_Subj_042_fa.mat'));     DT=dt.dt;

% for now, remove the b2000 result
DV = DV(DV.AcquMD.scanbValue~=2000,:);
DT = DT(DT.AcquMD.scanbValue~=2000,:);

mrvNewGraphWin('VOLUME profiles'); 
% figHdl = figure('Name','FA profiles', ...
%                 'NumberTitle','off', ...
%                 'visible',   'on', ...
%                 'color','w', ...
%                 'Units','pixel', ...
%                 'Position',[0 0 1900 1100]);
dt = DV;
tn = 'LeftThalamicRadiation';
tracts      = unique(dt.Struct);
repetitions = height(dt(dt.Struct==tn,:));
colors   = distinguishable_colors(repetitions,[0,0,0;1,1,1]);
% colors   = parula(repetitions);
for nt = 1: length(tracts)
    subplot(4,5,nt)
    tn = tracts(nt);
    values = dt{dt.Struct==tn,'Val'};
    for nc=1:repetitions
        plot([1:100],values(nc,:),'-','linewidth',2,'color',colors(nc,:));hold on;
    end
    ylabel('Volume'); % xlabel('Divisions')
    set(gca,'xtick',[])
    title(sprintf('%s',tn))
end
set(gcf,'color','w');
suptitle({'WH\_042\_volume\_profiles'})
saveas(gcf,fullfile(saveItHere, 'WH_042_volume_profiles.png'),'png');
close(gcf)
% FA
mrvNewGraphWin('FA profiles'); 
% figHdl = figure('Name','FA profiles', ...
%                 'NumberTitle','off', ...
%                 'visible',   'on', ...
%                 'color','w', ...
%                 'Units','pixel', ...
%                 'Position',[0 0 1900 1100]);
dt = DT;
tn = 'LeftThalamicRadiation';
tracts      = unique(dt.Struct);
repetitions = height(dt(dt.Struct==tn,:));
for nt = 1: length(tracts)
    subplot(4,5,nt);
    tn = tracts(nt);
    values = dt{dt.Struct==tn,'Val'};
    for nc=1:repetitions
        plot([1:100],values(nc,:),'-','linewidth',2,'color',colors(nc,:));hold on;
    end
    ylabel('FA'); % xlabel('Divisions')
    set(gca,'xtick',[])
    title(sprintf('%s',tn))
end
set(gcf,'color','w');
suptitle({'WH\_042\_fa\_profiles'})
saveas(gcf,fullfile(saveItHere, 'WH_042_fa_profiles.png'),'png');
close(gcf)
% SCATTERPLOTS
% FA
mrvNewGraphWin('VOL-FA scatterplots'); 
% figHdl = figure('Name','FA profiles', ...
%                 'NumberTitle','off', ...
%                 'visible',   'on', ...
%                 'color','w', ...
%                 'Units','pixel', ...
%                 'Position',[0 0 1900 1100]);
dt = DT;
tn = 'LeftThalamicRadiation';
tracts      = unique(dt.Struct);
repetitions = height(dt(dt.Struct==tn,:));
for nt = 1: length(tracts)
    subplot(4,5,nt)
    tn = tracts(nt);
    favalues = DT{DT.Struct==tn,'Val'}';falong = favalues(:);
    volvalues= DV{DV.Struct==tn,'Val'}';volong = volvalues(:);
    colortable = [];
    for nt=1:size(colors,1)
        tmp = repmat(colors(nt,:),[100,1]);
        colortable = [colortable;tmp];
    end
    [rho,pval] = corr(falong, volong);
    scatter(falong,volong,20,colortable,'filled')
    xlabel('FA'); ylabel('Volume')
    text(min(falong),min(volong),sprintf('Correlation=%.2f (p=$.3f)',rho,pval));
    title(sprintf('%s',tn))
end
set(gcf,'color','w');
suptitle({'WH\_042\_fa-vol\_scatterplot'})
saveas(gcf,fullfile(saveItHere, 'WH_042_fa-vol_scatterplot.png'),'png');
close(gcf)

%% Volume-FA (mean) comparisons for the 042 subject
dv = load(fullfile(paperReprPath,'local','cache','WestonHavens_Subj_042_volume.mat')); DV=dv.dt;
dt = load(fullfile(paperReprPath,'local','cache','WestonHavens_Subj_042_fa.mat'));     DT=dt.dt;

% for now, remove the b2000 result
DV = DV(DV.AcquMD.scanbValue~=2000,:);
DT = DT(DT.AcquMD.scanbValue~=2000,:);


% SCATTERPLOTS
% FA
mrvNewGraphWin('VOL-FA scatterplots'); 
% figHdl = figure('Name','FA profiles', ...
%                 'NumberTitle','off', ...
%                 'visible',   'on', ...
%                 'color','w', ...
%                 'Units','pixel', ...
%                 'Position',[0 0 1900 1100]);
dt = DT;
tn = 'LeftThalamicRadiation';
tracts      = unique(dt.Struct);
repetitions = height(dt(dt.Struct==tn,:));
for nt = 1: length(tracts)
    subplot(4,5,nt)
    tn = tracts(nt);
    falong = mean(DT{DT.Struct==tn,'Val'}')';
    volong = mean(DV{DV.Struct==tn,'Val'}')';
    [rho,pval] = corr(falong, volong);
    scatter(falong,volong,30,'filled');hold on;lsline
    xlabel('FA'); ylabel('Volume')
    text(min(falong),min(volong),sprintf('Correlation=%.2f (p=%.3f)',rho,pval));
    title(sprintf('%s',tn))
end
set(gcf,'color','w');
suptitle({'WH\_042\_fa-vol\_scatterplot'})
saveas(gcf,fullfile(saveItHere, 'WH_042_meanFa-vol_scatterplot.png'),'png');
close(gcf)


%% EXIT PROCEDURES
% Remove the repository from the Matlab Path
rmpath(genpath(paperReprPath));
if isdeployed
    % TODO: copy the figures to the output part. 
end
