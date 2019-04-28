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
saveItHere = string(fullfile(paperPath, '01_REVIEW01/Figures/REV01_v01/sources'));


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

% TODO: Right now the TEST and RETESTS will have the same colours, add it as an option
[dtALL, unsProfALL, unsMeansALL] = dr_createWorkingTables(DT, ...
                                 'sliceBasedOn',           {'Proj', 'SHELL','TRT'}, ...
                                 'SHELL',                  {'1000', '2000', '3000'}, ...
                                 'removeTractsContaining', {'Callosum'}, ... 
                                 'useThisTractsOnly',       useTracts, ...
                                 'createTractPairs',        true, ...
                                 'AGE',                     [14,58], ...
                                 'GENDER',                 {'male', 'female'});
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

%% Fig4A (FA PROFILES: HCP TRT) 

% REPLICATION EXPERIMENT: HCP TEST-RETEST (FIGURE 4A)
fnameRoot = "FA_Profiles_HCP_TEST_RETEST";
createProfile(dtTRT,'fnameRoot',fnameRoot,'saveItHere',string(saveItHere), ...
                    'saveSvg'  ,true     ,'WahlOrder' ,true, 'HCPTRT',true)

%% Fig5A1-5A2 (FA PROFILES: b1000 and b2000)
% GENERALIZATION EXPERIMENT: b1000 and b2000 (FIGURE 5A)
bvals     = [1000   , 2000];
dts       = {dtt1000, dtt2000};
% Calculate for b1000 and b2000
for nb=1:length(bvals)
    bval      = bvals(nb);
    fnameRoot = strcat("FA_Profiles_b", num2str(bval));
    createProfile(dts{nb},'saveSvg',true,'saveItHere',string(saveItHere),'fnameRoot',fnameRoot,'WahlOrder',true)
end

%% Fig xxx (MEAN FA BARS: b1000 and b2000)
% Plot bars with covar
fnameRoot = "FA_Means_withCoVdots";
createBars(dtt,'fnameRoot',fnameRoot,'saveItHere',string(saveItHere), ...
                    'saveSvg'  ,false     ,'WahlOrder' ,true, 'HCPTRT',false)

%% Fig 4D (FA GI: HCP) 
fnameRoot = "GI_Distribution_FA_HCP_TRT_ALL";
includeExp =  {'HCP1000TEST','HCP1000RETEST', ...
               'HCP2000TEST','HCP2000RETEST', ...
               'HCP3000TEST','HCP3000RETEST'};
GIplot(unsMeansALL, includeExp, 'useDistribution',true, 'ResultType','FA', ...
       'winSizeInch',[0,0,10,14], 'nrowcol',[3,4], 'ylab', 'FA', ...
       'fnameRoot',fnameRoot, 'saveItHere',string(saveItHere), 'saveSvg',false)
close(gcf)
   
%% Fig xxx (MEAN FA BARS: HCP TRT)
fnameRoot = "FA_Means_withCoVdots_HCP_TRT";
createBars(dtTRT,'fnameRoot',fnameRoot,'saveItHere',string(saveItHere), ...
                    'saveSvg',true     ,'WahlOrder' ,true, 'HCPTRT',true)
                
%% Fig xxx (SD FA BARS: HCP TRT)
fnameRoot = "FA_SD_HCP_TRT";
createBars(dtTRT, 'fnameRoot',fnameRoot, 'saveItHere',string(saveItHere), ...
           'meanOrSd',"SD", 'saveSvg',true, 'WahlOrder',true, 'HCPTRT',true)

%% Fig 4B (FA SCATTERPLOT: HCP TRT)
fnameRoot = "FA_Scatterplots_rmse_CoV_HCP_TRT_tigthAxes";
[TRTtests,allCOVS] = createTRTscatterplots(dtTRT,unsMeansTRT, ...
                                           'fnameRoot',fnameRoot, ...
                                           'saveItHere',string(saveItHere), 'saveSvg',true,...
                                           'WahlOrder',true, 'HCPTRT',true)

%% Fig 5B (FA Normal Distributions: b1000) 
fnameRoot = "FA_Distributions_b1000";
createDistributionsPlots(unsMeanst1000,'fnameRoot',fnameRoot, ...
                        'saveItHere',string(saveItHere), 'plotSum',false, ...
                        'nrowcol',[3,4], 'normalKsdensity','normal', ...
                        'ylims'  ,[0, 20],'xlims', [0.3, 0.66], ...
                        'saveSvg',true, 'WahlOrder',true, 'HCPTRT',false)

%% Fig 5B (FA Normal Distributions: b2000) 
fnameRoot = "FA_Distributions_b2000";
createDistributionsPlots(unsMeanst2000,'fnameRoot',fnameRoot,...
                        'saveItHere',string(saveItHere), 'plotSum',false, ...
                         'saveSvg',true, 'WahlOrder',true, 'HCPTRT',false)
                 
%% Fig 4C (FA Normal Distributions: HCP TRT)                
fnameRoot = "FA_Distributions_HCP_TRT";
createDistributionsPlots(unsMeansTRT,'fnameRoot',fnameRoot,'saveItHere',string(saveItHere), ...
                         'saveSvg',false, 'WahlOrder',true, 'HCPTRT',true)
       
%% Fig 5C1 (FA GI: b1000) 
fnameRoot = "GI_Distribution_FA_b1000";
includeExp =  {'WHL1000TEST','YWM1000TEST', 'HCP1000TEST'};
GIplot(unsMeansALL, includeExp, 'useDistribution',true, 'ResultType','FA', ...
       'winSizeInch',[0,0,10,14], 'nrowcol',[3,4], 'ylab', 'FA', ...
       'fnameRoot',fnameRoot, 'saveItHere',string(saveItHere), 'saveSvg',false)  
close(gcf)
  
%% Fig 5C2 (FA GI: b2000) 
fnameRoot = "GI_Distribution_FA_b2000";
includeExp =  {'YWM2000TEST','HCP2000TEST'};
GIplot(unsMeansALL, includeExp, 'useDistribution',true, 'ResultType','FA', ...
       'winSizeInch',[0,0,10,14], 'nrowcol',[3,4], 'ylab', 'FA', ...
       'fnameRoot',fnameRoot, 'saveItHere',string(saveItHere), 'saveSvg',true)
close(gcf)

%% Fig 6A (GI Correlation FA: HCP ALL) 
fnameRoot = "GI_Correlation_FA_HCP_ALL";
includeExp =  {'HCP1000TEST','HCP1000RETEST', ...
               'HCP2000TEST','HCP2000RETEST', ...
               'HCP3000TEST','HCP3000RETEST'};
GIplot(unsMeansALL, includeExp, 'useDistribution',false, 'ResultType','Correlation', ...
       'winSizeInch',[0,0,16,4], 'nrowcol',[1,6], 'ylab', 'L-R Correlation', ...
       'calcType','correlation', 'nRep',5000, 'onlyBilateral',true,...
       'fnameRoot',fnameRoot, 'saveItHere',string(saveItHere), 'saveSvg',true, 'savePng',true)
close(gcf)
                     
%% Fig 6B (GI Correlation FA: HCP b1000) 
fnameRoot = "GI_Correlation_FA_HCP_b1000";
includeExp =  {'HCP1000TEST','HCP1000RETEST'};
GIplot(unsMeansALL, includeExp, 'useDistribution',false, 'ResultType','Correlation', ...
       'winSizeInch',[0,0,16,4], 'nrowcol',[1,6], 'ylab', 'L-R Correlation', ...
       'calcType','correlation', 'nRep',5000, 'onlyBilateral',true,...
       'fnameRoot',fnameRoot, 'saveItHere',string(saveItHere), 'saveSvg',true, 'savePng',true)
close(gcf)

%% Fig 6C (GI Correlation FA: HCP b2000) 
fnameRoot = "GI_Correlation_FA_HCP_b2000";
includeExp =  {'HCP2000TEST','HCP2000RETEST'};
GIplot(unsMeansALL, includeExp, 'useDistribution',false, 'ResultType','Correlation', ...
       'winSizeInch',[0,0,16,4], 'nrowcol',[1,6], 'ylab', 'L-R Correlation', ...
       'calcType','correlation', 'nRep',5000, 'onlyBilateral',true,...
       'fnameRoot',fnameRoot, 'saveItHere',string(saveItHere), 'saveSvg',true, 'savePng',true)
close(gcf)

%% Fig 6D (GI Correlation FA: HCP b3000) 
fnameRoot = "GI_Correlation_FA_HCP_b3000";
includeExp =  {'HCP3000TEST','HCP3000RETEST'};
GIplot(unsMeansALL, includeExp, 'useDistribution',false, 'ResultType','Correlation', ...
       'winSizeInch',[0,0,16,4], 'nrowcol',[1,6], 'ylab', 'L-R Correlation', ...
       'calcType','correlation', 'nRep',5000, 'onlyBilateral',true,...
       'fnameRoot',fnameRoot, 'saveItHere',string(saveItHere), 'saveSvg',true, 'savePng',true)
close(gcf)

%% Fig 8A (GI Correlation FA: ALL b1000) 
fnameRoot = "GI_Correlation_FA_ALL_b1000";
includeExp =  {'WHL1000TEST','YWM1000TEST','HCP1000TEST','HCP1000RETEST'};
GIplot(unsMeansALL, includeExp, 'useDistribution',false, 'ResultType','Correlation', ...
       'winSizeInch',[0,0,16,4], 'nrowcol',[1,6], 'ylab', 'L-R Correlation', ...
       'calcType','correlation', 'nRep',5000, 'onlyBilateral',true,...
       'fnameRoot',fnameRoot, 'saveItHere',string(saveItHere), 'saveSvg',true, 'savePng',true)
close(gcf)

%% Fig 8B (GI Correlation FA: ALL b2000) 
fnameRoot = "GI_Correlation_FA_ALL_b2000";
includeExp =  {'YWM2000TEST','HCP2000TEST','HCP2000RETEST'};
GIplot(unsMeansALL, includeExp, 'useDistribution',false, 'ResultType','Correlation', ...
       'winSizeInch',[0,0,16,4], 'nrowcol',[1,6], 'ylab', 'L-R Correlation', ...
       'calcType','correlation', 'nRep',5000, 'onlyBilateral',true,...
       'fnameRoot',fnameRoot, 'saveItHere',string(saveItHere), 'saveSvg',true, 'savePng',true)
close(gcf)

%% Fig 9 (FA LR SCATTERPLOT: HCP TRT)

fnameRoot = "FA_Scatterplots_LR_HCP_TRT";
includeExp =  {'WHL1000TEST','YWM1000TEST', 'YWM2000TEST', ...
               'HCP1000TEST','HCP1000RETEST', ...
               'HCP2000TEST','HCP2000RETEST', ...
               'HCP3000TEST','HCP3000RETEST'};
createLRscatterplots(unsMeansALL, includeExp, ...
                    'winSizeInch',[0,0,16,5], 'nrowcol',[1,6], ...
                    'ylab', 'Right FA', 'xlab', 'Left FA', 'WahlOrder',true, ...
                    'xlimits', [0.3,0.7], 'ylimits', [0.3,0.7],...
                    'fnameRoot',fnameRoot, 'saveItHere',string(saveItHere), 'saveSvg',false,'savePng',false)
 
%% Fig 9 (GI_RawResiduals_FA_ALL_NoBootstrap_Individually) 
fnameRoot = "GI_RawResiduals_FA_ALL_NoBootstrap_Individually";
includeExp =  {'WHL1000TEST','YWM1000TEST','YWM2000TEST', ...
               'HCP1000TEST','HCP1000RETEST', ...
               'HCP2000TEST','HCP2000RETEST', ...
               'HCP3000TEST','HCP3000RETEST'};
GIplot(unsMeansALL, includeExp, 'useDistribution',false, 'ResultType','FARegressionResiduals', ...
       'winSizePix',[0,0,1900,500], 'nrowcol',[1,6], 'ylab', 'L-R FA Resid', ...
       'showXnames',true, 'normResidual', false, 'Group2Individual',false,...
       'calcType','RegressionResiduals', 'nRep',5000, 'onlyBilateral',true,...
       'fnameRoot',fnameRoot,'saveItHere',saveItHere,'saveSvg',false,'savePng',true)
close(gcf)

%% Fig 9 (GI_NormalizedResiduals_FA_ALL_NoBootstrap_Individually) 
fnameRoot = "GI_NormalizedResiduals_FA_ALL_NoBootstrap_Individually";
includeExp =  {'WHL1000TEST','YWM1000TEST','YWM2000TEST', ...
               'HCP1000TEST','HCP1000RETEST', ...
               'HCP2000TEST','HCP2000RETEST', ...
               'HCP3000TEST','HCP3000RETEST'};
GIplot(unsMeansALL, includeExp, 'useDistribution',false, 'ResultType','FARegressionNormResiduals', ...
       'winSizePix',[0,0,1900,500], 'nrowcol',[1,6], 'ylab', 'L-R FA Resid', ...
       'showXnames',true, 'normResidual', true, ...
       'calcType','RegressionResiduals', 'nRep',5000, 'onlyBilateral',true,...
       'fnameRoot',fnameRoot,'saveItHere',saveItHere,'saveSvg',false,'savePng',false)
% close(gcf)

%% Fig 9 ("GI_RawResiduals_FA_OnlyHCP_NoBootstrap_IndividualExpAndALL") 
fnameRoot = "GI_RawResiduals_FA_ALL_NoBootstrap_IndividualExpAndALL_line";
includeExp =  {'WHL1000TEST','YWM1000TEST','YWM2000TEST', ...
                'HCP1000TEST', ...
               'HCP2000TEST',...
               'HCP3000TEST'};
groupStats = GIplot(unsMeansALL, includeExp, 'useDistribution',false, 'ResultType','FARegressionResiduals', ...
       'winSizePix',[0,0,1900,500], 'nrowcol',[1,6], 'ylab', 'L-R FA Resid', ...
       'showXnames',true, 'normResidual', false, 'Group2Individual',false,...
       'calcType','RegressionResiduals', 'nRep',5000, 'onlyBilateral',true,...
       'plotGroupBar', true, ...
       'fnameRoot',fnameRoot,'saveItHere',saveItHere,'saveSvg',false,'savePng',true)
fnameRoot = "IS_FA_ALL_line";
findSpace(groupStats.CI90, ...
       'winSizePix',[0,0,1900,500], 'nrowcol',[1,6], 'ylab', 'Right FA ', 'showLineForm',true, ...
       'xlab', 'Left FA ', 'LineColor',[0 0 0],'LineStyle','-', ...
       'useStdAllBS','BS', 'howManySD',1 ,'CIperc',90, ...
       'fnameRoot',fnameRoot,'saveItHere',saveItHere,'saveSvg',false,'savePng',true)

%% Fig 9 ("GI_RawResiduals_FA_ALL_NoBootstrap_IndividualExpAndALL") 
fnameRoot = "GI_RawResiduals_FA_ALL_NoBootstrap_IndividualExpAndALL_line";
includeExp =  {'WHL1000TEST','YWM1000TEST','YWM2000TEST', ...
               'HCP1000TEST','HCP1000RETEST', ...
               'HCP2000TEST','HCP2000RETEST', ...
               'HCP3000TEST','HCP3000RETEST'};
groupStats = GIplot(unsMeansALL, includeExp, 'useDistribution',false, 'ResultType','FARegressionResiduals', ...
       'winSizePix',[0,0,1900,500], 'nrowcol',[1,6], 'ylab', 'L-R FA Resid', ...
       'showXnames',true, 'normResidual', false, 'Group2Individual',false,...
       'calcType','RegressionResiduals', 'nRep',5000, 'onlyBilateral',true,...
       'plotGroupBar', true, ...
       'fnameRoot',fnameRoot,'saveItHere',saveItHere,'saveSvg',false,'savePng',true)
fnameRoot = "IS_FA_ALL_line";
findSpace(groupStats.CI90, ...
       'winSizePix',[0,0,1900,500], 'nrowcol',[1,6], 'ylab', 'Right FA ', 'showLineForm',true, ...
       'xlab', 'Left FA ', 'LineColor',[0 0 0],'LineStyle','-', ...
       'useStdAllBS','BS', 'howManySD',1 ,'CIperc',90, ...
       'fnameRoot',fnameRoot,'saveItHere',saveItHere,'saveSvg',false,'savePng',true)


   
fnameRoot = "GI_RawResiduals_FA_b1000_NoBootstrap_IndividualExpAndALL_line";
includeExp =  {'WHL1000TEST','YWM1000TEST', ...
               'HCP1000TEST','HCP1000RETEST'};
groupStats = GIplot(unsMeansALL, includeExp, 'useDistribution',false, 'ResultType','FARegressionResiduals', ...
       'winSizePix',[0,0,1900,500], 'nrowcol',[1,6], 'ylab', 'L-R FA Resid', ...
       'showXnames',true, 'normResidual', false, 'Group2Individual',false,...
       'calcType','RegressionResiduals', 'nRep',5000, 'onlyBilateral',true,...
       'plotGroupBar', true, ...
       'fnameRoot',fnameRoot,'saveItHere',saveItHere,'saveSvg',false,'savePng',true)
fnameRoot = "IS_FA_b1000_line";
findSpace(groupStats.CI90, ...
       'winSizePix',[0,0,1900,500], 'nrowcol',[1,6], 'ylab', 'Right FA ', 'showLineForm',true, ...
       'xlab', 'Left FA ', 'LineColor',[.5 .5 .5 ],'LineStyle','-.', ...
       'useStdAllBS','BS', 'howManySD',1 ,'CIperc',90, ...
       'fnameRoot',fnameRoot,'saveItHere',saveItHere,'saveSvg',false,'savePng',true)
   
   
   
fnameRoot  = "GI_RawResiduals_FA_b2000_NoBootstrap_IndividualExpAndALL_line";
includeExp =  {'YWM2000TEST', 'HCP2000TEST'};
groupStats = GIplot(unsMeansALL, includeExp, 'useDistribution',false, 'ResultType','FARegressionResiduals', ...
       'winSizePix',[0,0,1900,500], 'nrowcol',[1,6], 'ylab', 'L-R FA Resid', ...
       'showXnames',true, 'normResidual', false, 'Group2Individual',true,...
       'calcType','RegressionResiduals', 'nRep',5000, 'onlyBilateral',true,...
       'plotGroupBar', true,  ...
       'fnameRoot',fnameRoot,'saveItHere',saveItHere,'saveSvg',false,'savePng',true)
fnameRoot = "IS_FA_b2000_line";
findSpace(groupStats.CI90, ...
       'winSizePix',[0,0,1900,500], 'nrowcol',[1,6], 'ylab', 'Right FA ', 'showLineForm',true, ...
       'xlab', 'Left FA ', 'LineColor',[.1 .1 .1 ],'LineStyle','-', ...
       'useStdAllBS','SD', 'howManySD',[2,1] ,'CIperc',[95,68], ...
       'plotSubj',true, 'SubjValues', SubjValues,...
       'fnameRoot',fnameRoot,'saveItHere',saveItHere,'saveSvg',false,'savePng',true) 
% close all

%% Fig 10 (FA Normal Distributions of Residuals HCP TRT)   
fnameRoot  = "GI_RawResiduals_HCP_TRT_Distributions";
includeExp =  {'HCP1000TEST', 'HCP2000TEST', 'HCP3000TEST', ...
               'HCP1000RETEST', 'HCP2000RETEST', 'HCP3000RETEST'};
groupStats = GIplot(unsMeansALL, includeExp, 'useDistribution',false, ...
       'ResultType','FARegressionResiduals', 'HCPTRT', false, ...
       'winSizePix',[0,0,1900,500], 'nrowcol',[1,6], 'ylab', 'Count', ...
       'showXnames',false, 'normResidual', false, 'Group2Individual',true,...
       'calcType','RegressionResiduals', 'nRep',5000, 'onlyBilateral',true,...
       'plotGroupBar', true, 'GIcoloringMethod', 'distributions', ...
       'fnameRoot',fnameRoot,'saveItHere',saveItHere,'saveSvg',false,'savePng',false)
   
%% Fig 10 (FA Normal Distributions of Residuals)  
fnameRoot  = "GI_RawResiduals_ALL_Distributions";
includeExp =  {'WHL1000TEST','YWM1000TEST' , 'YWM2000TEST', ...
               'HCP1000TEST', 'HCP2000TEST', 'HCP3000TEST'};
tractsOrder = { 'LeftArcuate'            , 'RightArcuate'            , ...
                'LeftIFOF'               , 'RightIFOF' ,...
                'LeftCorticospinal'      , 'RightCorticospinal'};
groupStats = GIplot(unsMeansALL, includeExp, 'tractsOrder',tractsOrder, ...
       'useDistribution',false, 'ResultType','FARegressionResiduals', ...
       'winSizePix',[0,0,1900,500], 'nrowcol',[1,6], 'ylab', 'Count', ...
       'showXnames',false, 'normResidual', false, 'Group2Individual',true,...
       'calcType','RegressionResiduals', 'nRep',5000, 'onlyBilateral',true,...
       'plotGroupBar', true, 'GIcoloringMethod', 'distributions', 'plotIt', true, ...
       'fnameRoot',fnameRoot,'saveItHere',saveItHere,'saveSvg',false,'savePng',false)
fnameRoot = "IS_FA_b2000_line";
findSpace(groupStats.CI90, ...
       'winSizePix',[0,0,1900,500], 'nrowcol',[1,6], 'ylab', 'Right FA ', 'showLineForm',false, ...
       'xlab', 'Left FA ', 'LineColor',[.1 .1 .1 ],'LineStyle','-', ...
       'ylims', [.3,.7], 'xlims',[.3,.7], ...
       'useStdAllBS','SD', 'howManySD',[2,1] ,'CIperc',[95,68], ...
       'plotSubj',true, 'SubjValues', unsMeansALL, 'findBand',true ,...
       'fnameRoot',fnameRoot,'saveItHere',saveItHere,'saveSvg',false,'savePng',false)   
   
%% Plot all the GI plots to have the GI for all the exp that we are interested

includeExps = {...
               {'HCP1000TEST','HCP1000RETEST'},...
               {'HCP2000TEST','HCP2000RETEST'},...
               {'HCP3000TEST','HCP3000RETEST'},...
               {'HCP1000TEST','HCP1000RETEST','HCP2000TEST','HCP2000RETEST','HCP3000TEST','HCP3000RETEST'},...
               {'WHL1000TEST','YWM1000TEST'  ,'HCP1000TEST'}, ...
               {'YWM2000TEST','HCP2000TEST'} , ...
               {'WHL1000TEST','YWM1000TEST'  , 'YWM2000TEST', 'HCP1000TEST', 'HCP2000TEST', 'HCP3000TEST'}...
              };
tractsOrder = { 'LeftCingulumCingulate'  , 'RightCingulumCingulate'  , ...
                'LeftArcuate'            , 'RightArcuate'            , ...
                'LeftIFOF'               , 'RightIFOF'               , ...
                'LeftILF'                , 'RightILF'                , ...
                'LeftUncinate'           , 'RightUncinate'           , ...
                'LeftCorticospinal'      , 'RightCorticospinal'      };
% tractsOrder = { 'LeftArcuate'            , 'RightArcuate'            , ...
%                 'LeftIFOF'               , 'RightIFOF'               };

% Create table to store all the ranges
% toStr      = @(x) string(strrep(strjoin(includeExp),' ','_'));
% allExpStrs = cellfun(toStr, includeExps);
% VariableNames = ["Struct", allExpStrs];
% RangeTable = array2table(nan(length(tractsOrder), length(VariableNames)), ...
%                         'VariableNames', VariableNames);

for nie = 1:length(includeExps)
    includeExp =  includeExps{nie};
    fnameRoot = string(strrep(strjoin(includeExp),' ','_'));
    [groupStats, Range] = GIplot(unsMeansALL, includeExp, 'tractsOrder',tractsOrder, ...
        'useDistribution',true, 'ResultType','FA', 'CIrange',95,'GIcoloringMethod', 'barsdots',...
        'winSizeInch',[0,0,10,10], 'nrowcol',[2,3], 'ylab', 'FA', ...
        'plotIt', true, 'showLineForm', false, 'plotIndex' , true, ...
        'fnameRoot',fnameRoot, 'saveItHere',string(saveItHere), 'saveSvg',false)  
end   
   
   
% Plot the GI for the Right FA residuals after fitting with the left

for nie = 1:length(includeExps)
    % includeExp =  includeExps{nie};
    includeExp = {'WHL1000TEST','YWM1000TEST'  , 'YWM2000TEST', 'HCP1000TEST', 'HCP2000TEST', 'HCP3000TEST'}
    % includeExp = {'HCP1000TEST','HCP1000RETEST','HCP2000TEST','HCP2000RETEST','HCP3000TEST','HCP3000RETEST'};
    fnameRoot = string(strrep(strjoin(includeExp),' ','_'));
    groupStats = GIplot(unsMeansALL, includeExp, 'tractsOrder',tractsOrder, ...
           'useDistribution',false, 'ResultType','FARegressionResiduals', ...
           'winSizePix',[0,0,1900,500], 'nrowcol',[1,6], 'ylab', 'Right FA', ...
           'showXnames',false, 'normResidual', false, 'Group2Individual',true,...
           'calcType','RegressionResiduals', 'nRep',5000, 'onlyBilateral',true,...
           'plotGroupBar', true, 'GIcoloringMethod', 'barsdots','CIrange',95, 'plotIt', false,...
           'fnameRoot',fnameRoot,'saveItHere',saveItHere,'saveSvg',false,'savePng',false) 
    groupStats = groupStats.CI95; 
    allSubjValues = unsMeansALL(ismember(unsMeansALL.SliceCats,includeExp),:);
    % Plot for testing something
    fnameRoot = strcat("ContourPlot_2x3_",fnameRoot)
    % 'winSizePix',[0,0,1900,500], 'nrowcol',[1,6],
    groupStats = findSpace(groupStats, ...
       'winSizePix',[0,0,1000,1000], 'nrowcol',[2,3], 'ylab', 'Right FA ', ... 
       'showLineForm',false, ...
       'xlab', 'Left FA ', 'LineColor',[.1 .1 .1 ],'LineStyle','-', ...
       'ylims', [.25,.75], 'xlims',[.25,.75], 'plotIt',true, 'addEllipse',true,...
       'useStdAllBS','BS', 'howManySD',[2,1] ,'CIperc',[95,68], ...
       'plotSubj',true, 'SubjValues', allSubjValues, 'findBand',true ,...
       'fnameRoot',fnameRoot,'saveItHere',saveItHere,'saveSvg',true,'savePng',false)  
   % {
   % Obtain rmse for test-retest
   trt = sortrows(allSubjValues,{'SliceCats','SubjID'});
   % Obtain LeftRight, TestRetest IFOF
   for bval={"1000","2000","3000"}
       for tn={"IFOF"}
           LT   = trt.(strcat("Left",tn{:}))(trt.SliceCats==strcat("HCP",bval{:},"TEST"));
           RT   = trt.(strcat("Right",tn{:}))(trt.SliceCats==strcat("HCP",bval{:},"TEST"));
           LRT  = trt.(strcat("Left",tn{:}))(trt.SliceCats==strcat("HCP",bval{:},"RETEST"));
           RRT  = trt.(strcat("Right",tn{:}))(trt.SliceCats==strcat("HCP",bval{:},"RETEST"));
           sprintf('Left %s%4.0f, RMSE: %.3f',tn{:},bval{:},mean(sqrt((LT-LRT).^2)))
           sprintf('%s%4.0f, RMSE: %.3f',tn{:},bval{:},mean(sqrt((LT-LRT).^2 + (RT-RRT).^2)))
           
           sprintf('Left %s%4.0f, SD: %.3f', tn{:}, bval{:}, std(abs(LT-LRT)))
           sprintf('Left TEST %s%4.0f, SD: %.3f', tn{:}, bval{:}, std(LT))
           sprintf('Left RETEST %s%4.0f, SD: %.3f', tn{:}, bval{:}, std(LRT))
           sprintf('%s%4.0f, SD: %.3f',tn{:}, bval{:}, std(sqrt((LT-LRT).^2 + (RT-RRT).^2))) 
       end
   end
       
    
   
   % Correlations
       for tn=1:length(groupStats.CorName)
           tract = string(groupStats.CorName(tn));
           Left  = allSubjValues.(strcat("Left",tract));
           Right = allSubjValues.(strcat("Right",tract));
           sprintf('%s corr: %.2f',tract, corr(Left,Right,'rows','pairwise')) 
       end

   
   
   
   
   
   %}
end

%% Plot the test retest in the residuals
% Fit all with the same line and then plot one residuals agains the others.
includeExp  = {'HCP1000TEST','HCP1000RETEST','HCP2000TEST','HCP2000RETEST','HCP3000TEST','HCP3000RETEST'};
tractsOrder = { 'LeftArcuate'            , 'RightArcuate'            , ...
                'LeftIFOF'               , 'RightIFOF'               };
fnameRoot = "IsoResiduals_withAll_TEST-RETEST_HCP_datapoints";
groupStatsTESTRETEST = GIplot(unsMeansALL, includeExp, 'tractsOrder',tractsOrder, ...
       'useDistribution',false, 'ResultType','FARegressionResiduals', ...
       'winSizePix',[0,0,1900,500], 'nrowcol',[1,6], 'ylab', 'Right FA', ...
       'showXnames',false, 'normResidual', false, 'Group2Individual',true,...
       'calcType','RegressionResiduals', 'nRep',5000, 'onlyBilateral',true,...
       'plotGroupBar', true, 'GIcoloringMethod', 'barsdots','CIrange',95, 'plotIt', false,...
       'fnameRoot',fnameRoot,'saveItHere',saveItHere,'saveSvg',false,'savePng',false) 
    groupStatsTESTRETEST = groupStatsTESTRETEST.CI95; 
    allSubjValues = unsMeansALL(ismember(unsMeansALL.SliceCats,includeExp),:);
    % Plot for testing something
    fnameRoot = strcat("ContourPlot_2x3_",fnameRoot)
    % 'winSizePix',[0,0,1900,500], 'nrowcol',[1,6],
    groupStats = findSpace(groupStatsTESTRETEST, ...
       'winSizePix',[0,0,1000,1000], 'nrowcol',[2,3], 'ylab', 'Right FA ', ... 
       'showLineForm',false, ...
       'xlab', 'Left FA ', 'LineColor',[.1 .1 .1 ],'LineStyle','-', ...
       'ylims', [.25,.75], 'xlims',[.25,.75], 'plotIt',true, ...
       'useStdAllBS','BS', 'howManySD',[2,1] ,'CIperc',[95,68], ...
       'plotSubj',true, 'SubjValues', allSubjValues, 'findBand',true ,...
       'fnameRoot',fnameRoot,'saveItHere',saveItHere,'saveSvg',true,'savePng',false)  


% plot only the IFOF but be sure that subject names are the same
residuals = table(groupStatsTEST{groupStatsTEST.CorName=="IFOF",'residuals'}{:}','VariableNames',{'Residuals'});
residuals.SliceCats = groupStatsTEST{groupStatsTEST.CorName=="IFOF",'SliceCats'}{:}';
residuals.SubjID    = groupStatsTEST{groupStatsTEST.CorName=="IFOF",'SubjID'}{:}';
residuals = sortrows(residuals,{'SubjID','SliceCats'});
plotIt=true;
fnameRoot = "Residuals_TEST-RETEST_HCP";
if plotIt
    b1000T    = residuals(residuals.SliceCats=='HCP1000TEST',{'SubjID'});
    b1000RT   = residuals(residuals.SliceCats=='HCP1000RETEST',{'SubjID'});
    b1000T    = sortrows(b1000T,'SubjID'); b1000RT=sortrows(b1000RT,'SubjID');
    if isequal(b1000T.SubjID,b1000RT.SubjID)
        b1000T = residuals{residuals.SliceCats=='HCP1000TEST',{'Residuals'}};
        b1000RT = residuals{residuals.SliceCats=='HCP1000RETEST',{'Residuals'}};
    else
        error('The test-retest is not comparing the same subjects')
    end

    b2000T    = residuals(residuals.SliceCats=='HCP2000TEST',{'SubjID'});
    b2000RT   = residuals(residuals.SliceCats=='HCP2000RETEST',{'SubjID'});
    b2000T    = sortrows(b2000T,'SubjID'); b2000RT=sortrows(b2000RT,'SubjID');
    if isequal(b2000T.SubjID,b2000RT.SubjID)
        b2000T = residuals{residuals.SliceCats=='HCP2000TEST',{'Residuals'}};
        b2000RT = residuals{residuals.SliceCats=='HCP2000RETEST',{'Residuals'}};
    else
        error('The test-retest is not comparing the same subjects')
    end    
    
    b3000T    = residuals(residuals.SliceCats=='HCP3000TEST',{'SubjID'});
    b3000RT   = residuals(residuals.SliceCats=='HCP3000RETEST',{'SubjID'});
    b3000T    = sortrows(b3000T,'SubjID'); b3000RT=sortrows(b3000RT,'SubjID');
    if isequal(b3000T.SubjID,b3000RT.SubjID)
        b3000T = residuals{residuals.SliceCats=='HCP3000TEST',{'Residuals'}};
        b3000RT = residuals{residuals.SliceCats=='HCP3000RETEST',{'Residuals'}};
    else
        error('The test-retest is not comparing the same subjects')
    end   
    
    bigfig = figure('Name',fnameRoot, ...
                            'NumberTitle','off', ...
                            'visible',   'on', ...
                            'color','w', ...
                            'WindowStyle','normal', ...
                            'Units','pixel', ...
                            'Position',[0,0,400,400])
    scatter(b1000T,b1000RT,50,[0, 0, 128]/255    ,'filled'); hold on;
    scatter(b2000T,b2000RT,50,[0, 130, 200]/255  ,'filled');
    scatter(b3000T,b3000RT,50,[70, 240, 240]/255 ,'filled');
    % Fit ellipses
    fitEllipse(b1000T,b1000RT,[0, 0, 128]/255)
    fitEllipse(b2000T,b2000RT,[0, 130, 200]/255)
    fitEllipse(b3000T,b3000RT,[70, 240, 240]/255)
    % xlim([0.3, 0.66]);ylim([0.3, 0.66]);
    % xticks([0.3:0.1:0.6]);yticks([0.3:0.1:0.6]);
    identityLine(gca);
    axis equal
    xlim([-0.065, 0.065]);ylim([-0.065, 0.065]);
    xticks(-0.06:0.02:0.06);
    yticks(-0.06:0.02:0.06);
    xlabel('\Delta FA (TEST)'  ,'FontWeight','bold');
    ylabel('\Delta FA (RETEST)','FontWeight','bold');
    set(gca,'FontSize',18)
    grid off
    saveas(gcf,fullfile(saveItHere, strcat(fnameRoot,'.svg')),'svg');
    
end
                   
%% Plot a dot in the prediction plots
% Prepare one subject
SubjValues = unsMeansALL(unsMeansALL.SubjID=='NIH-TBI0F', contains(unsMeansALL.Properties.VariableNames, {'Left','Right'}));
SubjValues = [SubjValues;SubjValues];
SubjValues{2,contains(SubjValues.Properties.VariableNames,{'Left'})} = ...
     0.075 + SubjValues{2,contains(SubjValues.Properties.VariableNames,{'Left'})};
SubjValues{2,contains(SubjValues.Properties.VariableNames,{'Right'})} = ...
     -0.075 + SubjValues{2,contains(SubjValues.Properties.VariableNames,{'Right'})};
 
fnameRoot  = "GI_RawResiduals_FA_allTEST_Bootstrap_IndividualExpAndALL_lineAndDot";
includeExp =  {'WHL1000TEST','YWM1000TEST','YWM2000TEST', ...
               'HCP1000TEST', 'HCP2000TEST', 'HCP3000TEST'};
groupStats = GIplot(unsMeansALL, includeExp, 'useDistribution',false, 'ResultType','FARegressionResiduals', ...
       'winSizePix',[0,0,1900,500], 'nrowcol',[1,6], 'ylab', 'L-R FA Resid', ...
       'showXnames',true, 'normResidual', false, 'Group2Individual',false,...
       'calcType','RegressionResiduals', 'nRep',5000, 'onlyBilateral',true,...
       'plotGroupBar', true, ...
       'fnameRoot',fnameRoot,'saveItHere',saveItHere,'saveSvg',false,'savePng',true)
fnameRoot = "IS_ExampleFig3_lineAndDot";
findSpace(groupStats.CI90, ...
       'winSizePix',[0,0,1900,500], 'nrowcol',[1,6], 'ylab', 'Right FA ', 'showLineForm',true, ...
       'xlab', 'Left FA ', 'LineColor',[.1 .1 .1 ],'LineStyle','-', ...
       'useStdAllBS','BS', 'howManySD',[2,1] ,'CIperc',[95,68], ...
       'plotSubj',true, 'SubjValues', SubjValues, ...
       'fnameRoot',fnameRoot,'saveItHere',saveItHere,'saveSvg',true,'savePng',true) 

%% Fig 9 ("GI_RawResiduals_FA_ALL_NoBootstrap_Group2IndividualExpAndALL") 
fnameRoot = "GI_RawResiduals_FA_ALL_NoBootstrap_Group2IndividualExpAndALL";
includeExp =  {'WHL1000TEST','YWM1000TEST','YWM2000TEST', ...
               'HCP1000TEST','HCP1000RETEST', ...
               'HCP2000TEST','HCP2000RETEST', ...
               'HCP3000TEST','HCP3000RETEST'};
groupStats = GIplot(unsMeansALL, includeExp, 'useDistribution',false, 'ResultType','FARegressionResiduals', ...
       'winSizePix',[0,0,1900,500], 'nrowcol',[1,6], 'ylab', 'L-R FA Resid', ...
       'showXnames',true, 'normResidual', false, 'Group2Individual',true,...
       'calcType','RegressionResiduals', 'nRep',5000, 'onlyBilateral',true,...
       'plotGroupBar', true, ...
       'fnameRoot',fnameRoot,'saveItHere',saveItHere,'saveSvg',false,'savePng',false)
close(gcf)

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
