%% INSTRUCTIONS


% Download or clone the repository
%       !git clone https://github.com/garikoitz/paper-reproducibility.git
% 
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
    clear all; close all; clc;
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

% CompRepCheck (tests I did to see how minDist and minLen and using new wmMask)
%{
dt = load(fullfile(paperReprPath,'local','cache', ...
              sprintf('CompRepTest_FSL_%s.mat',measure)));
dt = dt.dt;
% Rename project names
dt.Proj = renamecats(dt.Proj,'HCP_preproc','HCP');
dt.Proj = renamecats(dt.Proj,'PRATIK','WHL');
dt.Proj = renamecats(dt.Proj,'Weston Havens','YWM');
% Change the names of the TRT fields
dt.TRT(dt.Proj=="WHL" & dt.AnalysisMD.mrtrix_useACT==0)   = 'TEST_FSL';
dt.TRT(dt.Proj=="YWM" & dt.AnalysisMD.mrtrix_useACT==0)   = 'TEST_FSL';
dt.TRT(dt.Proj=="HCP" & dt.AnalysisMD.mrtrix_useACT==0 & dt.TRT=="TEST")   = 'TEST_FSL';
dt.TRT(dt.Proj=="HCP" & dt.AnalysisMD.mrtrix_useACT==0 & dt.TRT=="RETEST")   = 'RETEST_FSL';
dt.TRT = removecats(dt.TRT);


% Just maintain the same subjects I am using in the test
DT = DT(ismember(DT.SubjID,unique(dt.SubjID)),:);

% Add the missing fields
DT.AnalysisMD.mrtrix_useACT   = false(height(DT),1);
DT.AnalysisMD.mrtrix_autolmax = true(height(DT),1);
DT.AnalysisMD.mrtrix_lmax     = 8 * ones(height(DT),1);



dT = [DT;dt];
tnames = unique(dT.Struct);
for nt=1:20
    nt=10
    % subplot(4,5,nt)
    tname  = tnames(nt);
    before = dT(dT.Struct==tname & ...
                dT.TRT=="TEST" & ...
                dT.AcquMD.scanbValue==1000,:);
    before = before(before.SubjID=="042_CB",:);
    after  = dT(dT.Struct==tname & ...
                dT.TRT=="TEST_FSL" & ...
                dT.AcquMD.scanbValue==1000,:);
    after  = after(after.SubjID=="042_CB",:);
    plot(before{:,'Val'}', '-'); hold on;
    plot( after{:,'Val'}', '--'); hold off;
    % plot(dT{dT.Struct==tname & dT.TRT=="RETEST",'Val'}', '--'); hold off;
    title(string(tname))
end
%} 
% END OF CompRepCheck

% VOLUME CHECK (instead of FA)
%{
VOL_CR = load(fullfile(paperReprPath, 'local','cache','ComputationalReproducibility_volume.mat'));
measure2 = 'volume';
DV = VOL_CR.dt;
% Rename project names
DV.Proj = renamecats(DV.Proj,'HCP_preproc','HCP');
DV.Proj = renamecats(DV.Proj,'PRATIK','WHL');
DV.Proj = renamecats(DV.Proj,'Weston Havens','YWM');
% Change the names of the TRT fields
DV.TRT(DV.Proj=="WHL") = 'TEST';
DV.TRT(DV.Proj=="YWM")   = 'TEST';
% Clean
allstructs = unique(DV.Struct);
% for ns=1:length(unique(DV.Struct))
%     subplot(5,5,ns)
%     tn  = allstructs(ns);
%     COM1000 = DV(DV.SubjID=='042_CB' & DV.Struct==tn & DV.AcquMD.scanbValue==1000,:);
%     scatter(COM1000.Val(1,:),COM1000.Val(2,:));hold on; identityLine(gca);
%     title(string(tn));xlabel('TEST1');ylabel('TEST2');
% end
% I am going to delete the second instance of the subject

DV.ROW   = [1:height(DV)]';
toDelete = [];
for ns=1:length(unique(DV.Struct))
    tn  = allstructs(ns);
    chooseOne = DV.ROW(DV.SubjID=='042_CB' & DV.Struct==tn & DV.AcquMD.scanbValue==1000 ,:);
    toDelete = [toDelete, chooseOne(2)];
end
DV(toDelete,:) = [];
%}
% END OF VOLUME CHECK

% Check the older WH dataset
%{
A       = load('~/Downloads/WH_database_full_metadata.mat');
fadt    = cell2table(A.afq.vals.fa);
fadt.Properties.VariableNames = strrep(A.afq.fgnames,' ','');
tracts  = [{'CingulumCingulate'}, {'Arcuate'}, {'IFOF'}, {'ILF'}, {'Uncinate'}, {'Corticospinal'}];
fadt    = fadt(:, contains(fadt.Properties.VariableNames,tracts));
FADT    = table();
meanFADT= table();
for nf=1:width(fadt)
    tname = fadt.Properties.VariableNames{nf};
    FADT.(tname) = cell2table(num2cell(fadt{:,nf}{1}, 2));
    % FADT.(strcat('mean',tname)) = mean(FADT.(tname){:,:},2);
    meanFADT.(tname) = mean(FADT.(tname){:,:},2);
end
meanFADT.Properties.VariableNames = strrep(meanFADT.Properties.VariableNames,'mean','');
% See what's going on with the nans
onans = @(x) nnz(isnan(x));
nansByProject = varfun(onans,meanFADT);
nansByProject.Properties.VariableNames = strrep(nansByProject.Properties.VariableNames,'Fun_','')





%}
% End of check the older WH dataset


% legend(cellstr(allstructs),'location','southeast')
% identityLine(gca);

% Check if everything looks ok
% summary(dt)
cmapname = 'copper';
if ~exist(saveItHere); mkdir(saveItHere); end

%% Create working tables
% Select only some of the tracts
    % tracts = [{'IFOF'}, {'ILF'}, {'SLF'}, {'Arcuate'}];
    % tracts = [{'Thalamic'}, {'ILF'}, {'Uncinate'}, {'Corticospinal'}];
    % Wahl tracts
    tracts = [{'CingulumCingulate'}, {'Arcuate'}, {'IFOF'}, {'ILF'}, {'Uncinate'}, {'Corticospinal'}];

% Do the filtering and obtain the unstacked tract profiles
[dtt, unsProf, unsMeans, indivTracts, pairedTracts] = dr_createWorkingTables(DT, ...
                                 'sliceBasedOn',           {'Proj', 'SHELL'}, ...
                                 'SHELL',                  {'1000', '2000', '3000'}, ...
                                 'removeTractsContaining', {'Callosum'}, ... 
                                 'useThisTractsOnly',       tracts, ...
                                 'createTractPairs',        true, ... 
                                 'TRT',                    'TEST', ...
                                 'AGE',                     [14,58], ...
                                 'GENDER',                 {'male', 'female'});
                             
[dtrt, unsProfrt, unsMeansrt] = dr_createWorkingTables(DT, ...
                                 'sliceBasedOn',           {'Proj', 'SHELL'}, ...
                                 'SHELL',                  {'1000', '2000', '3000'}, ...
                                 'removeTractsContaining', {'Callosum'}, ... 
                                 'useThisTractsOnly',       tracts, ...
                                 'createTractPairs',        true, ...
                                 'TRT',                    'RETEST', ...
                                 'AGE',                     [14,58], ...
                                 'GENDER',                 {'male', 'female'});
                             
[dtTRT, unsProfTRT, unsMeansTRT] = dr_createWorkingTables(DT(DT.Proj=='HCP',:), ...
                                 'sliceBasedOn',           {'Proj', 'SHELL','TRT'}, ...
                                 'SHELL',                  {'1000', '2000', '3000'}, ...
                                 'removeTractsContaining', {'Callosum'}, ... 
                                 'useThisTractsOnly',       tracts, ...
                                 'createTractPairs',        true, ...
                                 'TRT',                    'notTRAIN', ...
                                 'AGE',                     [14,58], ...
                                 'GENDER',                 {'male', 'female'});
                             
                             
% See what' going on with the nans
fprintf('Total number of NaN: %d, total values: %d, %% of NaNs: %2.2f%%\n', sum(isnan(dtt.meanVal)), length(dtt.meanVal), 100*(sum(isnan(dtt.meanVal))/length(dtt.meanVal)))
fprintf('\n\nSame thing but per project: \n')
onans = @(x) nnz(isnan(x));
nansByProject = varfun(onans,dtt,'GroupingVariables',{'Proj','SHELL'},'InputVariables',{'meanVal'});
nansByProject.Percentage = 100*(nansByProject.Fun_meanVal./nansByProject.GroupCount)
fprintf('\n\nSame thing but per tract and project: \n')
nansByTractProject = varfun(onans,dtt,'GroupingVariables',{'Proj','SHELL','Struct'},'InputVariables',{'meanVal'});
nansByTractProject.Percentage = 100*(nansByTractProject.Fun_meanVal./nansByTractProject.GroupCount)
                             
% Do the same for the VOLUMES
%{
[dvv, dvunsProf, dvunsMeans, dvindivTracts, dvpairedTracts] = dr_createWorkingTables(DV, ...
                                 'sliceBasedOn',           {'Proj', 'SHELL'}, ...
                                 'SHELL',                  {'1000', '2000', '3000'}, ...
                                 'removeTractsContaining', {'Callosum'}, ... 
                                 'useThisTractsOnly',       tracts, ...
                                 'createTractPairs',        true, ... 
                                 'TRT',                    'TEST', ...
                                 'AGE',                     [14,58], ...
                                 'GENDER',                 {'male', 'female'});
                             
[dvtrt, dvunsProfrt, dvunsMeansrt] = dr_createWorkingTables(DV, ...
                                 'sliceBasedOn',           {'Proj', 'SHELL'}, ...
                                 'SHELL',                  {'1000', '2000', '3000'}, ...
                                 'removeTractsContaining', {'Callosum'}, ... 
                                 'useThisTractsOnly',       tracts, ...
                                 'createTractPairs',        true, ...
                                 'TRT',                    'RETEST', ...
                                 'AGE',                     [14,58], ...
                                 'GENDER',                 {'male', 'female'});
                             
[dvTRT, dvunsProfTRT, dvunsMeansTRT] = dr_createWorkingTables(DV(DV.Proj=='HCP',:), ...
                                 'sliceBasedOn',           {'Proj', 'SHELL','TRT'}, ...
                                 'SHELL',                  {'1000', '2000', '3000'}, ...
                                 'removeTractsContaining', {'Callosum'}, ... 
                                 'useThisTractsOnly',       tracts, ...
                                 'createTractPairs',        true, ...
                                 'TRT',                    'notTRAIN', ...
                                 'AGE',                     [14,58], ...
                                 'GENDER',                 {'male', 'female'});
                           
                             
% Create version of the FA corrected by the volumes. Same names, we will use it
% calculate the plots and then deactivate this part.
normFaWithVolumeAcrossBvalues = false;
normFaWithVolumeIndepBvalues  = false;
if normFaWithVolumeAcrossBvalues
    % I need unsMeans and unsMeansrt to have the FA corrected by volume and by
    % project, but NOT by b value, at first. 
    tracts      = unsMeans.Properties.VariableNames(contains(unsMeans.Properties.VariableNames,{'Left','Right'}));
    projects    = unique(unsMeans.Proj);
    
    
%     isequal(unique(FAs.SubjID), unique(VOLs.SubjID))
%     FAs         = sortrows(FA , 'SubjID');   
%     VOLs        = sortrows(VOL, 'SubjID');
%     bvalues     = {1000,2000,3000};
%     bluecolors{1000}=[0         0    0.5625];
%     bluecolors{2000}=[0     0.447     0.741];
%     bluecolors{3000}=[0    0.8125         1];
    
    for nt = 1: length(tracts)
        for np = 1: length(projects)
            % unsMeans
            favalues  = unsMeans(unsMeans.Proj==projects(np),{tracts{nt},'SubjID'});
            volvalues = dvunsMeans(dvunsMeans.Proj==projects(np),{tracts{nt},'SubjID'});
            if ~isequal(favalues(:,2), volvalues(:,2)); error('Subject names/order do not coincide'); end
            favalues  = favalues{:,1}; volvalues = volvalues{:,1};
            favaluesN = dr_normRemoveEffect(favalues,volvalues);
            unsMeans{unsMeans.Proj==projects(np),{tracts{nt}}} = favaluesN;
            
            % unsMeansrt
            if strcmp(projects(np),'HCP')
                favalues  = unsMeansrt(unsMeansrt.Proj==projects(np),{tracts{nt},'SubjID'});
                volvalues = dvunsMeansrt(dvunsMeansrt.Proj==projects(np),{tracts{nt},'SubjID'});
                if ~isequal(favalues(:,2), volvalues(:,2)); error('Subject names/order do not coincide'); end
                favalues  = favalues{:,1}; volvalues = volvalues{:,1};
                favaluesN = dr_normRemoveEffect(favalues,volvalues);
                unsMeansrt{unsMeansrt.Proj==projects(np),{tracts{nt}}} = favaluesN;
            end
            
        end
    end
end
if normFaWithVolumeIndepBvalues
    % I need unsMeans and unsMeansrt to have the FA corrected by volume and by
    % project, but NOT by b value, at first. 
    tracts      = unsMeans.Properties.VariableNames(contains(unsMeans.Properties.VariableNames,{'Left','Right'}));
    projects    = unique(unsMeans.Proj);

    
%     isequal(unique(FAs.SubjID), unique(VOLs.SubjID))
%     FAs         = sortrows(FA , 'SubjID');   
%     VOLs        = sortrows(VOL, 'SubjID');
%     bvalues     = {1000,2000,3000};
%     bluecolors{1000}=[0         0    0.5625];
%     bluecolors{2000}=[0     0.447     0.741];
%     bluecolors{3000}=[0    0.8125         1];
    
    for nt = 1: length(tracts)
        for np = 1: length(projects)
            shells      = unique(unsMeans.SHELL(unsMeans.Proj==projects(np)));
            for ns = 1: length(shells)
                % unsMeans
                favalues  = unsMeans(unsMeans.Proj==projects(np) & ...
                                     unsMeans.SHELL==shells(ns), ...
                                     {tracts{nt},'SubjID'});
                volvalues = dvunsMeans(dvunsMeans.Proj==projects(np) & ...
                                      dvunsMeans.SHELL==shells(ns), ...
                                      {tracts{nt},'SubjID'});
                if ~isequal(favalues(:,2), volvalues(:,2)); error('Subject names/order do not coincide'); end
                favalues  = favalues{:,1}; volvalues = volvalues{:,1};
                favaluesN = dr_normRemoveEffect(favalues,volvalues);
                unsMeans{unsMeans.Proj==projects(np) & unsMeans.SHELL==shells(ns),{tracts{nt}}} = favaluesN;

                % unsMeansrt
                if strcmp(projects(np),'HCP')
                    favalues  = unsMeansrt(unsMeansrt.Proj==projects(np) & ...
                                     unsMeans.SHELL==shells(ns), ...
                                     {tracts{nt},'SubjID'});
                    volvalues = dvunsMeansrt(dvunsMeansrt.Proj==projects(np) & ...
                                      dvunsMeans.SHELL==shells(ns), ...
                                      {tracts{nt},'SubjID'});
                    if ~isequal(favalues(:,2), volvalues(:,2)); error('Subject names/order do not coincide'); end
                    favalues  = favalues{:,1}; volvalues = volvalues{:,1};
                    favaluesN = dr_normRemoveEffect(favalues,volvalues);
                    unsMeansrt{unsMeansrt.Proj==projects(np) & unsMeansrt.SHELL==shells(ns),{tracts{nt}}} = favaluesN;
                end
            end  
        end
    end
end  
%} 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ALL CALCULATIONS BELOW NEED TO BE INDEPENDENT AND STAND-ALONE                             
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                              
                             
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

%% FA PROFILES 
% VISUALIZE RESULTS

% TEST
dt        = dtt;
cats      = categories(dt.SliceCats);
catcolors = unique(dt.SliceCatsRGB);
tracts    = unique(dt.Struct);
% Wahl Order
tracts    = tracts([3,4,  11,12,  5,6,   7,8,   9,10,   1,2]);
mrvNewGraphWin('FA profiles'); 
for nt = 1: length(tracts)
    subplot(3,4,nt)
    tn = tracts(nt);
    a = [];
    for nc=1:length(cats)
        cat      = cats{nc}; 
        catColors= dt{dt.SliceCats==cat,'SliceCatsRGB'}; catColor=catColors{1,:};
        values   = mean(dt{dt.Struct==tn & dt.SliceCats==cat,'Val'},1,'omitnan');
        sdvalues = std(dt{dt.Struct==tn & dt.SliceCats==cat,'Val'},1,'omitnan');
        upper    = values + sdvalues; lower = values - sdvalues;
        a = [a; plot([1:100],values,'color',catColor,'linewidth',2)]; hold on;
        jbfill([1:100],upper,lower,catColor,catColor,0,0.1);
    end
    if (nt > 6 && nt < 11)
        legend(a, strrep(cats,'_','\_'), 'location','northwest'); 
    else
        legend(a, strrep(cats,'_','\_'), 'location','southwest'); 
    end
    ylabel('FA', 'FontWeight','bold'); % xlabel('Divisions')
    set(gca,'xtick',[])
    ylim([0.1, 0.75]); yticks([0.2,.4,.6])
    set(gca,'FontSize',18)
    title(sprintf('%s',tn))
end
set(gcf,'color','w');
suptitle({'FA profiles.', 'Group average (1 SD)'})
saveas(gcf,fullfile(saveItHere, 'FA_Profiles.svg'),'svg');
saveas(gcf,fullfile(saveItHere, 'FA_Profiles.png'),'png');
close(gcf)

% HCP TEST-RETEST
dt = dtTRT;
cats      = categories(dt.SliceCats);
catcolors = unique(dt.SliceCatsRGB);
catcolors = catcolors([1,1,2,2,3,3],:);
linestyles= {':','-',':','-',':','-'};
tracts    = unique(dt.Struct);
% Wahl Order
tracts    = tracts([3,4,  11,12,  5,6,   7,8,   9,10,   1,2]);
mrvNewGraphWin('FA profiles'); 
for nt = 1: length(tracts)
    subplot(3,4,nt)
    tn = tracts(nt);
    a = [];myCat = {};
    for nc=1:length(cats)
        cat      = cats{nc}; catColor = catcolors{nc,:}; lineStyle = linestyles{nc}; 
        values   = mean(dt{dt.Struct==tn & dt.SliceCats==cat,'Val'},1,'omitnan');
        sdvalues = std(dt{dt.Struct==tn & dt.SliceCats==cat,'Val'},1,'omitnan');
        upper    = values + sdvalues; lower = values - sdvalues;
        a = [a; plot([1:100],values,'color',catColor,'linewidth',2,'LineStyle',lineStyle)]; hold on; 
        myCat = {myCat; cat};
        jbfill([1:100],upper,lower,catColor,catColor,0,0.1);
    end

    if (nt > 6 && nt < 11)
        legend(a, strrep(cats,'_','\_'), 'location','northwest'); 
    else
        legend(a, strrep(cats,'_','\_'), 'location','southwest'); 
    end
    ylabel('FA', 'FontWeight','bold'); % xlabel('Divisions')
    set(gca,'xtick',[])
    ylim([0.1, 0.75]); yticks([0.2,.4,.6])
    set(gca,'FontSize',18)
    title(sprintf('%s',tn))
end
set(gcf,'color','w');
suptitle({'HCP TEST-RETEST FA profiles.', 'Group average (1 SD)'})
saveas(gcf,fullfile(saveItHere, 'FA_Profiles_HCP_TEST_RETEST.svg'),'svg');
saveas(gcf,fullfile(saveItHere, 'FA_Profiles_HCP_TEST_RETEST.png'),'png');
close(gcf)

%% TEST: MEAN FA BARS/GENERALIZATION PLOTS
% Duplicate the first table in Wahl 2010
omean   = @(x) mean(x,'omitnan');
ostd    = @(x) std(x,'omitnan');
left_color = [0 0 .55]; right_color = [0.3 0.3 0.3];


% ALL
tmp     = varfun(ostd,dtt,'GroupingVariables',{'Struct','SliceCats'},'InputVariables',{'meanVal'});
replRes = varfun(omean,dtt,'GroupingVariables',{'Struct','SliceCats'},'InputVariables',{'meanVal'});
replRes.Properties.VariableNames{'Fun_meanVal'} = 'Mean';
replRes.SD  = tmp.Fun_meanVal;
replRes     = replRes(:,[1,2,4,5]);
% Unstack it to show it as consecutive tables. Reorder to replicate
replRes     = unstack(replRes,{'Mean','SD'},{'SliceCats'});
replRes     = replRes([3,4,  11,12,   5,6,   7,8,   9,10,   1,2],:);
replRes.Mean_WHLorig = [0.565,0.522,0.517,0.490,0.549,0.534,0.510,0.497,0.489,0.470,0.587,0.575]';
replRes.SD_WHLorig   = [0.028,0.028,0.023,0.03,0.026,0.024,0.028,0.026,0.022,0.024,0.023,0.021]';
replResSorted         = replRes(:,[1, 14,15,  5,11,   6,12,   7,13,   2,8,  3,9,   4,10]);

% Plot barplots with error bars
dt = dtt;
cats      = categories(dt.SliceCats);
catcolors = unique(dt.SliceCatsRGB);
catcolors.cat = categorical(cats);
catcolors     = [table(0,0,0,{'WHLorig'},'VariableNames',{'R','G','B','cat'}); ...
                 catcolors];
catcolors     = catcolors([1,5,6,7,2,3,4],:);

bigfig = mrvNewGraphWin('Mean FA per project'); 
for nt = 1: height(replResSorted)
    set(bigfig,'defaultAxesColorOrder',[left_color; right_color]);
    subplot(3,4,nt)
    tn = string(replResSorted{nt,'Struct'});
    Names = replResSorted.Properties.VariableNames(contains(replResSorted.Properties.VariableNames, 'Mean_'));
    Names = strrep(Names,'Mean_','');
    means = replResSorted{nt,contains(replResSorted.Properties.VariableNames, 'Mean_')};
    stds  = replResSorted{nt,contains(replResSorted.Properties.VariableNames, 'SD_')};
    CV    = 100*(stds ./ means);
    inBarColor = {'w','k','w','k','w','k','k'};
    yyaxis right; ylim([0 20]); ylabel('CoV (%)');
    plot([1:length(CV)],CV,'--','Color',right_color,'linewidth',2);
    for nc=1:length(means)
        yyaxis left; ylabel('Mean FA'); ylim([0, .8])
        bar(nc, means(nc), 'FaceColor', catcolors{nc,{'R','G','B'}},'EdgeColor','none','BarWidth',.9); hold on;
        set(gca,'xtick',[]); set(gca,'xtickLabel',[]);
        % T = text(nc, -0.05, string(catcolors.cat(nc)),'Color',[.0 .0 .0 ], 'Rotation', 45, 'HorizontalAlignment', 'Right');
        d  = errorbar(nc,means(nc),stds(nc),'.', 'color',[.5 .5 .5 ]); d.LineWidth = 2;
        T1 = text(nc, -0.05, string(catcolors.cat(nc)),'Color',[.0 .0 .0 ], 'Rotation', 45, 'HorizontalAlignment', 'Right');
        % T2 = text(nc, 0.05, sprintf('CoV=%1.2f%%',CV(nc)),'Color',inBarColor{nc}, 'Rotation', 90, 'HorizontalAlignment', 'Left','fontweight','bold','fontsize',12);
        yyaxis right; ylim([0 20]);ylabel('CoV (%)');
        plot(nc,CV(nc),'o', 'MarkerEdgeColor',right_color,'MarkerFaceColor','w','MarkerSize',8);hold on;
    end
    title(sprintf('%s',tn))
end
set(bigfig,'color','w');
suptitle({'Mean FA per project.'})
saveas(bigfig,fullfile(saveItHere, 'FA_Means_withCoVdots.png'),'png');
close(bigfig)










% Prepare the data to be able to plot later on. 
% This has bee copied from below, were the correlations were calculated
recalculate  = true;
nRep = 500;
CIrangeOrVals = 90;
useDistribution = true;
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

%% HCP TEST-RETEST: MEAN FA BARS/GENERALIZATION PLOTS
omean   = @(x) mean(x,'omitnan');
ostd    = @(x) std(x,'omitnan');
left_color = [0 0 .55]; right_color = [0.3 0.3 0.3];

tmp     = varfun(ostd,dtTRT,'GroupingVariables',{'Struct','SliceCats'},'InputVariables',{'meanVal'});
replRes = varfun(omean,dtTRT,'GroupingVariables',{'Struct','SliceCats'},'InputVariables',{'meanVal'});
replRes.Properties.VariableNames{'Fun_meanVal'} = 'Mean';
replRes.SD  = tmp.Fun_meanVal;
replRes     = replRes(:,[1,2,4,5]);
% Unstack it to show it as consecutive tables. Reorder to replicate
replRes     = unstack(replRes,{'Mean','SD'},{'SliceCats'});
replRes     = replRes([3,4,  11,12,   5,6,   7,8,   9,10,   1,2],:);

% This code is with WHLorig at the beginning
replRes.Mean_WHLorig = [0.565,0.522,0.517,0.490,0.549,0.534,0.510,0.497,0.489,0.470,0.587,0.575]';
replRes.SD_WHLorig   = [0.028,0.028,0.023,0.03,0.026,0.024,0.028,0.026,0.022,0.024,0.023,0.021]';
% replResSortedTRT      = replRes(:,[1, 14,15,  3,9,   2,8,   5,11,   4,10,  7,13,   6,12]);
% This code is without WHLorig at the beginning
% replRes.Mean_WHLorig = [0.565,0.522,0.517,0.490,0.549,0.534,0.510,0.497,0.489,0.470,0.587,0.575]';
% replRes.SD_WHLorig   = [0.028,0.028,0.023,0.03,0.026,0.024,0.028,0.026,0.022,0.024,0.023,0.021]';
replResSortedTRT      = replRes(:,[1, 3,9,   2,8,   5,11,   4,10,  7,13,   6,12]);

% Plot barplots with error bars
dt = dtTRT;
cats      = categories(dt.SliceCats);
catcolors = unique(dt.SliceCatsRGB);
% Do not add black for WHLorig, and remove the first column then
% catcolors     = [table(0,0,0,{'WHLorig'},'VariableNames',{'R','G','B','cat'}); ...
%                  catcolors];
% catcolors     = catcolors([1,2,2,3,3,4,4],:);
catcolors     = catcolors([1,1,2,2,3,3],:);

bigfig = mrvNewGraphWin('Mean FA per project'); 
for nt = 1: height(replResSortedTRT)
    % left_color = [0 0 .55]; right_color = [0.2 0.2 0.2];
    set(bigfig,'defaultAxesColorOrder',[left_color; right_color]);
    subplot(3,4,nt);
    tn = string(replResSortedTRT{nt,'Struct'});
    Names = replResSortedTRT.Properties.VariableNames(contains(replResSortedTRT.Properties.VariableNames, 'Mean_'));
    Names = strrep(Names,'Mean_','');
    means = replResSortedTRT{nt,contains(replResSortedTRT.Properties.VariableNames, 'Mean_')};
    stds  = replResSortedTRT{nt,contains(replResSortedTRT.Properties.VariableNames, 'SD_')};
    CV    = 100*(stds ./ means);
    inBarColor = {'w','w','k','k','k','k'};
    yyaxis right; ylim([0 20]); ylabel('CoV (%)');
    % plot([1:length(CV)],CV,'--','Color',right_color,'linewidth',2);
    for nc=1:length(means)
        yyaxis left; ylabel('Mean FA');ylim([0, .8])
        bar(nc, means(nc), 'FaceColor', catcolors{nc,{'R','G','B'}},'EdgeColor','none','BarWidth',.9); hold on;
        set(gca,'xtick',[]); set(gca,'xtickLabel',[]);
        T = text(nc, -0.03, Names{nc},'Color',[.0 .0 .0 ], 'Rotation', 45, 'HorizontalAlignment', 'Right');
        d = errorbar(nc,means(nc),stds(nc),'.', 'color',[.5 .5 .5 ]); d.LineWidth = 2;
        % T2 = text(nc, 0.05, sprintf('CoV=%1.2f%%',CV(nc)),'Color',inBarColor{nc}, 'Rotation', 90, 'HorizontalAlignment', 'Left','fontweight','bold','fontsize',12);
        yyaxis right; ylim([0 20]);ylabel('CoV (%)');
        plot(nc,CV(nc),'o', 'MarkerEdgeColor',right_color,'MarkerFaceColor','w','MarkerSize',8);hold on;
    end
    title(sprintf('%s',tn))
end
set(bigfig,'color','w');
suptitle({'Mean FA per project.', 'HCP TEST-RETEST'})
saveas(bigfig,fullfile(saveItHere, 'FA_Means_withCoVdots_HCP_TRT.png'),'png');
close(bigfig)


mrvNewGraphWin('SD FA per project'); 
for nt = 1: height(replResSortedTRT)
    subplot(3,4,nt)
    tn = string(replResSortedTRT{nt,'Struct'});
    Names = replResSortedTRT.Properties.VariableNames(contains(replResSortedTRT.Properties.VariableNames, 'Mean_'));
    Names = strrep(Names,'Mean_','');
    means = replResSortedTRT{nt,contains(replResSortedTRT.Properties.VariableNames, 'Mean_')};
    stds  = replResSortedTRT{nt,contains(replResSortedTRT.Properties.VariableNames, 'SD_')};
    for nc=1:length(means)
        bar(nc, stds(nc), 'FaceColor', catcolors{nc,{'R','G','B'}},'EdgeColor','none','BarWidth',.9); hold on;
        set(gca,'xtick',[]); set(gca,'xtickLabel',[]);
        T = text(nc, -0.03, Names{nc},'Color',[.0 .0 .0 ], 'Rotation', 45, 'HorizontalAlignment', 'Right');
        % d = errorbar(nc,means(nc),stds(nc),'.', 'color',[.5 .5 .5 ]); d.LineWidth = 2;
    end
    ylabel('SD of FA');
    ylim([0, .05])
    title(sprintf('%s',tn))
end
set(gcf,'color','w');
suptitle({'SD FA per project.', 'HCP TEST-RETEST'})
saveas(gcf,fullfile(saveItHere, 'FA_SD_HCP_TRT.png'),'png');
close(gcf)


% TRT scatterplots
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

mrvNewGraphWin('Mean FA TEST-RETEST Scatterplots'); 
allrmse1000  = []; allrmse2000  = []; allrmse3000  = [];
allrrmse1000 = []; allrrmse2000 = []; allrrmse3000 = [];
allicc1000   = []; allicc2000   = []; allicc3000   = [];
allcov1000   = []; allcov2000   = []; allcov3000   = [];
Tstd1000     = []; Tstd2000     = []; Tstd3000     = [];
RTstd1000    = []; RTstd2000    = []; RTstd3000    = [];
for nt = 1: height(replResSortedTRT)
    subplot(3,4,nt)
    tn = string(replResSortedTRT{nt,'Struct'});
    b1000T    = unsMeansTRT(unsMeansTRT.TRT=='TEST' & unsMeansTRT.SHELL=='1000',{char(tn),'SubjID'});
    b1000RT   = unsMeansTRT(unsMeansTRT.TRT=='RETEST' & unsMeansTRT.SHELL=='1000', {char(tn),'SubjID'});
    b1000T=sortrows(b1000T,'SubjID'); b1000RT=sortrows(b1000RT,'SubjID');
    if isequal(b1000T.SubjID,b1000RT.SubjID)
        b1000T = b1000T.(tn); b1000RT = b1000RT.(tn);
    else
        error('The test-retest is not comparing the same subjects')
    end

    b2000T    = unsMeansTRT(unsMeansTRT.TRT=='TEST' & unsMeansTRT.SHELL=='2000', {char(tn),'SubjID'});
    b2000RT   = unsMeansTRT(unsMeansTRT.TRT=='RETEST' & unsMeansTRT.SHELL=='2000', {char(tn),'SubjID'});
    b2000T=sortrows(b2000T,'SubjID'); b2000RT=sortrows(b2000RT,'SubjID');
    if isequal(b2000T.SubjID,b2000RT.SubjID)
        b2000T = b2000T.(tn); b2000RT = b2000RT.(tn);
    else
        error('The test-retest is not comparing the same subjects')
    end
    
    b3000T    = unsMeansTRT(unsMeansTRT.TRT=='TEST' & unsMeansTRT.SHELL=='3000', {char(tn),'SubjID'});
    b3000RT   = unsMeansTRT(unsMeansTRT.TRT=='RETEST' & unsMeansTRT.SHELL=='3000', {char(tn),'SubjID'});
    b3000T=sortrows(b3000T,'SubjID'); b3000RT=sortrows(b3000RT,'SubjID');
    if isequal(b3000T.SubjID,b3000RT.SubjID)
        b3000T = b3000T.(tn); b3000RT = b3000RT.(tn);
    else
        error('The test-retest is not comparing the same subjects')
    end
    
    scatter(b1000T,b1000RT,20,catcolors{1,:},'filled');hold on;
    scatter(b2000T,b2000RT,18,catcolors{3,:});
    scatter(b3000T,b3000RT,20,catcolors{5,:},'filled');
    xlim([0.14, 0.66]);ylim([0.14, 0.66]);
    xticks([0.2:0.1:0.6]);yticks([0.2:0.1:0.6]);
    identityLine(gca);axis equal
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
    T1 = text(.20, 0.625, sprintf('CoV (b=1000)=%1.2f%%',CoVm1),'Color',catcolors{1,:}, 'Rotation', 0, 'HorizontalAlignment', 'Left','fontweight','normal','fontsize',16);
    T2 = text(.15, 0.575, sprintf('CoV (b=2000)=%1.2f%%',CoVm2),'Color',catcolors{3,:}, 'Rotation', 0, 'HorizontalAlignment', 'Left','fontweight','normal','fontsize',16);
    T3 = text(.1, 0.525, sprintf('CoV (b=3000)=%1.2f%%',CoVm3),'Color',catcolors{5,:}, 'Rotation', 0, 'HorizontalAlignment', 'Left','fontweight','normal','fontsize',16);
    xlabel('FA (TEST)','FontWeight','bold');
    ylabel('FA (RETEST)','FontWeight','bold');
    set(gca,'FontSize',18)
    grid off

end
set(gcf,'color','w');
suptitle({'Mean FA scatterplots', 'HCP TEST-RETEST'})
saveas(gcf,fullfile(saveItHere, 'FA_Scatterplots_rmse_CoV_HCP_TRT.png'),'png');
saveas(gcf,fullfile(saveItHere, 'FA_Scatterplots_rmse_CoV_HCP_TRT.svg'),'svg');
close(gcf)   
    
% Create the table with the avsvgges
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

% Plot to illustrate the point in the discussion
% The point is that I do not have a point: sd is no lower.
%{
correlation1000 = [];correlation2000 = [];correlation3000 = [];
std1000 = [];std2000 = [];std3000 = [];
for nt=[1:2:11]
    subplot(3,4,nt)
    tn = string(replResSortedTRT{nt,'Struct'});
    L1000    = unsMeansTRT{unsMeansTRT.TRT=='TEST' & unsMeansTRT.SHELL=='1000',{char(tn)}};
    L2000    = unsMeansTRT{unsMeansTRT.TRT=='TEST' & unsMeansTRT.SHELL=='2000',{char(tn)}};
    L3000    = unsMeansTRT{unsMeansTRT.TRT=='TEST' & unsMeansTRT.SHELL=='3000',{char(tn)}};
    tn = string(replResSortedTRT{nt+1,'Struct'});
    R1000    = unsMeansTRT{unsMeansTRT.TRT=='TEST' & unsMeansTRT.SHELL=='1000',{char(tn)}};
    R2000    = unsMeansTRT{unsMeansTRT.TRT=='TEST' & unsMeansTRT.SHELL=='2000',{char(tn)}};
    R3000    = unsMeansTRT{unsMeansTRT.TRT=='TEST' & unsMeansTRT.SHELL=='3000',{char(tn)}};
    correlation1000 = [correlation1000, corr(L1000,R1000)];
    std1000 = [std1000, (std(L1000)+std(R1000))/2];
    correlation2000 = [correlation2000, corr(L2000,R2000)];
    std2000 = [std2000, (std(L2000)+std(R2000))/2];
    correlation3000 = [correlation3000, corr(L3000,R3000)];
    std3000 = [std3000, (std(L3000)+std(R3000))/2];
    
    scatter((std(L1000)+std(R1000))/2,corr(L1000,R1000),20,catcolors{1,:},'filled');hold on;
    scatter((std(L2000)+std(R2000))/2,corr(L2000,R2000),18,catcolors{3,:});
    scatter((std(L3000)+std(R3000))/2,corr(L3000,R3000),20,catcolors{5,:},'filled');
    
%     scatter(std1000,correlation1000,20,catcolors{1,:},'filled');hold on;
%     scatter(std2000,correlation2000,18,catcolors{3,:});
%     scatter(std3000,correlation3000,20,catcolors{5,:},'filled');

end
%}

% Same plot with volumes instead of FA
%{
omean   = @(x) mean(x,'omitnan');
ostd    = @(x) std(x,'omitnan');

% ALL
dvtmp     = varfun(ostd,dvv,'GroupingVariables',{'Struct','SliceCats'},'InputVariables',{'meanVal'});
dvreplRes = varfun(omean,dvv,'GroupingVariables',{'Struct','SliceCats'},'InputVariables',{'meanVal'});
dvreplRes.Properties.VariableNames{'Fun_meanVal'} = 'Mean';
dvreplRes.SD  = dvtmp.Fun_meanVal;
dvreplRes     = dvreplRes(:,[1,2,4,5]);
% Unstack it to show it as consecutive tables. Reorder to replicate
dvreplRes     = unstack(dvreplRes,{'Mean','SD'},{'SliceCats'});
dvreplRes     = dvreplRes([3,4,  11,12,   5,6,   7,8,   9,10,   1,2],:);
dvreplResSorted         = dvreplRes(:,[1,  5,11,   6,12,   7,13,   2,8,  3,9,   4,10]);

% Plot barplots with error bars
dv              = dvv;
dvcats          = categories(dv.SliceCats);
dvcatcolors     = unique(dv.SliceCatsRGB);
dvcatcolors.cat = categorical(dvcats);
dvcatcolors     = dvcatcolors([4,5,6,1,2,3],:);

mrvNewGraphWin('Mean VOL per project'); 
for nt = 1: height(dvreplResSorted)
    subplot(3,4,nt)
    tn = string(dvreplResSorted{nt,'Struct'});
    Names = dvreplResSorted.Properties.VariableNames(contains(dvreplResSorted.Properties.VariableNames, 'Mean_'));
    Names = strrep(Names,'Mean_','');
    means = dvreplResSorted{nt,contains(dvreplResSorted.Properties.VariableNames, 'Mean_')};
    stds  = dvreplResSorted{nt,contains(dvreplResSorted.Properties.VariableNames, 'SD_')};
    for nc=1:length(means)
        bar(nc, means(nc), 'FaceColor', dvcatcolors{nc,{'R','G','B'}},'EdgeColor','none','BarWidth',.9); hold on;
        set(gca,'xtick',[]); set(gca,'xtickLabel',[]);
        T = text(nc, -0.1, string(dvcatcolors.cat(nc)),'Color',[.0 .0 .0 ], 'Rotation', 45, 'HorizontalAlignment', 'Right');
        d = errorbar(nc,means(nc),stds(nc),'.', 'color',[.5 .5 .5 ]); d.LineWidth = 2;
    end
    ylabel('Mean FA');
    % ylim([0, .8])
    title(sprintf('%s',tn))
end
set(gcf,'color','w');
suptitle({'Mean VOL per project.'})
saveas(gcf,fullfile(saveItHere, 'VOL_Means.png'),'png');
close(gcf)
%}
% End Volumes

%% VOLUME vs FA
%{
% SCATTERPLOTS
% FA
mrvNewGraphWin('VOL-FA scatterplots'); 
tn = 'LeftThalamicRadiation';
tracts      = unique(DV.Struct);
repetitions = height(DV(DV.Struct==tn,:));
FA          = DT(DT.Proj=='HCP' & DT.TRT=='TEST',:);
VOL         = DV(DV.Proj=='HCP' & DV.TRT=='TEST',:);
bvalues     = {1000,2000,3000};
bluecolors{1000}=[0         0    0.5625];
bluecolors{2000}=[0     0.447     0.741];
bluecolors{3000}=[0    0.8125         1];
for nt = 1: length(tracts)
    subplot(4,5,nt)
    tn = tracts(nt);
    for nb=bvalues    
        falong = mean(FA{FA.Struct==tn & FA.AcquMD.scanbValue==nb{:},'Val'},2,'omitnan');
        volong = mean(VOL{VOL.Struct==tn & VOL.AcquMD.scanbValue==nb{:},'Val'},2,'omitnan');
        if nb{:}==2000
            scatter(volong,falong,30,bluecolors{nb{:}});hold on;
        else
            scatter(volong,falong,30,bluecolors{nb{:}},'filled');hold on;
        end
        xlabel('Volume'); ylabel('FA');
        
        title(sprintf('%s',tn));
    end
    falong = mean(FA{FA.Struct==tn,'Val'},2,'omitnan');
    volong = mean(VOL{VOL.Struct==tn,'Val'},2,'omitnan');
    [rho,pval] = corr(falong, volong,'rows','complete')
    text(min(falong),max(volong)-5,sprintf('Correlation=%.2f (p=%.3f)',rho,pval));
    legend({'b1000','b2000','b3000'})
end
set(gcf,'color','w');
suptitle({'mean fa - mean vol scatterplot','HCP TEST dataset'})
saveas(gcf,fullfile(saveItHere, 'meanFa-volume_scatterplot.png'),'png');
close(gcf)


% REMOVE THE EFFECT OF THE VOLUME IN THE FA AND PLOT AGAIN

mrvNewGraphWin('VOL-FA(norm) scatterplots'); 
tn = 'LeftThalamicRadiation';
tracts      = unique(DV.Struct);
repetitions = height(DV(DV.Struct==tn,:));
FA          = DT(DT.Proj=='HCP' & DT.TRT=='TEST',:);
VOL         = DV(DV.Proj=='HCP' & DV.TRT=='TEST',:);
isequal(unique(FAs.SubjID), unique(VOLs.SubjID))
FAs         = sortrows(FA , 'SubjID');   
VOLs        = sortrows(VOL, 'SubjID');
bvalues     = {1000,2000,3000};
bluecolors{1000}=[0         0    0.5625];
bluecolors{2000}=[0     0.447     0.741];
bluecolors{3000}=[0    0.8125         1];
for nt = 1: length(tracts)
    subplot(4,5,nt)
    tn = tracts(nt);
    % Obtain the values for all b values
    % Calculate the b for all of them, we assume the slope is the same
    % If we want to remove the effect of the b value, we need to calculate all
    % of them
    Falong  = mean(FA{FA.Struct==tn,'Val'},2,'omitnan');
    Volong  = mean(VOL{VOL.Struct==tn,'Val'},2,'omitnan');
    [FalongN,slope] = dr_normRemoveEffect(Falong, Volong);
    [rho,pval] = corr(FalongN, Volong,'rows','complete');
    
    for nb=bvalues    
        % falong = mean(FA{FA.Struct==tn & FA.AcquMD.scanbValue==nb{:},'Val'},2,'omitnan');
        % volong = mean(VOL{VOL.Struct==tn & VOL.AcquMD.scanbValue==nb{:},'Val'},2,'omitnan');
        if nb{:}==1000
            volong  = Volong(1:44);
            falongN = FalongN(1:44);
            scatter(volong,falongN,30,bluecolors{nb{:}},'filled');hold on;
        end
        if nb{:}==2000
            volong  = Volong(44+1:44+44);
            falongN = FalongN(44+1:44+44);
            scatter(volong,falongN,30,bluecolors{nb{:}});hold on;
        end
        if nb{:}==3000
            volong  = Volong(44+44+1:44+44+44);
            falongN = FalongN(44+44+1:44+44+44);
            scatter(volong,falongN,30,bluecolors{nb{:}},'filled');hold on;
        end
        
        xlabel('Volume'); ylabel('FA');
        
        title(sprintf('%s',tn));
    end
    
    text(min(falong),max(volong)-5,sprintf('Correlation=%.2f (p=%.3f)',rho,pval));
    legend({'b1000','b2000','b3000'})
end
set(gcf,'color','w');
suptitle({'normalized mean fa - mean vol scatterplot','HCP TEST dataset'})
saveas(gcf,fullfile(saveItHere, 'nomalized_meanFa-volume_scatterplot.png'),'png');
close(gcf)
%}

%% MEAN FA DISTRIBUTIONS (normality, Cohens'd)
nBoot     = 5000;
normalKsdensity = 'normal';
ylims = [0, 10];
xlims = [0.1, .8];

left_color = [0 0 .55]; right_color = [0.3 0.3 0.3];
% ALL
tmp     = varfun(ostd,dtt,'GroupingVariables',{'Struct','SliceCats'},'InputVariables',{'meanVal'});
replRes = varfun(omean,dtt,'GroupingVariables',{'Struct','SliceCats'},'InputVariables',{'meanVal'});
replRes.Properties.VariableNames{'Fun_meanVal'} = 'Mean';
replRes.SD  = tmp.Fun_meanVal;
replRes     = replRes(:,[1,2,4,5]);
% Unstack it to show it as consecutive tables. Reorder to replicate
replRes     = unstack(replRes,{'Mean','SD'},{'SliceCats'});
replRes     = replRes([3,4,  11,12,   5,6,   7,8,   9,10,   1,2],:);
replRes.Mean_WHLorig = [0.565,0.522,0.517,0.490,0.549,0.534,0.510,0.497,0.489,0.470,0.587,0.575]';
replRes.SD_WHLorig   = [0.028,0.028,0.023,0.03,0.026,0.024,0.028,0.026,0.022,0.024,0.023,0.021]';
replResSorted         = replRes(:,[1, 14,15,  5,11,   6,12,   7,13,   2,8,  3,9,   4,10]);

% Plot barplots with error bars
dt = dtt;
cats      = categories(dt.SliceCats);
catcolors = unique(dt.SliceCatsRGB);
catcolors.cat = categorical(cats);
catcolors     = [table(0,0,0,{'WHLorig'},'VariableNames',{'R','G','B','cat'}); ...
                 catcolors];
catcolors     = catcolors([1,5,6,7,2,3,4],:);

mrvNewGraphWin('Mean FA distributions'); 
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
set(gcf,'color','w');
suptitle({mySupTitle})
saveas(gcf,fullfile(saveItHere, 'FA_Distributions.svg'),'svg');
saveas(gcf,fullfile(saveItHere, 'FA_Distributions.png'),'png');
close(gcf)





% Same plot for TEST-RETEST
HCPunsMeans=unsMeans(unsMeans.Proj=='HCP',:);
HCPunsMeans.SliceCats = removecats(HCPunsMeans.SliceCats);
HCPunsMeansrt=unsMeansrt(unsMeansrt.Proj=='HCP',:);
HCPunsMeansrt.SliceCats = removecats(HCPunsMeansrt.SliceCats);

nBoot     = 5000;
normalKsdensity = 'normal';
mrvNewGraphWin('T-RT Mean FA distributions'); 
for nt = 1: height(replResSorted)
    subplot(3,4,nt)
    tn = char(replResSorted{nt,'Struct'});
    cats      = categories(HCPunsMeans.SliceCats);
    catcolors = unique(HCPunsMeans.SliceCatsRGB);
    catsrt      = categories(HCPunsMeans.SliceCats);
    catcolorsrt = catcolors;
    a = [];
    l = {};
    for sc=1:length(cats)
        x   = HCPunsMeans.(tn)(HCPunsMeans.SliceCats==cats{sc});
        xrt = HCPunsMeansrt.(tn)(HCPunsMeansrt.SliceCats==cats{sc});
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
set(gcf,'color','w');
suptitle({mySupTitle})
saveas(gcf,fullfile(saveItHere, 'FA_Distributions_HCP_TRT.svg'),'svg');
saveas(gcf,fullfile(saveItHere, 'FA_Distributions_HCP_TRT.png'),'png');
close(gcf)

%% WAHL 2010: CORRELATION TABLE FA's REPLICATION (NOT USED)
%{
% Duplicate the correlation table
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
    
WahlorigCorrTable = dr_createCorrelationMatrix(...
                               unsMeans(unsMeans.SliceCats=="WHL1000",:), ... 
                               'tractsOrder',tractsOrder, ...
                               'useBootstrap',true, 'nRep',10000, 'CIrange',95);
WahlorigCorrTable{:,:} = WahlorigCorr;

for ns = 1:3
    % Do Wahl and WH
    if ns ==1
        corrMatTable = dr_createCorrelationMatrix(...
                               unsMeans(unsMeans.SliceCats==string(cats(ns)),:), ... 
                               'tractsOrder',tractsOrder, ...
                               'useBootstrap',true, 'nRep',10000, 'CIrange',95);
        corrMatrices{ns} = corrMatTable;
        display_matrix([{WahlorigCorrTable},corrMatrices(ns)], 'title', strcat(string(cats(ns))," (TEST|RETEST)"), ...
                         'rowheader',WahlTractNames,'colheader',WahlTractNames);
    else
        corrMatTable = dr_createCorrelationMatrix(...
                               unsMeans(unsMeans.SliceCats==string(cats(ns)),:), ... 
                               'tractsOrder',tractsOrder, ...
                               'useBootstrap',true, 'nRep',10000, 'CIrange',95);
        corrMatrices{ns} = corrMatTable;

        display_matrix(corrMatTable{:,:}, 'title', string(cats(ns)), ...
                             'rowheader',WahlTractNames,'colheader',WahlTractNames);
    end
    
    % Do HCP-T
    corrMatTable = dr_createCorrelationMatrix(...
                               unsMeans(unsMeans.SliceCats==string(cats(ns+3)),:), ... 
                               'tractsOrder',tractsOrder, ...
                               'useBootstrap',true, 'nRep',10000, 'CIrange',95);
    corrMatrices{ns+3} = corrMatTable;
    
    % Do HCP-RT
    corrMatTablert = dr_createCorrelationMatrix(...
                               unsMeansrt(unsMeansrt.SliceCats==string(cats(ns+3)),:), ... 
                               'tractsOrder',tractsOrder, ...
                               'useBootstrap',true, 'nRep',10000, 'CIrange',95);
    corrMatrices{ns+6} = corrMatTablert;
    % Plot both together
    display_matrix(corrMatrices([ns+3,ns+6]), 'title', strcat(string(cats(ns+3))," (TEST|RETEST)"), ...
                         'rowheader',WahlTractNames,'colheader',WahlTractNames);
end
                   
% Now plot a summary plot + a ranking
AllCorrMatrices = [{WahlorigCorrTable}, corrMatrices]';
Experiments     = [{'WHLorig'};cellstr(cats(1:3));...
                   strcat(cellstr(cats(4:6)),'TEST');strcat(cellstr(cats(4:6)),'RETEST')];
sumPlotTable    = table(AllCorrMatrices, Experiments);
Unstacked       = {};
for nt=1:height(sumPlotTable)
    sumPlotTable.AllCorrMatrices{nt}.Properties.VariableNames = strrep(WahlTractNames,' ','');
    sumPlotTable.AllCorrMatrices{nt}.Properties.RowNames = strrep(WahlTractNames,' ','');
    sumPlotTable.AllCorrMatrices{nt} = sumPlotTable.AllCorrMatrices{nt}(2:end,1:end-1);
    Unstacked = [Unstacked; {dr_stackTable(sumPlotTable.AllCorrMatrices{nt})}];
end
sumPlotTable.Unstacked = Unstacked;
sumPlotTable = sumPlotTable(:,[2,3]);
% Now make the unstacked col numeric
long = sumPlotTable{1,2}{1};
long.Properties.VariableNames = {'CorName', sumPlotTable{1,1}{1}};
for na = 2:height(sumPlotTable)
    tmp = sumPlotTable{na,2}{1};
    long.(sumPlotTable{na,1}{1}) = tmp.VAL;
end

% Decide: if 1 is cero, the rest are 0, it is like a confidence interval

% RANKING
longNo0     = long(:,[1:6,9, 7,10, 8,11]);
Experiments = longNo0.Properties.VariableNames(2:end)';
% Create index that if any is 0, delete the whole correlation
longNo0(logical(sum(longNo0{:,2:end}==0,2)~=0),:) = [];
ranking = longNo0;
for nr=2:size(ranking,2)
    [~,I] = sort(ranking{:,nr},'Descend');
    [~,r] = sort(I);
    ranking{:,nr} = r;
end
ranking.sum   = sum(ranking{:,2:end},2);
ranking.mean  = mean(longNo0{:,2:end},2);
ranking.sd    = std(longNo0{:,2:end}')';
sortedRanking = sortrows(ranking,'WHLorig');

mrvNewGraphWin('Ranking over Experiments'); 
nrow = 1; ncol = 2;
C    = colormap(jet);
CMAP = C(1:size(C,1)/8:end, :);
colormap(CMAP)
% Plot correlation values
% Select eight to eight, using the ranking in Wahlorig
CorrelationValues = longNo0;
CorrelationValues.RankingWahlOrig = ranking.WHLorig;
CorrelationValues = sortrows(CorrelationValues, 'RankingWahlOrig');
OnlyCorrNames     =  CorrelationValues{:,1};
OnlyCorrValues    =  CorrelationValues{:,2:11};


% HIGHEST RANKS AND CORRELATIONS
selRows = [1:length(OnlyCorrNames)];
subplot(nrow,ncol,1)
for nt=selRows
    if nt <= 8
        plot(sortedRanking{nt,2:length(Experiments)+1}', ...
             'LineStyle','-','LineWidth',3);hold on;
    end
    if nt > 8 && nt <= 16
        plot(sortedRanking{nt,2:length(Experiments)+1}', ...
             'LineStyle','-.','LineWidth',2);hold on;
    end
    if nt > 16 
        plot(sortedRanking{nt,2:length(Experiments)+1}', ...
             'LineStyle',':','LineWidth',2);hold on;
    end
    ylabel('Tract-pair Names');title('Ranking of correlations per Experiment')
    set(gca,'Ydir','reverse','xticklabels',Experiments,'xtick', [1:length(Experiments)],...
            'ytick', [1:height(sortedRanking)], 'XTickLabelRotation', 45, ...
            'XLim', [1,length(Experiments)], 'YAxisLocation', 'left',...
            'yticklabels', strrep(sortedRanking.CorName(selRows),'_','\_'), ...
            'FontName', 'times', 'FontSize', 14,'FontWeight', 'bold')
end

subplot(nrow,ncol,2)
for nt=selRows
    if nt <= 8
        plot(OnlyCorrValues(nt,:)','LineStyle','-','LineWidth',3);hold on;
    end
    if nt > 8 && nt <= 16
        plot(OnlyCorrValues(nt,:)','LineStyle','-.','LineWidth',2);hold on;
    end
    if nt > 16 
        plot(OnlyCorrValues(nt,:)','LineStyle',':','LineWidth',2);hold on;
    end
    ylabel('Correlation'); title('Tract-pair correlation values per experiment')
    set(gca,'xticklabels',Experiments,'xtick', [1:length(Experiments)],...
        'ytick', [0:0.1:1], 'XTickLabelRotation', 45, ...
        'YLim', [0.27,.93], 'XLim', [1,length(Experiments)], ...
        'FontName', 'times', 'FontSize', 14,'FontWeight', 'bold')
end
% legend(strrep(OnlyCorrNames,'_','\_'), 'location', 'SouthWest'), % 'NorthEastOutside')


% Now plot the averages across projects for all correlations
mrvNewGraphWin('AvgOverEXPs'); 
mc       = CorrelationValues;
mc.mean  = mean(mc{:,2:11},2);
mc.sd    = std(mc{:,2:11}')';
mc.upper = mc.mean + mc.sd;
mc.lower = mc.mean - mc.sd;
mc       = sortrows(mc,'mean','descend');
pairNames= strrep(mc.CorName, '_', '\_');
plot(mc.mean,'LineStyle','-' ,'LineWidth',3,'color','k');hold on;
plot(mc.upper ,'LineStyle','-.','LineWidth',2,'color','k');
plot(mc.lower ,'LineStyle','-.','LineWidth',2,'color','k');
ylabel('Correlation'); title('Tract-pair correlations  (averaged experiments, 1 SD)')

set(gca,'xticklabels',pairNames,'xtick', [1:length(pairNames)],...
        'ytick', [0:0.1:1], 'XTickLabelRotation', 45, ...
        'YLim' , [0.3,.9], 'XLim', [1,length(pairNames)], ...
        'FontName', 'times', 'FontSize', 14,'FontWeight', 'bold'); hold off;



% SUMMARY PLOT
tmp      = corrMatTable;
tmp{:,:} = NaN(size(tmp{:,:}));
tmp.Properties.VariableNames = strrep(WahlTractNames,' ','');
tmp.Properties.RowNames     = strrep(WahlTractNames,' ','');
tmp{:,:}(ismissing(tmp)) = 0;
for ii=1:height(mc)
    cornames = split(char(mc.CorName(ii)),'_');
    bat = cornames{1};
    bi  = cornames{2};
    % fprintf('%s_%s :: %f\n',bat, bi, mc.mean(ii))
    tmp{bat,bi} = mc.mean(ii);
end
%}
           
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

% Create the algorithm and plot that compares the means and the CI-s
% PLOT 6 x 11
%{
if plotIt; mrvNewGraphWin('FA Correlation CI Plots'); end
ncol=6;nrow=11;
% Figure paper plots 4 and 9 
for nc=1:length(unique(long.CorName))
    if plotIt; subplot(ncol,nrow,nc); end
    corname  = string(long.CorName(nc));
    % Find the range. 

    
    % DO WHLorig agains the rest
    referenceColumn = 'WHLorig';
    versus          = 'REST';
    takeout = {};
    [WahlInsideCI, isOverlap, minUpper, maxLower, midCI]=dr_compareCI(long, corname,...
                                                    'referenceColumn',referenceColumn,...
                                                    'refLine' , refLine, ...
                                                    'plotIt',plotIt,'showLegend',showLegend, ...
                                                    'takeout',takeout, 'showXnames',showXnames);  



%     % WHLorig vs WHL1000
%     referenceColumn = 'WHLorig';
%     versus          = 'WHL1000';    
%     takeout =  {'YWM1000'     ,'YWM2000'       ,...
%                 'HCP1000TEST','HCP1000RETEST',...
%                 'HCP2000TEST','HCP2000RETEST',...
%                 'HCP3000TEST','HCP3000RETEST'};
%     [WahlInsideCI, isOverlap, minUpper, maxLower, midCI]=dr_compareCI(long, corname,...
%                                                     'referenceColumn',referenceColumn,...
%                                                     'refLine' , refLine, ...
%                                                     'plotIt',plotIt,'showLegend',showLegend, ...
%                                                     'takeout',takeout, 'showXnames',showXnames);  



%     % WHLorig vs b1000s
%     referenceColumn = 'WHLorig';
%     versus          = 'b1000s';    
%     takeout =  {'YWM2000'       , ...
%                 'HCP2000TEST','HCP2000RETEST',...
%                 'HCP3000TEST','HCP3000RETEST'};
%     [WahlInsideCI, isOverlap, minUpper, maxLower, midCI]=dr_compareCI(long, corname,...
%                                                     'referenceColumn',referenceColumn,...
%                                                     'refLine' , refLine, ...
%                                                     'plotIt',plotIt,'showLegend',showLegend, ...
%                                                     'takeout',takeout, 'showXnames',showXnames);  




    % Rotate all the HCPs to obtain them all, do not loop...
%     referenceColumn = 'HCP3000TEST';
%     versus          = 'HCP3000RETEST';    
%     takeout =  {'WHLorig'   ,'WHL1000'     ,...
%                 'YWM1000'     ,'YWM2000'       ,...
%                 'HCP1000TEST','HCP1000RETEST',...
%                 'HCP2000TEST','HCP2000RETEST'};
%     [WahlInsideCI, isOverlap, minUpper, maxLower, midCI]=dr_compareCI(long, corname,...
%                                                     'referenceColumn',referenceColumn,...
%                                                     'refLine' , refLine, ...
%                                                     'plotIt',plotIt,'showLegend',showLegend, ...
%                                                     'takeout',takeout, 'showXnames',showXnames);  
                                                
                                                
    % Now populate the tables
    cornames = split(corname,'_');
    bat = cornames{1};
    bi  = cornames{2};
    % fprintf('%s_%s :: %f\n',bat, bi, mc.mean(ii))
    TWahlInsideCI{bat,bi} = WahlInsideCI;
    TisOverlap{bat,bi}    = isOverlap;
    TminUpper{bat,bi}     = minUpper;
    TmaxLower{bat,bi}     = maxLower;
    TmidCI{bat,bi}        = midCI;
end
if plotIt; 
    set(gcf,'color','w');
    saveas(gcf,fullfile(saveItHere, sprintf('SM2_ALLb1000s_%s_CI%d.png', referenceColumn, CIrange)),'png');
    close(gcf)
end
% Plot the correlation matrices graphically with a colorbar
TReference = TmidCI;
TReference{2:end,1:end-1} = sumPlotTable.AllCorrMatrices{sumPlotTable.Experiments==string(referenceColumn)}{:,:};
allMats = [{TReference{:,:}},{TWahlInsideCI{:,:}},{TisOverlap{:,:}},{TminUpper{:,:}},{TmaxLower{:,:}},{TmidCI{:,:}}];
display_matrix_graphic(allMats, ...
               'title', sprintf('FA Correlation Replication (C.I.: %d%%)',CIrange), ...
               'rowheader',WahlTractNames,'colheader',WahlTractNames, ...
               'asterisk', false);
set(gcf, 'PaperOrientation','portrait')
set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'color','w');
saveas(gcf,fullfile(saveItHere, sprintf('FA_REST_%svs%s_CI%d.png', referenceColumn, versus, CIrange)),'png');
close(gcf)

%}

% PLOT LOWER DIAGONAL (substitutes both figures)
if plotIt
    mrvNewGraphWin('Lower Diag. FA Correlation CI Plots');
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
