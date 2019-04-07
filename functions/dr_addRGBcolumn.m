function T = dr_addRGBcolumn(T)




% Create color map by group
% Projects will be asigned main colors
% Slice categories under the groups, with shades of the color
% Projects are sorted by alphabet
% BLACK
mainColor{1} = [0,0,0];
% BLUE
mainColor{100} = [0  , 76 , 153]/255;
mainColor{101} = [0  , 102, 204]/255;
mainColor{102} = [0  , 128, 255]/255;
mainColor{103} = [51 , 153, 255]/255;
mainColor{104} = [102, 178, 255]/255;
mainColor{105} = [153, 204, 255]/255;
mainColor{106} = [0, 0.4470, 0.7410];
mainColor{107} = [0, 0.4470, 0.7410];
mainColor{108} = [0, 0.4470, 0.7410];
mainColor{109} = [0, 0.4470, 0.7410];
mainColor{110} = [0, 0.4470, 0.7410];
mainColor{111} = [0, 0.4470, 0.7410];
mainColor{112} = [0, 0.4470, 0.7410];
mainColor{113} = [0, 0.4470, 0.7410];
mainColor{114} = [0, 0.4470, 0.7410];
mainColor{115} = [0, 0.4470, 0.7410];
mainColor{116} = [0, 0.4470, 0.7410];
mainColor{117} = [0, 0.4470, 0.7410];

% GREEN
mainColor{200} = [50,205,50]/255;
mainColor{201} = [0,100,0]/255;
mainColor{202} = [152,251,152]/255;
mainColor{203} = [0.9290, 0.6940, 0.1250];
mainColor{204} = [0.9290, 0.6940, 0.1250];
mainColor{205} = [0.9290, 0.6940, 0.1250];
mainColor{206} = [0.9290, 0.6940, 0.1250];
mainColor{207} = [0.9290, 0.6940, 0.1250];
mainColor{208} = [0.9290, 0.6940, 0.1250];
mainColor{209} = [0.9290, 0.6940, 0.1250];
mainColor{210} = [0.9290, 0.6940, 0.1250];

% RED
mainColor{300} = [0.7410, 0, 0];
mainColor{301} = [1,   0.6471,  0];
mainColor{302} = [1.0000    0.1875         0];
mainColor{303} = [1.0000    0.5000         0];
mainColor{304} = [0.8500, 0.3250, 0.0980];
mainColor{305} = [255,99,71];
mainColor{306} = [255,127,80];
mainColor{307} = [0.8500, 0.3250, 0.0980];
mainColor{308} = [0.8500, 0.3250, 0.0980];
mainColor{309} = [0.8500, 0.3250, 0.0980];
mainColor{310} = [0.8500, 0.3250, 0.0980];


mainColor{400} = [0.4940, 0.1840, 0.5560];
mainColor{401} = [0.4940, 0.1840, 0.5560];
mainColor{402} = [0.4940, 0.1840, 0.5560];
mainColor{403} = [0.4940, 0.1840, 0.5560];
mainColor{404} = [0.4940, 0.1840, 0.5560];
mainColor{405} = [0.4940, 0.1840, 0.5560];
mainColor{406} = [0.4940, 0.1840, 0.5560];
mainColor{407} = [0.4940, 0.1840, 0.5560];
mainColor{408} = [0.4940, 0.1840, 0.5560];
mainColor{409} = [0.4940, 0.1840, 0.5560];
mainColor{410} = [0.4940, 0.1840, 0.5560];

mainColor{500} = [0.4660, 0.6740, 0.1880];
mainColor{501} = [0.4660, 0.6740, 0.1880];
mainColor{502} = [0.4660, 0.6740, 0.1880];
mainColor{503} = [0.4660, 0.6740, 0.1880];
mainColor{504} = [0.4660, 0.6740, 0.1880];
mainColor{505} = [0.4660, 0.6740, 0.1880];
mainColor{506} = [0.4660, 0.6740, 0.1880];
mainColor{507} = [0.4660, 0.6740, 0.1880];
mainColor{508} = [0.4660, 0.6740, 0.1880];
mainColor{509} = [0.4660, 0.6740, 0.1880];
mainColor{510} = [0.4660, 0.6740, 0.1880];

mainColor{600} = [0.3010, 0.7450, 0.9330];
mainColor{601} = [0.3010, 0.7450, 0.9330];
mainColor{602} = [0.3010, 0.7450, 0.9330];
mainColor{603} = [0.3010, 0.7450, 0.9330];
mainColor{604} = [0.3010, 0.7450, 0.9330];
mainColor{605} = [0.3010, 0.7450, 0.9330];
mainColor{606} = [0.3010, 0.7450, 0.9330];
mainColor{607} = [0.3010, 0.7450, 0.9330];
mainColor{608} = [0.3010, 0.7450, 0.9330];
mainColor{609} = [0.3010, 0.7450, 0.9330];
mainColor{610} = [0.3010, 0.7450, 0.9330];

mainColor{700} = [0.6350, 0.0780, 0.1840];
mainColor{701} = [0.6350, 0.0780, 0.1840];
mainColor{702} = [0.6350, 0.0780, 0.1840];
mainColor{703} = [0.6350, 0.0780, 0.1840];
mainColor{704} = [0.6350, 0.0780, 0.1840];
mainColor{705} = [0.6350, 0.0780, 0.1840];
mainColor{706} = [0.6350, 0.0780, 0.1840];
mainColor{707} = [0.6350, 0.0780, 0.1840];
mainColor{708} = [0.6350, 0.0780, 0.1840];
mainColor{709} = [0.6350, 0.0780, 0.1840];
mainColor{710} = [0.6350, 0.0780, 0.1840];

% Create a new variable that will be a RGB colour to be applied to 
% each sliced categories
T.SliceCatsRGB = array2table(repmat(NaN(height(T),1),[1,3]));
T.SliceCatsRGB.Properties.VariableNames = {'R', 'G', 'B'};
maincats = sort(categories(T.Proj));
spcats   = categories(T.SliceCats);


% Create the index per every category
colorIndT     = cell2table(spcats);
colorIndT.RGB = 999*ones(height(colorIndT),1);
if length(maincats)==1
    colorNumber = 100;
    for ns = 1:length(colorIndT.RGB)
        % colorIndT.RGB(colorIndT.spcats==string(spMatches.spcats{ns})) = colorNumber;
        colorIndT.RGB(colorIndT.spcats==string(colorIndT.spcats{ns})) = colorNumber;
        colorNumber = colorNumber+1;
    end    
else
    for nm = 1:length(maincats)
        colorNumber = nm*100;
        spMatches   = colorIndT(contains(colorIndT.spcats,maincats{nm}),:);
        for ns = 1:height(spMatches)
            colorIndT.RGB(colorIndT.spcats==string(spMatches.spcats{ns})) = colorNumber;
            colorNumber = colorNumber+1;
        end
    end
end
    
% Now asign all values
for nc=1:length(spcats)
        colorIndex = colorIndT.RGB(colorIndT.spcats==string(spcats{nc}));
        howMany = height(T{T.SliceCats==spcats{nc}, 'SliceCatsRGB'});
        T{T.SliceCats==spcats{nc}, 'SliceCatsRGB'} = ...
                                       array2table(repmat(mainColor{colorIndex},[howMany,1]));
end
    
    
    
    
    
end