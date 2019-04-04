function [structPairs, newStructs, LIndex, RIndex] = dr_obtainPairs(T, longWide)





switch longWide
    case {'long'}
        structPairs = categories(T.Struct);
        LIndex      = contains(structPairs,{'Left'});
        RIndex      = contains(structPairs,{'Right'});
        structPairs = structPairs(LIndex | RIndex); 
        newStructs  = structPairs(find(contains(structPairs,'Left')));
        newStructs  = strrep(newStructs, 'Left', '');
        LIndex      = NaN;
        RIndex      = NaN;
        
    case {'wide'}
        structPairs = T.Properties.VariableNames;
        Llogic      = contains(structPairs,{'Left'});
        Rlogic      = contains(structPairs,{'Right'});
        structPairs = structPairs(Llogic | Rlogic); 

        % Check even number of structures
        if mod(length(structPairs),2); 
          error('There are not a even number of structures')
        end

        % Obtain left and right indexes
        LIndex = find(contains(structPairs,'Left'));
        RIndex = find(contains(structPairs,'Right'));
        %Check they are in the same order
        sStructure = strrep(structPairs, 'Left', '');
        sStructure = strrep(sStructure, 'Right', '');
        if ~isequal(sStructure(LIndex), sStructure(RIndex))
            error('Left-Right structures dont match')
        end
        newStructs = sStructure(LIndex);
    otherwise
        error(sprintf('The format %s is not known', longWide))
end        
   





end



% afqPairs = { 'LeftArcuate'             , 'RightArcuate'            , ...
%              'LeftCingulumCingulate'  , 'RightCingulumCingulate'  , ...
%              'LeftCingulumHippocampus', 'RightCingulumHippocampus', ...
%              'LeftCorticospinal'      , 'RightCorticospinal'      , ...
%              'LeftIFOF'               , 'RightIFOF'               , ...
%              'LeftILF'                , 'RightILF'                , ...
%              'LeftSLF'                , 'RightSLF'                , ...
%              'LeftThalamicRadiation'  , 'RightThalamicRadiation'  , ...
%              'LeftUncinate'           , 'RightUncinate'           };

