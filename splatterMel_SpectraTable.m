function splatterMel_SpectraTable(ppsRawDataDir, analysisDir)

fprintf('> Running %s\n', mfilename);

% Define the paths
outDir = fullfile(analysisDir, 'tables', 'spectra');
if ~isdir(outDir)
    mkdir(outDir);
end

theDataPaths = {'MaxLMS400pct/HERO_asb1/040716/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
    'MaxLMS400pct/HERO_aso1/033016/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
    'MaxLMS400pct/HERO_gka1/040116/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
    'MaxLMS400pct/HERO_mxs1/040816/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
    'MaxLMSCRF/HERO_asb1/060816/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MaxLMSCRF/HERO_aso1/060116/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MaxLMSCRF/HERO_gka1/060616/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MaxLMSCRF/HERO_mxs1/061016/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MaxLMSCRF/HERO_mxs1/062816/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableBStubby1_ND03CassetteB/20-Jun-2016_15_17_02/validation' ...
    'MaxMel400pct/HERO_asb1/032416/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
    'MaxMel400pct/HERO_aso1/032516/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
    'MaxMel400pct/HERO_gka1/033116/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
    'MaxMel400pct/HERO_mxs1/040616/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
    'MaxMelCRF/HERO_asb1/060716/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MaxMelCRF/HERO_aso1/053116/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MaxMelCRF/HERO_gka1/060216/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MaxMelCRF/HERO_mxs1/060916/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MaxMelCRF/HERO_mxs1/061016/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'SplatterControlCRF/HERO_asb1/051016/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation' ...
    'SplatterControlCRF/HERO_aso1/042916/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation' ...
    'SplatterControlCRF/HERO_gka1/050616/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation' ...
    'SplatterControlCRF/HERO_mxs1/050916/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation'};

theStimuli = {'LMS 400%' 'LMS 400%'  'LMS 400%'  'LMS 400%'  ...
    'LMS CRF'  'LMS CRF'  'LMS CRF'  'LMS CRF' 'LMS CRF'...
    'Mel 400%' 'Mel 400%'  'Mel 400%'  'Mel 400%'  ...
    'Mel CRF'  'Mel CRF'  'Mel CRF'  'Mel CRF' 'Mel CRF' ...
    'Splatter CRF' 'Splatter CRF' 'Splatter CRF'  'Splatter CRF'};
theObservers = {'ASB' 'ASO' 'GKA' 'MXS' ...
    'ASB' 'ASO' 'GKA' 'MXS [1]' 'MXS [2]'...
    'ASB' 'ASO' 'GKA' 'MXS' ...
    'ASB' 'ASO' 'GKA' 'MXS [1]' 'MXS [2]' ...
    'ASB' 'ASO' 'GKA' 'MXS'};
theContrastLevels = {[400] [400] [400] [400] ...
    [25 50 100 200 400] [25 50 100 200 400] [25 50 100 200 400] [25 50 100 200 400] [25 50 100 200 400] ...
    [400] [400] [400] [400] ...
    [25 50 100 200 400] [25 50 100 200 400] [25 50 100 200 400] [25 50 100 200 400] [25 50 100 200 400] ...
    [0.25 0.5 1 2] [0.25 0.5 1 2] [0.25 0.5 1 2] [0.25 0.5 1 2]};
theContrastScalars = {[1] [1] [1] [1] ...
    [1/16 1/8 1/4 1/2 1] [1/16 1/8 1/4 1/2 1] [1/16 1/8 1/4 1/2 1] [1/16 1/8 1/4 1/2 1] [1/16 1/8 1/4 1/2 1] ...
    [1] [1] [1] [1] ...
    [1/16 1/8 1/4 1/2 1] [1/16 1/8 1/4 1/2 1] [1/16 1/8 1/4 1/2 1] [1/16 1/8 1/4 1/2 1] [1/16 1/8 1/4 1/2 1] ...
    [0.25 0.5 1 2] [0.25 0.5 1 2] [0.25 0.5 1 2] [0.25 0.5 1 2]};
theObserverAge = [32 28 44 29 32 28 44 29 29 32 28 44 29 32 28 44 29 29 32 28 44 29];
theValidatedObserverAge = [32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32];

wls = SToWls([380 2 201]);

currDir = pwd;
Mc = [];
for d = 1:length(theDataPaths)
    outFile = fullfile(outDir, [theStimuli{d} '_' theObservers{d} '.csv']);
    fid = fopen(outFile, 'w');
    dataPath = theDataPaths{d};
    
    % Find the folders
    theFolders = dir(fullfile(ppsRawDataDir, dataPath));
    
    % Increment the counter
    clear bgSpdVal;
    clear modSpdVal;
    clear modSpdValMean;
    clear modSpdValSD;
    
    for k = length(theFolders):-1:1
        % remove non-folders
        if ~theFolders(k).isdir
            theFolders(k) = [ ];
            continue;
        end
        
        % remove folders starting with .
        fname = theFolders(k).name;
        if fname(1) == '.'
            theFolders(k) = [ ];
        end
    end
    
    % Iterate over the folders
    for f = 1:length(theFolders)
        % Go to the folder
        if isdir(fullfile(ppsRawDataDir, dataPath, theFolders(f).name))
            cd(fullfile(ppsRawDataDir, dataPath, theFolders(f).name));
        end
        
        % Find the only MAT file there is going to be
        theMATFile = dir([pwd '/*.mat']);
        
        if ~isempty(theMATFile)
            % Load the MAT file
            tmp = load(theMATFile.name);
            
            % Get the background spectrum
            bgSpdVal(:, f) = tmp.cals{1}.modulationBGMeas.meas.pr650.spectrum;
            bgSpdNom = tmp.cals{1}.describe.cache.data(theObserverAge(d)).backgroundSpd;
            
            % Calculate the contrast
            NContrastLevels = size(tmp.cals{end}.modulationAllMeas, 2)-1;
            for kk = 2:NContrastLevels+1
                modSpdVal{kk-1}(:, f) = tmp.cals{1}.modulationAllMeas(1, kk).meas.pr650.spectrum;
                modSpdNom(:, kk-1) = tmp.cals{1}.describe.cache.data(theObserverAge(d)).backgroundSpd+theContrastScalars{d}(kk-1)*tmp.cals{1}.describe.cache.data(theObserverAge(d)).differenceSpd;
            end
        end
        
        bgSpdValMean = mean(bgSpdVal, 2);
        bgSpdValSD = std(bgSpdVal, [], 2);
        
        for kk = 1:NContrastLevels
            modSpdValMean(:, kk) = mean(modSpdVal{kk}, 2);
            modSpdValSD(:, kk) = std(modSpdVal{kk}, [], 2);
        end
    end
    
    % Gather all the spectra
    theHeader = [theObservers{d} ', ' theStimuli{d}];
    theString = sprintf('Wavelength [nm],Background nominal [%i yrs],', theObserverAge(d));
    Mc = [wls];
    % First the nominal ones
    Mc = [Mc bgSpdNom];
    for ii = 1:NContrastLevels
        if strcmp(theStimuli{d}, 'Splatter CRF')
            theString = [theString sprintf('Modulation nominal [x%.2f %i yrs],', ...
                theContrastLevels{d}(ii), theObserverAge(d))];
        else
            theString = [theString sprintf('Modulation nominal [%i%s %i yrs],', ...
                theContrastLevels{d}(ii), '%', theObserverAge(d))];
        end
        Mc = [Mc modSpdNom(:, ii)];
    end
    
    % Then the validated ones
    theString = [theString sprintf('Background validated [%i yrs],', theValidatedObserverAge(d))];
    Mc = [Mc bgSpdValMean];
    for ii = 1:NContrastLevels
        if strcmp(theStimuli{d}, 'Splatter CRF')
            theString = [theString sprintf('Modulation validated [x%.2f %i yrs]', ...
                theContrastLevels{d}(ii), theValidatedObserverAge(d))];
        else
            theString = [theString sprintf('Modulation validated [%i%s %i yrs]', ...
                theContrastLevels{d}(ii), '%', theValidatedObserverAge(d))];
        end
        if ii < NContrastLevels
            theString = [theString ','];
        end
        Mc = [Mc modSpdValMean(:, ii)];
    end
    fprintf(fid, '%s', theString);
    fprintf(fid, '\n');
    fclose(fid);
    dlmwrite(outFile, Mc, '-append');
end
cd(currDir);