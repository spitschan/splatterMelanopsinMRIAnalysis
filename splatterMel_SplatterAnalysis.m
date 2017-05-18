function splatterMel_SplatterAnalysis(ppsRawDataDir, analysisDir)

% Define the paths
outDir = fullfile(analysisDir, 'tables');

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
theObserverAge = [32 28 44 29 32 28 44 29 29 32 28 44 29 32 28 44 29 29 32 28 44 29];
theValidatedObserverAge = [32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32];
theContrastLevels = {[400] [400] [400] [400] ...
    [25 50 100 200 400] [25 50 100 200 400] [25 50 100 200 400] [25 50 100 200 400] [25 50 100 200 400] ...
    [400] [400] [400] [400] ...
    [25 50 100 200 400] [25 50 100 200 400] [25 50 100 200 400] [25 50 100 200 400] [25 50 100 200 400] ...
    [0.25 0.5 1 2] [0.25 0.5 1 2] [0.25 0.5 1 2] [0.25 0.5 1 2]};

fid = fopen(fullfile(outDir, 'TableX_Splatter.csv'), 'w');
fprintf(fid, 'Stimulus,Observer,Actual observer age,Validated observer age,Nominal contrast [%s],Luminance [cd/m2],SD,Irradiance [sc td],[log10 sc td],SD,Irradiance [ph td],[log10 ph td],SD,x chromaticity,SD,y chromaticity,SD,L contrast [%s],SD,M contrast [%s],SD,S contrast [%s],SD,Melanopsin contrast [%s],SD,Rod contrast [%s],SD,LMS contrast [%s],SD,L-M contrast [%s],SD,S-[L+M] contrast [%s],SD,\n', '%', '%', '%', '%', '%', '%', '%', '%', '%');

currDir = pwd;
Mc = [];
for d = 1:length(theDataPaths)
    dataPath = theDataPaths{d};
    
    % Find the folders
    theFolders = dir(fullfile(ppsRawDataDir, dataPath));
    
    % Increment the counter
    clear contrasts;
    clear postRecepContrasts;
    clear luminance;
    clear chromaticity;
    clear irradianceScotTrolands;
    clear irradiancePhotTrolands;
    
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
            
            % Extract the infomation
            observerAgeInYrs = tmp.cals{1}.describe.cache.REFERENCE_OBSERVER_AGE;
            fractionBleached = tmp.cals{1}.describe.cache.data(observerAgeInYrs).describe.fractionBleached;
            pupilDiameterMm = tmp.cals{1}.describe.cache.data(observerAgeInYrs).describe.params.pupilDiameterMm;
            fieldSizeDegrees = tmp.cals{1}.describe.cache.data(observerAgeInYrs).describe.params.fieldSizeDegrees;
            S = tmp.cals{1}.describe.cache.data(observerAgeInYrs).describe.S;
            wls = SToWls(S);
            
            % Calculate luminance and chromaticity
            % Load the CIE functions
            load T_xyz1931
            T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,WlsToS(wls));
            
            % Get the background spectrum
            bgSpd = tmp.cals{1}.modulationBGMeas.meas.pr650.spectrum;
            
            % Calculate luminance and chromaticy
            luminance(f) = T_xyz(2, :)*bgSpd;
            chromaticity(:, f) = (T_xyz([1 2], :)*bgSpd)/sum((T_xyz*bgSpd));
            
            % Calculate irradiance
            pupilAreaMm2 = pi*((pupilDiameterMm/2)^2);
            eyeLengthMm = 17;
            degPerMm = RetinalMMToDegrees(1,eyeLengthMm);
            irradianceWattsPerUm2 = RadianceToRetIrradiance(bgSpd,S,pupilAreaMm2,eyeLengthMm);
            irradianceScotTrolands(f) = RetIrradianceToTrolands(irradianceWattsPerUm2, S, 'Scotopic', [], num2str(eyeLengthMm));
            irradiancePhotTrolands(f) = RetIrradianceToTrolands(irradianceWattsPerUm2, S, 'Photopic', [], num2str(eyeLengthMm));
            
            % Set up the receptor object
            receptorObj = SSTReceptorHuman('obsAgeInYrs', observerAgeInYrs, 'fieldSizeDeg', fieldSizeDegrees, 'obsPupilDiameterMm', pupilDiameterMm);
            T_rec = receptorObj.T.T_energyNormalized;
            T_val = tmp.cals{1}.describe.cache.data(observerAgeInYrs).describe.T_receptors;
            
            % Calculate the numerical difference between the assumed and the
            % reconstructed receptor sensitivities
            for ii = 1:size(T_rec, 1)-1
                tmp1 = sum(T_rec(1, :) - T_val(1, :));
                if tmp1 > 0
                    error('Error: Couldn''t reconstruct receptor sensitivities...');
                end
            end
            
            % Calcualte the contrast
            NContrastLevels = size(tmp.cals{end}.modulationAllMeas, 2)-1;
            for kk = 2:NContrastLevels+1
                modSpd = tmp.cals{1}.modulationAllMeas(1, kk).meas.pr650.spectrum;
                
                % Calculate the nominal contrast
                for jj = 1:size(T_rec, 1)
                    contrasts{kk-1}(:, f) = (T_rec*(modSpd-bgSpd))./(T_rec*bgSpd);
                end
                postRecepContrasts{kk-1}(:, f) = [1 1 1 0 0 ; 1 -1 0 0 0 ; 0 0 1 0 0]' \ contrasts{kk-1}(:, f);
            end
        end
    end
    % Take the average
    lumMean = median(luminance);
    lumSD = std(luminance);
    chromMean = median(chromaticity, 2);
    chromSD = std(chromaticity, [], 2);
    scotTdMean = median(irradianceScotTrolands);
    scotTdSD = std(irradianceScotTrolands);
    photTdMean = median(irradiancePhotTrolands);
    photTdSD = std(irradiancePhotTrolands);
    
    for ii = 1:NContrastLevels
        contrastsMean(:, ii) = median(contrasts{ii}, 2);
        contrastsSD(:, ii) = std(contrasts{ii}, [], 2);
        postRecepContrastsMean(:, ii) = median(postRecepContrasts{ii}, 2);
        postRecepContrastsSD(:, ii) = std(postRecepContrasts{ii}, [], 2);
    end
    
    % Assemble the data
    Mb = [];
    for ii = 1:NContrastLevels
        M = [];
        M = [M lumMean lumSD scotTdMean log10(scotTdMean) scotTdSD photTdMean log10(photTdMean) photTdSD chromMean(1) chromSD(1) chromMean(2) chromSD(2)];
        for m = 1:size(contrastsMean, 1)
            M = [M 100*contrastsMean(m, ii) 100*contrastsSD(m, ii)];
        end
        for m = 1:size(postRecepContrastsMean, 1)
            M = [M 100*postRecepContrastsMean(m, ii) 100*postRecepContrastsSD(m, ii)];
        end
        Mb = [Mb ; M];
    end
    
    % Write out the data
    for ii = 1:NContrastLevels
        fprintf(fid, '%s,%s,%i,%i,%.2f,', theStimuli{d}, theObservers{d}, theObserverAge(d), theValidatedObserverAge(d), theContrastLevels{d}(ii));
        for jj = 1:size(Mb, 2)
            if ii > 1 && jj < 14
                fprintf(fid, ',');
            else
                fprintf(fid, '%.2f,', Mb(ii, jj));
            end
        end
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n');
    fprintf(fid, '\n');
end
cd(currDir);
fclose(fid);