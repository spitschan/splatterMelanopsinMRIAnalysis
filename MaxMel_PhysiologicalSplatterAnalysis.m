% Define the paths
dropboxBasePath = '/Users/spitschan/Dropbox (Aguirre-Brainard Lab)';
outDir = fullfile(pwd, 'figures');
if ~isdir(outDir)
    mkdir(outDir);
end

theDataPaths = {'MELA_data/MelanopsinMR_fMRI/MaxMelCRF/HERO_asb1/060716/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxMelCRF/HERO_aso1/053116/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxMelCRF/HERO_gka1/060216/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxMelCRF/HERO_mxs1/060916/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MELA_data/MelanopsinMR_fMRI/SplatterControlCRF/HERO_asb1/051016/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation' ...
    'MELA_data/MelanopsinMR_fMRI/SplatterControlCRF/HERO_aso1/042916/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation' ...
    'MELA_data/MelanopsinMR_fMRI/SplatterControlCRF/HERO_gka1/050616/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation' ...
    'MELA_data/MelanopsinMR_fMRI/SplatterControlCRF/HERO_mxs1/050916/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation'};


theStimuli = {'Mel 400%' 'Mel 400%'  'Mel 400%'  'Mel 400%'  ...
    'Mel CRF'  'Mel CRF'  'Mel CRF'  'Mel CRF' 'Mel CRF' ...
    'Splatter CRF' 'Splatter CRF' 'Splatter CRF'  'Splatter CRF'};
theObservers = {'ASB' 'ASO' 'GKA' 'MXS' ...
    'ASB' 'ASO' 'GKA' 'MXS [1]' 'MXS [2]' ...
    'ASB' 'ASO' 'GKA' 'MXS'};
theContrastLevels = {[400] [400] [400] [400] ...
    [25 50 100 200 400] [25 50 100 200 400] [25 50 100 200 400] [25 50 100 200 400] [25 50 100 200 400] ...
    [0.25 0.5 1 2] [0.25 0.5 1 2] [0.25 0.5 1 2] [0.25 0.5 1 2]};
wls = SToWls([380 2 201]);

currDir = pwd;
Mc = [];
for d = 1:length(theDataPaths)
    outFile1 = fullfile(outDir, [theStimuli{d} '_' theObservers{d} '_LM.png']);
    outFile2 = fullfile(outDir, [theStimuli{d} '_' theObservers{d} '_S.png']);
    dataPath = theDataPaths{d};
    
    % Find the folders
    theFolders = dir(fullfile(dropboxBasePath, dataPath));
    
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
        if isdir(fullfile(dropboxBasePath, dataPath, theFolders(f).name))
            cd(fullfile(dropboxBasePath, dataPath, theFolders(f).name));
        end
        
        % Find the only MAT file there is going to be
        theMATFile = dir([pwd '/*.mat']);
        
        if ~isempty(theMATFile)
            % Load the MAT file
            tmp = load(theMATFile.name);
            
            % Get the background spectrum
            bgSpdVal(:, f) = tmp.cals{1}.modulationBGMeas.meas.pr650.spectrum;
            bgSpdNom = tmp.cals{1}.modulationBGMeas.predictedSpd;
            
            % Calculate the contrast
            NContrastLevels = size(tmp.cals{end}.modulationAllMeas, 2)-1;
            for kk = 2:NContrastLevels+1
                modSpdVal{kk-1}(:, f) = tmp.cals{1}.modulationAllMeas(1, kk).meas.pr650.spectrum;
            end
        end
    end
    bgSpdValMean = mean(bgSpdVal, 2);
    bgSpdValSD = std(bgSpdVal, [], 2);
    
    for kk = 1:NContrastLevels
        modSpdValMean(:, kk) = mean(modSpdVal{kk}, 2);
        modSpdValSD(:, kk) = std(modSpdVal{kk}, [], 2);
    end
    
    
    %% Make the receptor object
    observerAgeInYrs = tmp.cals{1}.describe.cache.REFERENCE_OBSERVER_AGE;
    fractionBleached = tmp.cals{1}.describe.cache.data(observerAgeInYrs).describe.fractionBleached;
    pupilDiameterMm = tmp.cals{1}.describe.cache.data(observerAgeInYrs).describe.params.pupilDiameterMm;
    fieldSizeDegrees = tmp.cals{1}.describe.cache.data(observerAgeInYrs).describe.params.fieldSizeDegrees;
    receptorObj = SSTReceptorHuman('obsAgeYrs', observerAgeInYrs, 'fieldSizeDeg', fieldSizeDegrees, 'obsPupilDiameterMm', pupilDiameterMm);
    NSamples = 1000;
    receptorObj.makeSpectralSensitivitiesStochastic('NSamples', NSamples);
    
    
    T_receptors = receptorObj.T.T_energyNormalized;
    for jj = 1:size(receptorObj.Ts{ii}.T_energyNormalized, 1)
        contrastsFixed(jj) = (T_receptors(jj, :)*(modSpdValMean(:, end)-bgSpdValMean))./(T_receptors(jj, :)*bgSpdValMean);
    end
    postRecepContrastsFixed = [1 1 0 ; 1 -1 0 ; 0 0 1]' \ contrastsFixed';
    
    for ii = 1:NSamples
        T_receptors = receptorObj.Ts{ii}.T_energyNormalized;
        for jj = 1:size(receptorObj.Ts{ii}.T_energyNormalized, 1)
            contrastsStochastic(jj, ii) = (T_receptors(jj, :)*(modSpdValMean(:, end)-bgSpdValMean))./(T_receptors(jj, :)*bgSpdValMean);
        end
        postRecepContrastsStochastic(:, ii) = [1 1 1 ; 1 -1 0 ; 0 0 1]' \ contrastsStochastic(:, ii);
    end
    
    %     theRGB = DefaultReceptorColors;
    %     fig1 = figure;
    %     XAxLims = [-0.1 0.1]; YAxLims = [-0.1 0.1];
    %     XNominalContrast = 0; YNominalContrast = 0;
    %     ScatterPlotNoHistogram(postRecepContrastsStochastic(1, :), postRecepContrastsStochastic(2, :), ...
    %         'XLim', XAxLims, 'YLim', YAxLims, 'XBinWidth', 0.01, 'YBinWidth', 0.01, ...
    %         'XLabel', 'L+M+S contrast', 'YLabel', 'L-M contrast', ...
    %         'XRefLines', [XAxLims ; YNominalContrast YNominalContrast], ...
    %         'YRefLines', [XNominalContrast XNominalContrast ; YAxLims], ...
    %         'XNominalContrast', XNominalContrast, ...
    %         'YNominalContrast', YNominalContrast, ...
    %         'Color', [theRGB(1, :) ; theRGB(2, :)]);
    %     plot(postRecepContrastsFixed(1), postRecepContrastsFixed(2), 'xg', 'MarkerSize', 10);
    %     title([theStimuli{d} '-' theObservers{d}]);
    %
    %     % Save the figure
    %     set(fig1, 'PaperPosition', [0 0 3 3]);
    %     set(fig1, 'PaperSize', [3 3]);
    %     set(fig1, 'Color', 'w');
    %     set(fig1, 'InvertHardcopy', 'off');
    %     cd(currDir);
    %     saveas(fig1, outFile1, 'png');
    %
    %     fig2 = figure;
    %     XAxLims = [-0.1 0.1]; YAxLims = [-0.3 0.3];
    %     XNominalContrast = 0; YNominalContrast = 0;
    %     ScatterPlotNoHistogram(postRecepContrastsStochastic(1, :), postRecepContrastsStochastic(3, :), ...
    %         'XLim', XAxLims, 'YLim', YAxLims, 'XBinWidth', 0.01, 'YBinWidth', 0.03, ...
    %         'XLabel', 'L+M+S contrast', 'YLabel', 'S contrast', ...
    %         'XRefLines', [XAxLims ; YNominalContrast YNominalContrast], ...
    %         'YRefLines', [XNominalContrast XNominalContrast ; YAxLims], ...
    %         'XNominalContrast', XNominalContrast, ...
    %         'YNominalContrast', YNominalContrast, ...
    %         'Color', [theRGB(1, :) ; theRGB(3, :)]);
    %     plot(postRecepContrastsFixed(1), postRecepContrastsFixed(3), 'xg', 'MarkerSize', 10);
    %     title([theStimuli{d} '-' theObservers{d}]);
    %
    %     % Save the figure
    %     set(fig2, 'PaperPosition', [0 0 3 3]);
    %     set(fig2, 'PaperSize', [3 3]);
    %     set(fig2, 'Color', 'w');
    %     set(fig2, 'InvertHardcopy', 'off');
    %     saveas(fig2, outFile2, 'png');
    %     close(fig1); close(fig2);
    
    % Absolute
    
    subplot(1, 2, 1);
    if d < 5
        plot(abs(postRecepContrastsStochastic(1, :)), abs(postRecepContrastsStochastic(2, :)), '.k'); hold on;
    else
        plot(abs(postRecepContrastsStochastic(1, :)), abs(postRecepContrastsStochastic(2, :)), '.r'); hold on;
    end
    xlim([-0.0025 0.1]); xlabel('L+M+S contrast');
    ylim([-0.0025 0.1]); ylabel('L-M contrast');
    pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off;
    
    subplot(1, 2, 2);
    if d < 5
        plot(abs(postRecepContrastsStochastic(1, :)), abs(postRecepContrastsStochastic(3, :)), '.k'); hold on;
    else
        plot(abs(postRecepContrastsStochastic(1, :)), abs(postRecepContrastsStochastic(3, :)), '.r'); hold on;
    end
    xlim([-0.0025 0.1]); xlabel('L+M+S contrast');
    ylim([-0.0075 0.3]); ylabel('S contrast');
    pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off;
    ls
end
        set(gcf, 'PaperPosition', [0 0 6 3]);
        set(gcf, 'PaperSize', [6 3]);
        set(gcf, 'Color', 'w');
        set(gcf, 'InvertHardcopy', 'off');
        saveas(gcf, '~/Desktop/Test.png', 'png');

ls