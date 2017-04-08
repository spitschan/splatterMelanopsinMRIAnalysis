% Define the paths
dropboxBasePath = '/Users/spitschan/Dropbox (Aguirre-Brainard Lab)';
outDir = fullfile(pwd, 'figures');
if ~isdir(outDir)
    mkdir(outDir);
end

whichDataSet = 'CRF';
% Define the data paths
switch whichDataSet
    case 'CRF'
        theDataPaths = {'MELA_data/MelanopsinMR_fMRI/MaxMelCRF/HERO_asb1/060716/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
            'MELA_data/MelanopsinMR_fMRI/MaxMelCRF/HERO_aso1/053116/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
            'MELA_data/MelanopsinMR_fMRI/MaxMelCRF/HERO_gka1/060216/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
            'MELA_data/MelanopsinMR_fMRI/MaxMelCRF/HERO_mxs1/060916/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
            'MELA_data/MelanopsinMR_fMRI/SplatterControlCRF/HERO_asb1/051016/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation' ...
            'MELA_data/MelanopsinMR_fMRI/SplatterControlCRF/HERO_aso1/042916/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation' ...
            'MELA_data/MelanopsinMR_fMRI/SplatterControlCRF/HERO_gka1/050616/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation' ...
            'MELA_data/MelanopsinMR_fMRI/SplatterControlCRF/HERO_mxs1/050916/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation'};
    case '400Only'
        theDataPaths = {'MELA_data/MelanopsinMR_fMRI/MaxMel400pct/HERO_asb1/032416/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
            'MELA_data/MelanopsinMR_fMRI/MaxMel400pct/HERO_aso1/032516/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
            'MELA_data/MelanopsinMR_fMRI/MaxMel400pct/HERO_gka1/033116/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
            'MELA_data/MelanopsinMR_fMRI/MaxMel400pct/HERO_mxs1/040616/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
            'MELA_data/MelanopsinMR_fMRI/SplatterControlCRF/HERO_asb1/051016/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation' ...
            'MELA_data/MelanopsinMR_fMRI/SplatterControlCRF/HERO_aso1/042916/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation' ...
            'MELA_data/MelanopsinMR_fMRI/SplatterControlCRF/HERO_gka1/050616/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation' ...
            'MELA_data/MelanopsinMR_fMRI/SplatterControlCRF/HERO_mxs1/050916/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation'};
end

% Define some meta data manually

wls = SToWls([380 2 201]);

% Make some empty variables
MelPostRecepContrasts = [];
SplatterPostRecepContrasts = [];

currDir = pwd;
indDiffParams = [];
adjIndDiffParams = [];
for d = 1:length(theDataPaths)
    outFile1 = fullfile(outDir, ['FigureX_PhysiologicalSplatter_' whichDataSet '.pdf']);
    outFile2 = fullfile(outDir, ['FigureX_PhysiologicalSplatterParam_' whichDataSet '.pdf']);
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
    for jj = 1:size(receptorObj.Ts{1}.T_energyNormalized, 1)
        contrastsFixed(jj) = (T_receptors(jj, :)*(modSpdValMean(:, end)-bgSpdValMean))./(T_receptors(jj, :)*bgSpdValMean);
    end
    postRecepContrastsFixed(:, d) = [1 1 1 ; 1 -1 0 ; 0 0 1]' \ contrastsFixed';
    
    for ii = 1:NSamples
        T_receptors = receptorObj.Ts{ii}.T_energyNormalized;
        for jj = 1:size(receptorObj.Ts{ii}.T_energyNormalized, 1)
            contrastsStochastic(jj, ii) = (T_receptors(jj, :)*(modSpdValMean(:, end)-bgSpdValMean))./(T_receptors(jj, :)*bgSpdValMean);
        end
        postRecepContrastsStochastic(:, ii) = [1 1 1 ; 1 -1 0 ; 0 0 1]' \ contrastsStochastic(:, ii);
        indDiffParams = [indDiffParams receptorObj.Ts{ii}.indDiffParams];
        adjIndDiffParams = [adjIndDiffParams receptorObj.Ts{ii}.adjIndDiffParams];
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
    %     fig1 = figure;
    %     XAxLims = [-0.1 0.1]; YAxLims = [-0.3 0.3];
    %     XNoeeminalContrast = 0; YNominalContrast = 0;
    %     edi(postRecepContrastsStochastic(1, :), postRecepContrastsStochastic(3, :), ...
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
    %     set(fig1, 'PaperPosition', [0 0 3 3]);
    %     set(fig1, 'PaperSize', [3 3]);
    %     set(fig1, 'Color', 'w');
    %     set(fig1, 'InvertHardcopy', 'off');
    %     saveas(fig1, outFile2, 'png');
    %     close(fig1); close(fig1);
    
    if d < 5
        MelPostRecepContrasts =  [MelPostRecepContrasts postRecepContrastsStochastic];
    else
        SplatterPostRecepContrasts = [SplatterPostRecepContrasts postRecepContrastsStochastic];
    end
    
    
end

fig1 = figure;
subplot(1, 2, 1);
ScatterplotWithHistogram(MelPostRecepContrasts(1, :), MelPostRecepContrasts(2, :), ...
    'XBinWidth', 0.005, 'YBinWidth', 0.005, 'XLim', [-0.1 0.1], 'YLim', [-0.1 0.1], ...
    'XLabel', '|L+M+S contrast|', 'YLabel', '|L-M contrast|', 'Color', [1 0 0 ; 1 0 0], ...
    'MaxP', 1, 'PlotMarginals', false);
ScatterplotWithHistogram(SplatterPostRecepContrasts(1, :), SplatterPostRecepContrasts(2, :), ...
    'XBinWidth', 0.005, 'YBinWidth', 0.005, 'XLim', [-0.1 0.1], 'YLim', [-0.1 0.1], ...
    'XLabel', '|L+M+S contrast|', 'YLabel', '|L-M contrast|', 'Color', [1 0.5 0 ; 1 0.5 0], ...
    'MaxP', 1, 'PlotMarginals', false, ...
    'XRefLines', [0 0; -0.1 0.1], 'YRefLines', [-0.1 0.1; 0 0]);
% Add contrast at target
plot(postRecepContrastsFixed(1, 1), postRecepContrastsFixed(2, 1), '+g', 'MarkerSize', 10);
plot(postRecepContrastsFixed(1, 2), postRecepContrastsFixed(2, 2), 'og', 'MarkerSize', 10);
plot(postRecepContrastsFixed(1, 3), postRecepContrastsFixed(2, 3), 'xg', 'MarkerSize', 10);
plot(postRecepContrastsFixed(1, 4), postRecepContrastsFixed(2, 4), '*g', 'MarkerSize', 10);

plot(postRecepContrastsFixed(1, 5), postRecepContrastsFixed(2, 5), '+k', 'MarkerSize', 10);
plot(postRecepContrastsFixed(1, 6), postRecepContrastsFixed(2, 6), 'ok', 'MarkerSize', 10);
plot(postRecepContrastsFixed(1, 7), postRecepContrastsFixed(2, 7), 'xk', 'MarkerSize', 10);
plot(postRecepContrastsFixed(1, 8), postRecepContrastsFixed(2, 8), '*k', 'MarkerSize', 10);

plot_error_ellipse([MelPostRecepContrasts(1, :) ; MelPostRecepContrasts(2, :)]')


% S splatter
subplot(1, 2, 2);
ScatterplotWithHistogram(MelPostRecepContrasts(1, :), MelPostRecepContrasts(3, :), ...
    'XBinWidth', 0.005, 'YBinWidth', 0.015, 'XLim', [-0.1 0.1], 'YLim', [-0.3 0.3], ...
    'XLabel', '|L+M+S contrast|', 'YLabel', '|S-(L+M+S)  contrast|', 'Color', [0 0 1 ; 0 0 1], ...
    'MaxP', 1, 'PlotMarginals', false);
ScatterplotWithHistogram(SplatterPostRecepContrasts(1, :), SplatterPostRecepContrasts(3, :), ...
    'XBinWidth', 0.005, 'YBinWidth', 0.015, 'XLim', [-0.1 0.1], 'YLim', [-0.3 0.3], ...
    'XLabel', '|L+M+S contrast|', 'YLabel', '|S-(L+M+S) contrast|', 'Color', [1 0.5 0 ; 1 0.5 0], ...
    'MaxP', 1, 'PlotMarginals', false, ...
    'XRefLines', [0 0; -0.3 0.3], 'YRefLines', [-0.1 0.1; 0 0]);
% Add contrast at target
plot(postRecepContrastsFixed(1, 1), postRecepContrastsFixed(3, 1), '+g', 'MarkerSize', 10);
plot(postRecepContrastsFixed(1, 2), postRecepContrastsFixed(3, 2), 'og', 'MarkerSize', 10);
plot(postRecepContrastsFixed(1, 3), postRecepContrastsFixed(3, 3), 'xg', 'MarkerSize', 10);
plot(postRecepContrastsFixed(1, 4), postRecepContrastsFixed(3, 4), '*g', 'MarkerSize', 10);

plot(postRecepContrastsFixed(1, 5), postRecepContrastsFixed(3, 5), '+k', 'MarkerSize', 10);
plot(postRecepContrastsFixed(1, 6), postRecepContrastsFixed(3, 6), 'ok', 'MarkerSize', 10);
plot(postRecepContrastsFixed(1, 7), postRecepContrastsFixed(3, 7), 'xk', 'MarkerSize', 10);
plot(postRecepContrastsFixed(1, 8), postRecepContrastsFixed(3, 8), '*k', 'MarkerSize', 10);

plot_error_ellipse([MelPostRecepContrasts(1, :) ; MelPostRecepContrasts(3, :)]')

% Save figure
set(fig1, 'PaperPosition', [0 0 6 3]);
set(fig1, 'PaperSize', [6 3]);
set(fig1, 'Color', 'w');
set(fig1, 'InvertHardcopy', 'off');
saveas(fig1, outFile1, 'pdf');
close(fig1);

%% "Parametric" analysis
fig2 = figure;
% Lens
subplot(3, 8, 1);
plot([indDiffParams(1:4000).dlens], MelPostRecepContrasts(1, :), '.k'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);
plot([0 0], [-0.1 0.1], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast');
title('Lens');

subplot(3, 8, 9);
plot([indDiffParams(1:4000).dlens], MelPostRecepContrasts(2, :), '.r'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);
plot([0 0], [-0.1 0.1], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast');

subplot(3, 8, 17);
plot([indDiffParams(1:4000).dlens], MelPostRecepContrasts(3, :), '.b'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.35 0.35]);
plot([0 0], [-0.35 0.35], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast'); xlabel('Deviation [%]');

% Macula
subplot(3, 8, 2);
plot([indDiffParams(1:4000).dmac], MelPostRecepContrasts(1, :), '.k'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);
plot([0 0], [-0.1 0.1], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast');
title('Macula');

subplot(3, 8, 10);
plot([indDiffParams(1:4000).dmac], MelPostRecepContrasts(2, :), '.r'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);
plot([0 0], [-0.1 0.1], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast');

subplot(3, 8, 18);
plot([indDiffParams(1:4000).dmac], MelPostRecepContrasts(3, :), '.b'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.35 0.35]);
plot([0 0], [-0.35 0.35], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast'); xlabel('Deviation [%]');

% Photopigment density
for ii = 1:4000
    LDensity(ii) = indDiffParams(ii).dphotopigment(1);
    MDensity(ii) = indDiffParams(ii).dphotopigment(2);
    SDensity(ii) = indDiffParams(ii).dphotopigment(3);
end
% L cone
subplot(3, 8, 3);
plot(LDensity, MelPostRecepContrasts(1, :), '.k'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);
plot([0 0], [-0.1 0.1], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast');
title('L cone density');

subplot(3, 8, 11);
plot(LDensity, MelPostRecepContrasts(2, :), '.r'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);
plot([0 0], [-0.1 0.1], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast');

subplot(3, 8, 19);
plot(LDensity, MelPostRecepContrasts(3, :), '.b'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.35 0.35]);
plot([0 0], [-0.35 0.35], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast'); xlabel('Deviation [%]');

% M cone
subplot(3, 8, 4);
plot(MDensity, MelPostRecepContrasts(1, :), '.k'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);
plot([0 0], [-0.1 0.1], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast');
title('M cone density');

subplot(3, 8, 12);
plot(MDensity, MelPostRecepContrasts(2, :), '.r'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);
plot([0 0], [-0.1 0.1], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast');

subplot(3, 8, 20);
plot(MDensity, MelPostRecepContrasts(3, :), '.b'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.35 0.35]);
plot([0 0], [-0.35 0.35], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast'); xlabel('Deviation [%]');

% S cone
subplot(3, 8, 5);
plot(SDensity, MelPostRecepContrasts(1, :), '.k'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);
plot([0 0], [-0.1 0.1], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast');
title('S cone density');

subplot(3, 8, 13);
plot(SDensity, MelPostRecepContrasts(2, :), '.r'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);
plot([0 0], [-0.1 0.1], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast');

subplot(3, 8, 21);
plot(SDensity, MelPostRecepContrasts(3, :), '.b'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.35 0.35]);
plot([0 0], [-0.35 0.35], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast'); xlabel('Deviation [%]');

% Lambda-max shifts
for ii = 1:4000
    LShift(ii) = indDiffParams(ii).lambdaMaxShift(1);
    MShift(ii) = indDiffParams(ii).lambdaMaxShift(2);
    SShift(ii) = indDiffParams(ii).lambdaMaxShift(3);
end
% L cone
subplot(3, 8, 6);
plot(LShift, MelPostRecepContrasts(1, :), '.k'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);
xlim([-10 10]);
plot([0 0], [-0.1 0.1], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast');
title('L cone shift');

subplot(3, 8, 14);
plot(LShift, MelPostRecepContrasts(2, :), '.r'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);
xlim([-10 10]);
plot([0 0], [-0.1 0.1], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast');

subplot(3, 8, 22);
plot(LShift, MelPostRecepContrasts(3, :), '.b'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.35 0.35]);
xlim([-10 10]);
plot([0 0], [-0.35 0.35], '-k');

% M cone
subplot(3, 8, 7);
plot(MShift, MelPostRecepContrasts(1, :), '.k'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);
xlim([-10 10]);
plot([0 0], [-0.1 0.1], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast');
title('M cone shift');

subplot(3, 8, 15);
plot(MShift, MelPostRecepContrasts(2, :), '.r'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);
xlim([-10 10]);
plot([0 0], [-0.1 0.1], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast');

subplot(3, 8, 23);
plot(MShift, MelPostRecepContrasts(3, :), '.b'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.35 0.35]);
xlim([-10 10]);
plot([0 0], [-0.35 0.35], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast'); xlabel('Deviation [nm]');

% S cone
subplot(3, 8, 8);
plot(SShift, MelPostRecepContrasts(1, :), '.k'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);
xlim([-10 10]);
plot([0 0], [-0.1 0.1], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast');
title('S cone shift');

subplot(3, 8, 16);
plot(SShift, MelPostRecepContrasts(2, :), '.r'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);
xlim([-10 10]);
plot([0 0], [-0.1 0.1], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast');

subplot(3, 8, 24);
plot(SShift, MelPostRecepContrasts(3, :), '.b'); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.35 0.35]);
xlim([-10 10]);
plot([0 0], [-0.35 0.35], '-k');
plot([-100 100], [0 0], '-k');
ylabel('Contrast'); xlabel('Deviation [nm]');

%% Save figure
set(fig2, 'PaperPosition', [0 0 15 6]);
set(fig2, 'PaperSize', [15 6]);
set(fig2, 'Color', 'w');
set(fig2, 'InvertHardcopy', 'off');
saveas(fig2, outFile2, 'pdf');

%% Do linear regression
for ii = 1:3
    x = [[indDiffParams(1:4000).dlens]' [indDiffParams(1:4000).dmac]' LDensity' MDensity' SDensity' LShift' MShift' SShift'];
    y = MelPostRecepContrasts(ii, :);
    fitlm(x, y')
end