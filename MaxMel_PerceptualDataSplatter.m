% Define the paths
dropboxBasePath = '/Users/mspits/Dropbox (Aguirre-Brainard Lab)/MELA_materials/MaxMelPulsePsychophysics';
outDir = fullfile(pwd, 'figures');
outTableDir = fullfile(pwd, 'tables', 'perceptualdata');
if ~isdir(outDir)
    mkdir(outDir);
end

if ~isdir(outTableDir)
    mkdir(outTableDir);
end

outFileSplatter = fullfile(pwd, 'tables', 'TableX_SplatterPerceptualData.csv');

% Define the # of subjects
NSubjects = 20;

subjectIDs={'MELA_0074',...
    'MELA_0087',...
    'MELA_0089',...
    'MELA_0026',...
    'MELA_0082',...
    'MELA_0038',...
    'MELA_0094',...
    'MELA_0096',...
    'MELA_0088',...
    'MELA_0079',...
    'MELA_0043',...
    'MELA_0073',...
    'MELA_0080',...
    'MELA_0090',...
    'MELA_0037',...
    'MELA_0049',...
    'MELA_0050',...
    'MELA_0075',...
    'MELA_0077',...
    'MELA_0081',...
    };

theMelData = {'040517/Cache-MelanopsinDirectedSuperMaxMel_MELA_0074_040517' ...
    '040417/Cache-MelanopsinDirectedSuperMaxMel_MELA_0087_040417' ...
    '033117/Cache-MelanopsinDirectedSuperMaxMel_MELA_0089_033117' ...
    '033117/Cache-MelanopsinDirectedSuperMaxMel_MELA_0026_033117' ...
    '032917/Cache-MelanopsinDirectedSuperMaxMel_MELA_0082_032917' ...
    '032917/Cache-MelanopsinDirectedSuperMaxMel_MELA_0038_032917' ...
    '021417/Cache-MelanopsinDirectedSuperMaxMel_MELA_0094_021417' ...
    '032817/Cache-MelanopsinDirectedSuperMaxMel_MELA_0096_032817' ...
    '032117/Cache-MelanopsinDirectedSuperMaxMel_MELA_0088_032117' ...
    '032017/Cache-MelanopsinDirectedSuperMaxMel_MELA_0079_032017' ...
    '022717/Cache-MelanopsinDirectedSuperMaxMel_MELA_0043_022717' ...
    '022717/Cache-MelanopsinDirectedSuperMaxMel_MELA_0073_022717' ...
    '022117/Cache-MelanopsinDirectedSuperMaxMel_MELA_0080_022117' ...
    '022117/Cache-MelanopsinDirectedSuperMaxMel_MELA_0090_022117' ...
    '021717/Cache-MelanopsinDirectedSuperMaxMel_MELA_0037_021717' ...
    '020317/Cache-MelanopsinDirectedSuperMaxMel_MELA_0049_020317' ...
    '020717/Cache-MelanopsinDirectedSuperMaxMel_MELA_0050_020717' ...
    '020617/Cache-MelanopsinDirectedSuperMaxMel_MELA_0075_020617' ...
    '020817/Cache-MelanopsinDirectedSuperMaxMel_MELA_0077_020817' ...
    '021017/Cache-MelanopsinDirectedSuperMaxMel_MELA_0081_021017'};
    %' ...
    
theLMSData = {'040517/Cache-LMSDirectedSuperMaxLMS_MELA_0074_040517' ...
    '040417/Cache-LMSDirectedSuperMaxLMS_MELA_0087_040417' ...
    '033117/Cache-LMSDirectedSuperMaxLMS_MELA_0089_033117' ...
    '033117/Cache-LMSDirectedSuperMaxLMS_MELA_0026_033117' ...
    '032917/Cache-LMSDirectedSuperMaxLMS_MELA_0082_032917' ...
    '032917/Cache-LMSDirectedSuperMaxLMS_MELA_0038_032917' ...
    '021417/Cache-LMSDirectedSuperMaxLMS_MELA_0094_021417' ...
    '032817/Cache-LMSDirectedSuperMaxLMS_MELA_0096_032817' ...
    '032117/Cache-LMSDirectedSuperMaxLMS_MELA_0088_032117' ...
    '032017/Cache-LMSDirectedSuperMaxLMS_MELA_0079_032017' ...
    '022717/Cache-LMSDirectedSuperMaxLMS_MELA_0043_022717' ...
    '022717/Cache-LMSDirectedSuperMaxLMS_MELA_0073_022717' ...
    '022117/Cache-LMSDirectedSuperMaxLMS_MELA_0080_022117' ...
    '022117/Cache-LMSDirectedSuperMaxLMS_MELA_0090_022117' ...
    '021717/Cache-LMSDirectedSuperMaxLMS_MELA_0037_021717' ...
    '020317/Cache-LMSDirectedSuperMaxLMS_MELA_0049_020317' ...
    '020717/Cache-LMSDirectedSuperMaxLMS_MELA_0050_020717' ...
    '020617/Cache-LMSDirectedSuperMaxLMS_MELA_0075_020617' ...
    '020817/Cache-LMSDirectedSuperMaxLMS_MELA_0077_020817' ...
    '021017/Cache-LMSDirectedSuperMaxLMS_MELA_0081_021017'};
    %' ...

% Turn off some warnings
warning('off', 'MATLAB:load:cannotInstantiateLoadedVariable');
warning('off', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
warning('off', 'MATLAB:class:EnumerableClassNotFound');

% Set up the wl vector
wls = SToWls([380 2 201]);

% Calculate luminance and chromaticity
% Load the CIE functions
load T_xyz1931
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,WlsToS(wls));

% Set up the header
headerSpectra = 'Wavelength [nm],Background (Mel),Modulation (Mel),Background (LMS),Modulation (LMS)\n';

% Save splatter table
fid = fopen(outFileSplatter, 'w');
fprintf(fid, 'Stimulus,Observer,Observer age,Nominal contrast [%s],Luminance [cd/m2],Irradiance [sc td],[log10 sc td],Irradiance [ph td],[log10 ph td],x chromaticity,y chromaticity,L contrast [%s],M contrast [%s],S contrast [%s],Melanopsin contrast [%s],Rod contrast [%s],LMS contrast [%s],L-M contrast [%s]\n', '%', '%', '%', '%', '%', '%', '%', '%');
fclose(fid);

% Load the files
for d = 1:NSubjects
    outFile = fullfile(outTableDir, ['Spectra_sub' num2str(d, '%03g') '.csv']);
    fid = fopen(outFile, 'w');
    fprintf(fid, headerSpectra);
    fclose(fid);
    
    M = [wls];
    % Clear data
    clear bgSpdVal;
    clear modSpdVal;
    
    % Load the mel data first
    dataPath = theMelData{d};
    
    % Find the folders
    theFolders = dir(fullfile(dropboxBasePath, dataPath));
    
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
            modSpdVal(:, f) = tmp.cals{1}.modulationMaxMeas.meas.pr650.spectrum;
            
            % Make the receptor object
            observerAgeInYrs = tmp.cals{1}.describe.cache.OBSERVER_AGE;
            fractionBleached = tmp.cals{1}.describe.cache.data(observerAgeInYrs).describe.fractionBleached;
            pupilDiameterMm = tmp.cals{1}.describe.cache.data(observerAgeInYrs).describe.params.pupilDiameterMm;
            fieldSizeDegrees = 64;
            receptorObj = SSTReceptorHuman('verbosity', 'none', 'obsAgeInYrs', observerAgeInYrs, 'fieldSizeDeg', fieldSizeDegrees, 'obsPupilDiameterMm', pupilDiameterMm);
        end
    end
    bgSpdValMean = median(bgSpdVal, 2);
    modSpdValMean = median(modSpdVal, 2);
    
    T_receptors = receptorObj.T.T_energyNormalized;
    for jj = 1:5
        contrastsFixed(jj) = (T_receptors(jj, :)*(modSpdValMean(:, end)-bgSpdValMean))./(T_receptors(jj, :)*bgSpdValMean);
    end
    postRecepContrastsFixedMel(:, d) = [1 1 1 0 0; 1 -1 0 0 0; 0 0 1 0 0]' \ contrastsFixed';
    
    M = [M bgSpdValMean modSpdValMean];
    
    % Calculate luminance and chromaticy
    luminance = T_xyz(2, :)*bgSpdValMean;
    chromaticity = (T_xyz([1 2], :)*bgSpdValMean)/sum((T_xyz*bgSpdValMean));
    
    % Calculate irradiance
    pupilAreaMm2 = pi*((pupilDiameterMm/2)^2);
    eyeLengthMm = 17;
    degPerMm = RetinalMMToDegrees(1,eyeLengthMm);
    irradianceWattsPerUm2 = RadianceToRetIrradiance(bgSpdValMean, WlsToS(wls),pupilAreaMm2,eyeLengthMm);
    irradianceScotTrolands = RetIrradianceToTrolands(irradianceWattsPerUm2, WlsToS(wls), 'Scotopic', [], num2str(eyeLengthMm));
    irradiancePhotTrolands = RetIrradianceToTrolands(irradianceWattsPerUm2, WlsToS(wls), 'Photopic', [], num2str(eyeLengthMm));

    fid = fopen(outFileSplatter, 'a');
    fprintf(fid, '%s,%s,%i,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n', 'Mel', ['sub' num2str(d, '%03g')], observerAgeInYrs, 400, luminance, irradianceScotTrolands, log10(irradianceScotTrolands), ...
        irradiancePhotTrolands, log10(irradiancePhotTrolands), chromaticity(1), chromaticity(2), 100*contrastsFixed(1), 100*contrastsFixed(2), 100*contrastsFixed(3), 100*contrastsFixed(4), 100*contrastsFixed(5), 100*postRecepContrastsFixedMel(1, d), 100*postRecepContrastsFixedMel(2, d));
    fclose(fid);
    
    % Clear data
    clear bgSpdVal;
    clear modSpdVal;
    
    % Load the LMS data first
    dataPath = theLMSData{d};
    
    % Find the folders
    theFolders = dir(fullfile(dropboxBasePath, dataPath));
       
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
     test(d) = length(theFolders);
     
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
            modSpdVal(:, f) = tmp.cals{1}.modulationMaxMeas.meas.pr650.spectrum;
            
            % Make the receptor object
            observerAgeInYrs = tmp.cals{1}.describe.cache.OBSERVER_AGE;
            fractionBleached = tmp.cals{1}.describe.cache.data(observerAgeInYrs).describe.fractionBleached;
            pupilDiameterMm = tmp.cals{1}.describe.cache.data(observerAgeInYrs).describe.params.pupilDiameterMm;
            fieldSizeDegrees = 64;
            receptorObj = SSTReceptorHuman('obsAgeInYrs', observerAgeInYrs, 'fieldSizeDeg', fieldSizeDegrees, 'obsPupilDiameterMm', pupilDiameterMm);
            
        end
    end
    bgSpdValMean = median(bgSpdVal, 2);
    modSpdValMean = median(modSpdVal, 2);
    
    M = [M bgSpdValMean modSpdValMean];
    
    T_receptors = receptorObj.T.T_energyNormalized;
    for jj = 1:5
        contrastsFixed(jj) = (T_receptors(jj, :)*(modSpdValMean(:, end)-bgSpdValMean))./(T_receptors(jj, :)*bgSpdValMean);
    end
    postRecepContrastsFixedLMS(:, d) = [1 1 1 0 0; 1 -1 0 0 0; 0 0 1 0 0]' \ contrastsFixed';
    
        % Calculate luminance and chromaticy
    luminance = T_xyz(2, :)*bgSpdValMean;
    chromaticity = (T_xyz([1 2], :)*bgSpdValMean)/sum((T_xyz*bgSpdValMean));
    
    % Calculate irradiance
    pupilAreaMm2 = pi*((pupilDiameterMm/2)^2);
    eyeLengthMm = 17;
    degPerMm = RetinalMMToDegrees(1,eyeLengthMm);
    irradianceWattsPerUm2 = RadianceToRetIrradiance(bgSpdValMean, WlsToS(wls),pupilAreaMm2,eyeLengthMm);
    irradianceScotTrolands = RetIrradianceToTrolands(irradianceWattsPerUm2, WlsToS(wls), 'Scotopic', [], num2str(eyeLengthMm));
    irradiancePhotTrolands = RetIrradianceToTrolands(irradianceWattsPerUm2, WlsToS(wls), 'Photopic', [], num2str(eyeLengthMm));

    fid = fopen(outFileSplatter, 'a');
    fprintf(fid, '%s,%s,%i,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n', 'LMS', ['sub' num2str(d, '%03g')], observerAgeInYrs, 400, luminance, irradianceScotTrolands, log10(irradianceScotTrolands), ...
        irradiancePhotTrolands, log10(irradiancePhotTrolands), chromaticity(1), chromaticity(2), 100*contrastsFixed(1), 100*contrastsFixed(2), 100*contrastsFixed(3), 100*contrastsFixed(4), 100*contrastsFixed(5), 100*postRecepContrastsFixedMel(1, d), 100*postRecepContrastsFixedMel(2, d));
    fclose(fid);
    
    % Write the table
    dlmwrite(outFile, M, '-append');
end

fig1 = figure;
% LMS spectra
subplot(1, 2, 1);
ScatterplotWithHistogram(postRecepContrastsFixedLMS(1, :), postRecepContrastsFixedLMS(2, :), ...
    'XBinWidth', 0.005, 'YBinWidth', 0.015, 'XLim', [3.9 4.1], 'YLim', [-0.05 0.05], ...
    'XLabel', 'L+M+S contrast', 'YLabel', 'L-M  contrast', 'Color', [1 0.5 0 ; 1 0.5 0], ...
    'XRefLines', [4 4 ; -0.05 0.05], 'YRefLines', [3.9 4.1 ; 0 0], ...
    'MaxP', 1, 'PlotMarginals', false);
% Get the error ellipse
plot(mean(postRecepContrastsFixedLMS(1, :)), mean(postRecepContrastsFixedLMS(2, :)), '+r');
[X, Y] = get_error_ellipse([postRecepContrastsFixedLMS(1, :) ; postRecepContrastsFixedLMS(2, :)]', 0.95); hold on;
plot(X, Y, '-k');

subplot(1, 2, 2);
ScatterplotWithHistogram(postRecepContrastsFixedLMS(1, :), postRecepContrastsFixedLMS(3, :), ...
    'XBinWidth', 0.005, 'YBinWidth', 0.015, 'XLim', [3.9 4.1], 'YLim', [-0.3 0.3], ...
    'XLabel', 'L+M+S contrast', 'YLabel', 'S-(L+M+S)  contrast', 'Color', [0 0 1 ; 0 0 1], ...
     'XRefLines', [4 4 ; -0.3 0.3], 'YRefLines', [3.9 4.1 ; 0 0], ...
    'MaxP', 1, 'PlotMarginals', false);
% Get the error ellipse
plot(mean(postRecepContrastsFixedLMS(1, :)), mean(postRecepContrastsFixedLMS(3, :)), '+r');
[X, Y] = get_error_ellipse([postRecepContrastsFixedLMS(1, :) ; postRecepContrastsFixedLMS(3, :)]', 0.95); hold on;
plot(X, Y, '-k');

% Save figure
set(fig1, 'PaperPosition', [0 0 8 3]);
set(fig1, 'PaperSize', [8 3]);
set(fig1, 'Color', 'w');
set(fig1, 'InvertHardcopy', 'off');
saveas(fig1, fullfile(outDir, 'FigureX_PerceptualDataValidationLMS'), 'pdf');
close(fig1);
%%
fig2 = figure;
% Mel spectra
subplot(1, 2, 1);
ScatterplotWithHistogram(postRecepContrastsFixedMel(1, :), postRecepContrastsFixedMel(2, :), ...
    'XBinWidth', 0.005, 'YBinWidth', 0.015, 'XLim', [-0.1 0.1], 'YLim', [-0.05 0.05], ...
    'XLabel', 'L+M+S contrast', 'YLabel', 'L-M  contrast', 'Color', [1 0.5 0 ; 1 0.5 0], ...
    'XRefLines', [0 0 ; -0.1 0.1], 'YRefLines', [-0.1 0.1 ; 0 0], ...
    'MaxP', 1, 'PlotMarginals', false);
% Get the error ellipse
plot(mean(postRecepContrastsFixedMel(1, :)), mean(postRecepContrastsFixedMel(2, :)), '+r');
[X, Y] = get_error_ellipse([postRecepContrastsFixedMel(1, :) ; postRecepContrastsFixedMel(2, :)]', 0.95); hold on;
plot(X, Y, '-k');

subplot(1, 2, 2);
ScatterplotWithHistogram(postRecepContrastsFixedMel(1, :), postRecepContrastsFixedMel(3, :), ...
    'XBinWidth', 0.005, 'YBinWidth', 0.015, 'XLim', [-0.1 0.1], 'YLim', [-0.3 0.3], ...
    'XLabel', 'L+M+S contrast', 'YLabel', 'S-(L+M+S)  contrast', 'Color', [0 0 1 ; 0 0 1], ...
    'XRefLines', [0 0 ; -0.3 0.3], 'YRefLines', [-0.1 0.1 ; 0 0], ...
    'MaxP', 1, 'PlotMarginals', false);
% Get the error ellipse
plot(mean(postRecepContrastsFixedMel(1, :)), mean(postRecepContrastsFixedMel(3, :)), '+r');
[X, Y] = get_error_ellipse([postRecepContrastsFixedMel(1, :) ; postRecepContrastsFixedMel(3, :)]', 0.95); hold on;
plot(X, Y, '-k');

% Save figure
set(fig2, 'PaperPosition', [0 0 8 3]);
set(fig2, 'PaperSize', [8 3]);
set(fig2, 'Color', 'w');
set(fig2, 'InvertHardcopy', 'off');
saveas(fig2, fullfile(outDir, 'FigureX_PerceptualDataValidationMel'), 'pdf');
close(fig2);

warning('on', 'MATLAB:load:cannotInstantiateLoadedVariable');
warning('on', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
warning('on', 'MATLAB:class:EnumerableClassNotFound');

%% See if there's a correlation between color rating and L-M contrast
for ii = 1:NSubjects
   colorRatingLMS(ii) = table2array(foldedDataTable{ii}(12, 4));
end