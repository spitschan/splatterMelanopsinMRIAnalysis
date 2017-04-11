% Define the paths
dropboxBasePath = '/Users/mspits/Dropbox (Aguirre-Brainard Lab)/MELA_materials/MaxMelPulsePsychophysics';
outDir = fullfile(pwd, 'figures');
if ~isdir(outDir)
    mkdir(outDir);
end

% Define the # of subjects
NSubjects = 20;

subjectIDs={'MELA_0090',...
    'MELA_0074',...
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
    'MELA_0037',...
    'MELA_0049',...
    'MELA_0050',...
    'MELA_0075',...
    'MELA_0077',...
    'MELA_0081',...
    };

theMelData = {'022117/Cache-MelanopsinDirectedSuperMaxMel_MELA_0090_022117' ...
    '040517/Cache-MelanopsinDirectedSuperMaxMel_MELA_0074_040517' ...
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
    '021717/Cache-MelanopsinDirectedSuperMaxMel_MELA_0037_021717' ...
    '020317/Cache-MelanopsinDirectedSuperMaxMel_MELA_0049_020317' ...
    '020717/Cache-MelanopsinDirectedSuperMaxMel_MELA_0050_020717' ...
    '020617/Cache-MelanopsinDirectedSuperMaxMel_MELA_0075_020617' ...
    '020817/Cache-MelanopsinDirectedSuperMaxMel_MELA_0077_020817' ...
    '021017/Cache-MelanopsinDirectedSuperMaxMel_MELA_0081_021017'};
    %' ...
    
theLMSData = {'022117/Cache-LMSDirectedSuperMaxLMS_MELA_0090_022117' ...
    '040517/Cache-LMSDirectedSuperMaxLMS_MELA_0074_040517' ...
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
    '021717/Cache-LMSDirectedSuperMaxLMS_MELA_0037_021717' ...
    '020317/Cache-LMSDirectedSuperMaxLMS_MELA_0049_020317' ...
    '020717/Cache-LMSDirectedSuperMaxLMS_MELA_0050_020717' ...
    '020617/Cache-LMSDirectedSuperMaxLMS_MELA_0075_020617' ...
    '020817/Cache-LMSDirectedSuperMaxLMS_MELA_0077_020817' ...
    '021017/Cache-LMSDirectedSuperMaxLMS_MELA_0081_021017'};
    %' ...


warning('off', 'MATLAB:load:cannotInstantiateLoadedVariable');
warning('off', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
warning('off', 'MATLAB:class:EnumerableClassNotFound');
% Load the files
for d = 1:NSubjects
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
            receptorObj = SSTReceptorHuman('obsAgeInYrs', observerAgeInYrs, 'fieldSizeDeg', fieldSizeDegrees, 'obsPupilDiameterMm', pupilDiameterMm);
            
        end
    end
    bgSpdValMean = median(bgSpdVal, 2);
    modSpdValMean = median(modSpdVal, 2);
    
    T_receptors = receptorObj.T.T_energyNormalized;
    for jj = 1:3
        contrastsFixed(jj) = (T_receptors(jj, :)*(modSpdValMean(:, end)-bgSpdValMean))./(T_receptors(jj, :)*bgSpdValMean);
    end
    postRecepContrastsFixedMel(:, d) = [1 1 1 ; 1 -1 0 ; 0 0 1]' \ contrastsFixed';
    
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
    
    T_receptors = receptorObj.T.T_energyNormalized;
    for jj = 1:3
        contrastsFixed(jj) = (T_receptors(jj, :)*(modSpdValMean(:, end)-bgSpdValMean))./(T_receptors(jj, :)*bgSpdValMean);
    end
    postRecepContrastsFixedLMS(:, d) = [1 1 1 ; 1 -1 0 ; 0 0 1]' \ contrastsFixed';
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

fig2 = figure;
% Mel spectra
subplot(1, 2, 1);
ScatterplotWithHistogram(postRecepContrastsFixedMel(1, :), postRecepContrastsFixedMel(2, :), ...
    'XBinWidth', 0.005, 'YBinWidth', 0.015, 'XLim', [-0.05 0.05], 'YLim', [-0.05 0.05], ...
    'XLabel', 'L+M+S contrast', 'YLabel', 'L-M  contrast', 'Color', [1 0.5 0 ; 1 0.5 0], ...
    'XRefLines', [0 0 ; -0.05 0.05], 'YRefLines', [-0.05 0.05 ; 0 0], ...
    'MaxP', 1, 'PlotMarginals', false);
% Get the error ellipse
plot(mean(postRecepContrastsFixedMel(1, :)), mean(postRecepContrastsFixedMel(2, :)), '+r');
[X, Y] = get_error_ellipse([postRecepContrastsFixedMel(1, :) ; postRecepContrastsFixedMel(2, :)]', 0.95); hold on;
plot(X, Y, '-k');

subplot(1, 2, 2);
ScatterplotWithHistogram(postRecepContrastsFixedMel(1, :), postRecepContrastsFixedMel(3, :), ...
    'XBinWidth', 0.005, 'YBinWidth', 0.015, 'XLim', [-0.05 0.05], 'YLim', [-0.3 0.3], ...
    'XLabel', 'L+M+S contrast', 'YLabel', 'S-(L+M+S)  contrast', 'Color', [0 0 1 ; 0 0 1], ...
    'XRefLines', [0 0 ; -0.3 0.3], 'YRefLines', [-0.05 0.05 ; 0 0], ...
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