% Define the paths
dropboxBasePath = '/Users/spitschan/Dropbox (Aguirre-Brainard Lab)';
outDir = fullfile(pwd, 'figures');
if ~isdir(outDir)
    mkdir(outDir);
end

theDataPaths = {'MELA_data/MelanopsinMR_fMRI/MaxLMS400pct/HERO_asb1/040716/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxLMSCRF/HERO_asb1/060816/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxMel400pct/HERO_asb1/032416/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxMelCRF/HERO_asb1/060716/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MELA_data/MelanopsinMR_fMRI/SplatterControlCRF/HERO_asb1/051016/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation'};

theStimuli = {'LMS 400%' 'LMS CRF' 'Mel 400%' 'Mel CRF' 'Splatter CRF'};
close all;

outFig1 = fullfile(outDir, 'FigureX_Spectra.pdf');
theFig = figure;

wls = SToWls([380 2 201]);
load T_xyz1931
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,WlsToS(wls));

currDir = pwd;
for d = 1:length(theDataPaths)
    subplot(3, 2, d);
    dataPath = theDataPaths{d};
    
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
    cd(fullfile(dropboxBasePath, dataPath, theFolders(1).name));
    
    % Find the only MAT file there is going to be
    theMATFile = dir([pwd '/*.mat']);
    
    if ~isempty(theMATFile)
        % Load the MAT file
        tmp = load(theMATFile.name);
        
        % Get the background spectrum
        bgSpd = tmp.cals{1}.modulationBGMeas.predictedSpd; hold on;
        
        % Extract the chromaticity
        
        
        % Get the background spectrum
        bgSpd = tmp.cals{1}.modulationBGMeas.meas.pr650.spectrum; hold on;
        chromaticity(:, d) = (T_xyz([1 2], :)*bgSpd)/sum((T_xyz*bgSpd));
        
        % Set up the color map
        NContrastLevels = size(tmp.cals{1}.modulationAllMeas, 2)-1;
        if NContrastLevels == 1
            theRGB = copper(NContrastLevels+1);
        else
            theRGB = copper(NContrastLevels+1);
        end
        for kk = 2:NContrastLevels+1
            % Get the spectrum
            modSpd = tmp.cals{1}.modulationAllMeas(1, kk).predictedSpd;
            
            % ... and plot it.
            plot(wls, modSpd, 'Color', theRGB(kk, :)); hold on;
        end
        % Plot the background spectrum
        plot(wls, bgSpd, '-k'); hold on;
    end
    
    % Set some plot properties
    xlabel('Wavelength [nm]');
    ylabel('Radiance [W/m2/sr/nm]');
    xlim([380 780]);
    if d < 5
        ylim([0 0.03]);
    end
    pbaspect([1 0.3 1]); box off; set(gca, 'TickDir', 'out');
    title(theStimuli{d});
end
set(theFig, 'PaperPosition', [0 0 10 10]);
set(theFig, 'PaperSize', [10 10]);
set(theFig, 'Color', 'w');
set(theFig, 'InvertHardcopy', 'off');
saveas(theFig, outFig1, 'pdf');
close(theFig);

%% Plot the chromaticity
outFig2 = fullfile(outDir, 'FigureX_Chromaticity.pdf');
theRGB = jet(length(chromaticity));
theFig = figure;
for ii = 1:length(chromaticity)
    h(ii) = plot(chromaticity(1, ii), chromaticity(2, ii), 'Marker', 's', 'Color', 'k', 'MarkerFaceColor', theRGB(ii, :)); hold on;
end
% Plot horseshoe
load T_xyz1931
out = SplineCmf(S_xyz1931, T_xyz1931, S_xyz1931);
x = out(1, :)./sum(out);
y = out(2, :)./sum(out);
plot([x(1:65) x(1)], [y(1:65) y(1)], '-k');
xlim([0 0.9]); ylim([0 0.9]);
legend(h, theStimuli);
set(theFig, 'PaperPosition', [0 0 5 5]);
set(theFig, 'PaperSize', [5 5]);
set(theFig, 'Color', 'w');
set(theFig, 'InvertHardcopy', 'off');
saveas(theFig, outFig2, 'pdf');
xlabel('x'); ylabel('y');

%close(theFig);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');

cd(currDir);