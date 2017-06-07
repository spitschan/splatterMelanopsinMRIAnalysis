function splatterMel_SpectralPlots(ppsRawDataDir, analysisDir)

fprintf('> Running %s\n', mfilename);

% Define the paths
outDir = fullfile(analysisDir, 'figures');

theDataPaths = {'MaxLMS400pct/HERO_asb1/040716/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
    'MaxLMSCRF/HERO_asb1/060816/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MaxMel400pct/HERO_asb1/032416/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
    'MaxMelCRF/HERO_asb1/060716/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'SplatterControlCRF/HERO_asb1/051016/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation'};

theStimuli = {'LMS 400%' 'LMS CRF' 'Mel 400%' 'Mel CRF' 'Splatter CRF'};
close all;


theFig = figure;

wls = SToWls([380 2 201]);
load T_xyz1931
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,WlsToS(wls));

currDir = pwd;
LOG_PLOT = false;
if LOG_PLOT
    outFig1 = fullfile(outDir, 'FigureX_Spectra_log10.pdf');
else
    outFig1 = fullfile(outDir, 'FigureX_Spectra.pdf');
end
for d = 1:length(theDataPaths)
    subplot(3, 2, d);
    dataPath = theDataPaths{d};
    
    % Find the folders
    theFolders = dir(fullfile(ppsRawDataDir, dataPath));
    
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
    cd(fullfile(ppsRawDataDir, dataPath, theFolders(1).name));
    
    % Find the only MAT file there is going to be
    theMATFile = dir([pwd '/*.mat']);
    
    theSpd = [];
    
    if ~isempty(theMATFile)
        % Load the MAT file
        tmp = load(theMATFile.name);
        
        % Get the background spectrum
        bgSpd = tmp.cals{1}.modulationBGMeas.predictedSpd; hold on;
        
        % Extract the chromaticity
        chromaticity{d}(1, :) = (T_xyz([1 2], :)*bgSpd)/sum((T_xyz*bgSpd));
        
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
            chromaticity{d}(kk, :) = (T_xyz([1 2], :)*modSpd)/sum((T_xyz*modSpd));
            % ... and plot it.
            if LOG_PLOT
                plot(wls, log10(modSpd), 'Color', theRGB(kk, :)); hold on;
            else
                plot(wls, modSpd, 'Color', theRGB(kk, :)); hold on;
            end
            
        end
        % Plot the background spectrum
        if LOG_PLOT
            plot(wls, log10(bgSpd), '-k'); hold on;
        else
            plot(wls, bgSpd, '-k'); hold on;
        end
    end
    
    % Set some plot properties
    xlabel('Wavelength [nm]');
    if LOG_PLOT
        ylabel('Radiance [log W/m2/sr/nm]');
    else
        ylabel('Radiance [W/m2/sr/nm]');
    end
    xlim([380 780]);
    if LOG_PLOT
        ylim([-6 0]);
    else
        if d < 5
            ylim([0 0.03]);
        end
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
c = 1;
for dd = 1:length(theDataPaths);
        h(c) = plot(chromaticity{dd}(1, 1), chromaticity{dd}(1, 2), 'Marker', 's', 'Color', 'k', 'MarkerFaceColor', theRGB(dd, :)); hold on;
    c = c+1;
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

outFig2 = fullfile(outDir, 'FigureX_ChromaticityCloseup.pdf');
theFig = figure;
c = 1;
for dd = [2 4 5]
    plot(chromaticity{dd}(:, 1), chromaticity{dd}(:, 2), '-k'); hold on;
    for ii = 1:length(chromaticity{dd})
        h(c) = plot(chromaticity{dd}(ii, 1), chromaticity{dd}(ii, 2), 'Marker', 's', 'Color', 'k', 'MarkerFaceColor', theRGB(dd, :)); hold on;
    end
    plot(chromaticity{dd}(1, 1), chromaticity{dd}(1, 2), 'sk', 'MarkerFaceColor', 'k');
    c = c+1;
end

% Plot horseshoe
load T_xyz1931
out = SplineCmf(S_xyz1931, T_xyz1931, S_xyz1931);
x = out(1, :)./sum(out);
y = out(2, :)./sum(out);
plot([x(1:65) x(1)], [y(1:65) y(1)], '-k');
xlim([0.4 0.7]); ylim([0.25 0.55]);
legend(h, theStimuli{[2 4 5]});
set(theFig, 'PaperPosition', [0 0 5 5]);
set(theFig, 'PaperSize', [5 5]);
set(theFig, 'Color', 'w');
set(theFig, 'InvertHardcopy', 'off');
saveas(theFig, outFig2, 'pdf');
xlabel('x'); ylabel('y');

%close(theFig);
pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');

cd(currDir);