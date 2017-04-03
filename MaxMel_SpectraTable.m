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


theFig = figure;

wls = SToWls([380 2 201]);
load T_xyz1931
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,WlsToS(wls));

currDir = pwd;
LOG_PLOT = true;
if LOG_PLOT
    outFig1 = fullfile(outDir, 'FigureX_Spectra_log10.pdf');
else
    outFig1 = fullfile(outDir, 'FigureX_Spectra.pdf');
end
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
    
    theSpd = [];
    
    if ~isempty(theMATFile)
        % Load the MAT file
        tmp = load(theMATFile.name);
        
        % Get the background spectrum
        bgSpd = tmp.cals{1}.modulationBGMeas.predictedSpd; hold on;ff
        % Get the background spectrum
        bgSpd = tmp.cals{1}.modulationBGMeas.meas.pr650.spectrum; hold on;
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
    end
end
cd(currDir);