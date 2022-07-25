
function plotIntensityBlebCompare(p)

% plotIntensityBlebCompare - plots some bleb paper-specific plots using data calculated by the IntensityBlebCompare Process


% check inputs
if ~isfield(p, 'analyzeOtherChannel'), p.analyzeOtherChannel = 0; end
assert(p.analyzeOtherChannel == 0 | p.analyzeOtherChannel == 1, 'p.analyzeOtherChannel must be either left unset or set to 0 or 1');

if ~isfield(p, 'analyzeDiffusion'), p.analyzeDiffusion = 0; end
assert(p.analyzeDiffusion == 0 | p.analyzeDiffusion == 1, 'p.analyzeDiffusion must be either left unset or set to 0 or 1');

if ~isfield(p, 'analyzeDistance'), p.analyzeDistance = 0; end
assert(p.analyzeDistance == 0 | p.analyzeDistance == 1, 'p.analyzeDistance must be either left unset or set to 0 or 1');

if ~isfield(p, 'analyzeForwardsMotion'), p.analyzeForwardsMotion = 0; end
assert(p.analyzeForwardsMotion == 0 | p.analyzeForwardsMotion == 1, 'p.analyzeForwardsMotion must be either left unset or set to 0 or 1');

if ~isfield(p, 'analyzeVonMises'), p.analyzeVonMises = 0; end
assert(p.analyzeVonMises == 0 | p.analyzeVonMises == 1, 'p.analyzeVonMises must be either left unset or set to 0 or 1');

% make a directory to save data in
if ~isfolder(p.savePath), system(['mkdir -p ' p.savePath]); end

% initialize variables
meanIntensityBlebCells = []; meanIntensityNotBlebCells = [];
maxIntensityBlebCells = []; maxIntensityNotBlebCells = [];
isProtrusionPatches = []; isCertainProtrusionPatches = [];
meanIntensityPatches = []; maxIntensityPatches = [];
meanMotionPatches = []; 
volumePatches = [];
isProtrusionFaces = []; isCertainProtrusionFaces = [];
intensityFaces = []; motionFaces = []; curvatureFaces = [];
motifVolumeFaces = [];
if p.analyzeOtherChannel == 1
    meanIntensityOtherPatches = []; maxIntensityOtherPatches = [];
    intensityOtherFaces = [];
end
if p.analyzeDiffusion == 1
    diffusedProtrusionsFaces = [];
    diffusedSVMFaces = [];
end
if p.analyzeDistance == 1
    distanceFaces = [];
end
if p.analyzeForwardsMotion == 1
    meanForwardsMotionPatches = [];
    forwardsMotionFaces = [];
end
if p.analyzeVonMises == 1
    vonMisesBlebs = [];
    vonMisesBlebCenters = [];
    vonMisesBlebCentersLargeBlebs = [];
    vonMisesBlebsRand = [];
    vonMisesSurface = [];
    vonMisesIntensity = [];
    vonMisesIntensityOther = [];
    vonMisesIntensityControl = [];
    vonMisesIntensityOtherControl = [];
    vonMisesHighIntensity = [];
    vonMisesIntensityPixel = [];
    vonMisesNegCurvature = [];
    vonMisesIntensityMin = [];
    vonMisesMotion = [];
    vonMisesPosMotion = [];
    vonMisesIntensityDiscrete = [];
    vonMisesIntensityRandDiscrete = [];
end
intensityDistIntensityMatrix = []; intensityDistIntensityBackgroundMatrix = [];intensityDistDistMatrix = [];
cellIndexPatches = []; cellIndexFaces = [];

% iterate through the list of cells
for s = 1:length(p.cellsList)
    
    % display progress
    disp(['Cell ' num2str(s) ' of ' num2str(length(p.cellsList))]);
    
    % load the bleb and intensity statistics
    analysisPath = fullfile(p.mainDirectory, p.cellsList{s}, 'Morphology', 'Analysis', 'IntensityBlebCompare');
    statsStruct = load(fullfile(analysisPath, 'stats.mat'));
    comparePatches = statsStruct.comparePatches;
    compareFaces = statsStruct.compareFaces;
    vonMises = statsStruct.vonMises;
    %intensityAsDist = statsStruct.intensityAsDist;
    convert = statsStruct.convert; 
    compareFrames = statsStruct.compareFrames;
    
    % append and subset data in the comparePatches variable
    isProtrusionPatches = [isProtrusionPatches, comparePatches.isProtrusion];
    isCertainProtrusionPatches = [isCertainProtrusionPatches, comparePatches.isCertainProtrusion];
    meanIntensityPatches = [meanIntensityPatches, comparePatches.meanIntensity];
    maxIntensityPatches = [maxIntensityPatches, comparePatches.maxIntensity];
    meanMotionPatches = [meanMotionPatches, convert.motion*comparePatches.meanMotion];
    volumePatches = [volumePatches, convert.volume*comparePatches.volume];
    if p.analyzeOtherChannel == 1
        meanIntensityOtherPatches = [meanIntensityOtherPatches, comparePatches.meanIntensityOther];
        maxIntensityOtherPatches = [maxIntensityOtherPatches, comparePatches.maxIntensityOther];
    end
    if p.analyzeForwardsMotion == 1
        meanForwardsMotionPatches = [meanForwardsMotionPatches, convert.motion*comparePatches.meanForwardsMotion];
    end
    
    % append and subset data in the compareFaces variable
    isProtrusionFaces = [isProtrusionFaces, compareFaces.isProtrusion];
    isCertainProtrusionFaces = [isCertainProtrusionFaces, compareFaces.isCertainProtrusion];
    intensityFaces = [intensityFaces, compareFaces.intensityNormal];
    %intensityFaces = [intensityFaces, zscore(compareFaces.intensity(:))']; % used for optical flow
    motionFaces = [motionFaces, convert.motion*compareFaces.motion];
    curvatureFaces = [curvatureFaces, convert.curvature*compareFaces.curvature];
    motifVolumeFaces = [motifVolumeFaces, convert.volume*compareFaces.motifVolume];
    if p.analyzeOtherChannel == 1
        intensityOtherFaces = [intensityOtherFaces, compareFaces.intensityOtherNormal];
    end
    if p.analyzeDiffusion == 1
        diffusedProtrusionsFaces = [diffusedProtrusionsFaces, compareFaces.diffusedProtrusions];
        diffusedSVMFaces = [diffusedSVMFaces, compareFaces.diffusedSVM];
    end
    if p.analyzeForwardsMotion == 1
        forwardsMotionFaces = [forwardsMotionFaces, convert.motion*compareFaces.forwardsMotion];
    end
    if p.analyzeDistance == 1
        distanceFaces = [distanceFaces, convert.edgeLength.*compareFaces.distanceTransformProtrusions];
    end
    
    % append and subset data in the compareFrames variable
    convert.motion = convert.motion*60/30.3; % should be temp !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     compareFrames.meanMotion = compareFrames.meanMotion*convert.motion;
     compareFrames.volBlebs = compareFrames.volBlebs*convert.volume;
    
    % append polarization data
    if p.analyzeVonMises == 1
        vonMisesBlebs = [vonMisesBlebs; vonMises.blebs];
        vonMisesBlebCenters = [vonMisesBlebCenters; vonMises.blebCenters];
%        vonMisesBlebCentersLargeBlebs = [vonMisesBlebCentersLargeBlebs; vonMises.blebCentersLargeBlebs];
        vonMisesBlebsRand = [vonMisesBlebsRand; vonMises.blebsRand];
        vonMisesSurface = [vonMisesSurface; vonMises.surface];
        vonMisesIntensity = [vonMisesIntensity; vonMises.intensity];
        vonMisesIntensityOther = [vonMisesIntensityOther; vonMises.intensityOther];
        vonMisesIntensityControl = [vonMisesIntensityControl; vonMises.intensityControl];
        vonMisesIntensityOtherControl = [vonMisesIntensityOtherControl; vonMises.intensityOtherControl];
        %vonMisesHighIntensity = [vonMisesHighIntensity; vonMises.highIntensity];
        vonMisesIntensityPixel = [vonMisesIntensityPixel; vonMises.intensityPixel];
        vonMisesNegCurvature = [vonMisesNegCurvature; vonMises.negCurvature];
        vonMisesIntensityMin = [vonMisesIntensityMin; vonMises.intensityMin];
        %vonMisesIntensityDiscrete = [vonMisesIntensityDiscrete; vonMises.intensityDiscrete];
        %vonMisesIntensityRandDiscrete = [vonMisesIntensityRandDiscrete; vonMises.intensityRandDiscrete];
        vonMisesMotion = [vonMisesMotion; vonMises.motion];
        %vonMisesPosMotion = [vonMisesPosMotion; vonMises.posMotion];
    end
    
%     % append intensityAsDist data
%     intensityDistIntensityMatrix = [intensityDistIntensityMatrix; intensityAsDist.intensity./mean(intensityAsDist.intensity, 'omitnan')];
%     intensityDistIntensityBackgroundMatrix = [intensityDistIntensityBackgroundMatrix; intensityAsDist.intensityBackgroundSubtract./mean(intensityAsDist.intensityBackgroundSubtract, 'omitnan')];
%     intensityDistDistMatrix = [intensityDistDistMatrix; intensityAsDist.dist];
    
    % keep track of the different cells
    cellIndexPatches = [cellIndexPatches, s.*ones(1,length(comparePatches.isProtrusion))];
    cellIndexFaces = [cellIndexFaces, s.*ones(1,length(compareFaces.isProtrusion))];
    
    % find the mean patch intensity for blebs and non-blebs
    meanIntensityBlebCells = [meanIntensityBlebCells, mean(comparePatches.meanIntensity(comparePatches.isProtrusion == 1))];
    meanIntensityNotBlebCells = [meanIntensityNotBlebCells, mean(comparePatches.meanIntensity(comparePatches.isProtrusion == 0))];
    maxIntensityBlebCells = [maxIntensityBlebCells, mean(comparePatches.maxIntensity(comparePatches.isProtrusion == 1))];
    maxIntensityNotBlebCells = [maxIntensityNotBlebCells, mean(comparePatches.maxIntensity(comparePatches.isProtrusion == 0))];

    sum(comparePatches.isProtrusion)
end
intensityFaces = intensityFaces - prctile(intensityFaces(:),1); % used for the optical flow, was 5

% check to see if the motion plots should be made
if sum(isfinite(motionFaces)) > 0
    makeMotionPlots = 1;
else
    makeMotionPlots = 0;
end

disp(['Total number of blebs '  num2str(sum(isProtrusionPatches))])

%% Find the mean intensity as a function of surface curvature and vice-versa (faces)
numGridPointsCurvature = 512; % 128 pip2
curvatureIndex = 1 + round((curvatureFaces - min(curvatureFaces)) ./ (max(curvatureFaces)-min(curvatureFaces))*(numGridPointsCurvature-1));
meanIntensity = nan(1,numGridPointsCurvature);
for i = 1:numGridPointsCurvature
    meanIntensity(1,i) = nanmean(intensityFaces(distanceFaces<0.5 & curvatureIndex==i));
end
f = figure;
plot(linspace(min(curvatureFaces), max(curvatureFaces), numGridPointsCurvature), meanIntensity, 'LineWidth', 2)
%axis([-3 3 -Inf Inf]);
xlabel('Curvature'); ylabel('Mean Intensity')
title('Mean Intensity as a Function of Surface Curvature For All Faces');
saveName = fullfile(p.savePath, 'intensityFunctionOfCurvatureFaces');
saveas(f, saveName, 'epsc'); savefig(f, saveName);

%% Find the mean intensity as a function of surface curvature with error (faces)
numCells = max(cellIndexFaces(:));
numGridPointsCurvature = 512; % 128 pip2
curvatureIndex = 1 + round((curvatureFaces - min(curvatureFaces)) ./ (max(curvatureFaces)-min(curvatureFaces))*(numGridPointsCurvature-1));
meanIntensity = nan(numCells,numGridPointsCurvature);
for i = 1:numGridPointsCurvature
    for nC = 1:numCells
        meanIntensity(nC,i) = nanmean(intensityFaces(cellIndexFaces==nC & curvatureIndex==i));
    end
end
meanIntensity_mean = nanmean(meanIntensity, 1);
%meanIntensity_mean = meanIntensity_mean(isfinite(meanIntensity_mean)); %this should maybe be removed from the other code too!
meanIntensity_ci = 1.96*nanstd(meanIntensity, 1)./sqrt(sum(~isnan(meanIntensity)));
f = figure;
xOffset = (1/(2*numGridPointsCurvature))*(max(curvatureFaces(isfinite(curvatureFaces))));
xPoints = linspace(min(curvatureFaces), max(curvatureFaces(isfinite(curvatureFaces))), numGridPointsCurvature) + xOffset;
plot(xPoints(1:length(meanIntensity_mean)), meanIntensity_mean, 'LineWidth', 2, 'Color', [0.5, 0, 0.5]); 
hold on
xPointsPatch = [xPoints(1:length(meanIntensity_mean)), fliplr(xPoints(1:length(meanIntensity_mean))), xPoints(1)];
patch = fill(xPointsPatch, [meanIntensity_mean + meanIntensity_ci, fliplr(meanIntensity_mean - meanIntensity_ci), meanIntensity_mean(1) + meanIntensity_ci(1)], [0.5, 0, 0.5]);
%set(patch, 'edgecolor', 'none');
set(patch, 'facecolor', [0.5,0,0.5]);
%set(patch, 'FaceAlpha', 0.5);
axis([-2 2 0 2]);
hold on;
plot(xPoints(1:length(meanIntensity_mean)), meanIntensity_mean, 'LineWidth', 2, 'Color', [0.5, 0, 0.5])
xlabel('Curvature'); ylabel('Mean Intensity')
title('Mean Intensity as a Function of Surface Curvature For All Faces');
saveName = fullfile(p.savePath, 'intensityFunctionOfCurvatureFacesError');
saveas(f, saveName, 'epsc'); savefig(f, saveName);
saveName = fullfile(p.savePath, 'meanIntensity_mean');
save([saveName '.txt'],'meanIntensity_mean','-ascii');
saveName = fullfile(p.savePath, 'meanIntensity_ci');
save([saveName '.txt'],'meanIntensity_ci','-ascii');
forXaxis = xPoints(1:length(meanIntensity_mean));
saveName = fullfile(p.savePath, 'forXaxisCurv');
save([saveName '.txt'],'forXaxis','-ascii');

% 
% numGridPointsIntensity = 128;
% intensityIndex = 1 + round((intensityFaces - min(intensityFaces)) ./ (max(intensityFaces)-min(intensityFaces))*(numGridPointsIntensity-1));
% meanCurvature = nan(1,numGridPointsIntensity);
% for i = 1:numGridPointsIntensity
%     meanCurvature(1,i) = mean(curvatureFaces(intensityIndex==i));
% end
% f = figure;
% plot(linspace(min(intensityFaces), max(intensityFaces), numGridPointsIntensity), meanCurvature, 'LineWidth', 2)
% xlabel('Normalized Intensity'); ylabel('Curvature')
% title('Curvature as a Function of Surface Intensity For All Faces');
% saveName = fullfile(p.savePath, 'curvatureFunctionOfIntensityFaces');
% saveas(f, saveName, 'epsc'); savefig(f, saveName);

% plot the distance transform data
if p.analyzeDistance == 1 
    
%     % plot distances on all faces
%     numGridPointsDistance = 64;
%     distanceIndex = 1 + round((distanceFaces) ./ (max(distanceFaces(isfinite(distanceFaces))))*(numGridPointsDistance-1));
%     distanceIndex(~isfinite(distanceIndex)) = NaN;
%     meanIntensity = nan(1,numGridPointsDistance);
%     for i = 1:numGridPointsDistance
%         meanIntensity(1,i) = nanmean(intensityFaces(distanceIndex==i));
%     end
%     f = figure;
%     xOffset = (1/(2*numGridPointsDistance) )*(max(distanceFaces(isfinite(distanceFaces))));
%     plot(linspace(0, max(distanceFaces(isfinite(distanceFaces))), numGridPointsDistance) + xOffset, meanIntensity, 'LineWidth', 2)
%     xlabel('Distance from bleb edge'); ylabel('Mean Intensity')
%     title('Mean Intensity as a Function of Distance from Bleb Edge For All Faces');
%     saveName = fullfile(p.savePath, 'intensityFunctionOfDistanceFaces');
%     saveas(f, saveName, 'epsc'); savefig(f, saveName);
%     
%     % plot distances on and off blebs
%     numGridPointsDistance = 64;
%     distanceIndex = 1 + round((distanceFaces) ./ (max(distanceFaces(isfinite(distanceFaces))))*(numGridPointsDistance-1));
%     distanceIndex(~isfinite(distanceIndex)) = NaN;
%     meanIntensityOnBlebs = nan(1,numGridPointsDistance);
%     meanIntensityOffBlebs = nan(1,numGridPointsDistance);
%     for i = 1:numGridPointsDistance
%         meanIntensityOnBlebs(1,i) = nanmean(intensityFaces(distanceIndex==i & isProtrusionFaces==1));
%         meanIntensityOffBlebs(1,i) = nanmean(intensityFaces(distanceIndex==i & isProtrusionFaces==0));
%     end
%     f = figure;
%     xOffset = (1/(2*numGridPointsDistance) )*(max(distanceFaces(isfinite(distanceFaces))));
%     plot(linspace(0, max(distanceFaces(isfinite(distanceFaces))), numGridPointsDistance) + xOffset, meanIntensityOnBlebs, 'LineWidth', 2)
%     hold on
%     plot(linspace(0, max(distanceFaces(isfinite(distanceFaces))), numGridPointsDistance) + xOffset, meanIntensityOffBlebs, 'LineWidth', 2)
%     xlabel('Distance from bleb edge'); ylabel('Mean Intensity')
%     title('Mean Intensity as a Function of Distance from Bleb Edge on and off Blebs');
%     legend('On blebs', 'Off blebs');
%     saveName = fullfile(p.savePath, 'intensityFunctionOfDistanceBlebsFaces');
%     saveas(f, saveName, 'epsc'); savefig(f, saveName);
    
    % plot intensity distances on and off blebs with confidence intervals
    numCells = max(cellIndexFaces(:));
    numGridPointsDistance = 64;
    distanceIndex = 1 + round((distanceFaces) ./ (max(distanceFaces(isfinite(distanceFaces))))*(numGridPointsDistance-1));
    distanceIndex(~isfinite(distanceIndex)) = NaN;
    meanIntensityOnBlebs = nan(1,numGridPointsDistance);
    meanIntensityOffBlebs = nan(1,numGridPointsDistance);
    %meanIntensityOnBlebs = nan(numCells,numGridPointsDistance);
    %meanIntensityOffBlebs = nan(numCells,numGridPointsDistance);
    for i = 1:numGridPointsDistance
        for nC = 1:numCells

            meanIntensityOnBlebs(nC,i) = nanmean(intensityFaces(distanceIndex==i & isProtrusionFaces==1 & cellIndexFaces==nC));
            meanIntensityOffBlebs(nC,i) = nanmean(intensityFaces(distanceIndex==i & isProtrusionFaces==0 & cellIndexFaces==nC));
        end
    end
    meanIntensityOnBlebs_mean = nanmean(meanIntensityOnBlebs, 1); %meanIntensityOnBlebs_mean = meanIntensityOnBlebs_mean(isfinite(meanIntensityOnBlebs_mean));
    meanIntensityOffBlebs_mean = nanmean(meanIntensityOffBlebs, 1); %meanIntensityOffBlebs_mean = meanIntensityOffBlebs_mean(isfinite(meanIntensityOffBlebs_mean));
    meanIntensityOnBlebs_ci = 1.96*nanstd(meanIntensityOnBlebs, 1)./sqrt(sum(~isnan(meanIntensityOnBlebs))); %meanIntensityOnBlebs_ci = meanIntensityOnBlebs_ci(isfinite(meanIntensityOnBlebs_ci));
    meanIntensityOffBlebs_ci = 1.96*nanstd(meanIntensityOffBlebs, 1)./sqrt(sum(~isnan(meanIntensityOffBlebs))); %meanIntensityOffBlebs_ci = meanIntensityOffBlebs_ci(isfinite(meanIntensityOffBlebs_ci));
    f = figure;
    xOffset = (1/(2*numGridPointsDistance) )*(max(distanceFaces(isfinite(distanceFaces))));
    xPoints = linspace(0, max(distanceFaces(isfinite(distanceFaces))), numGridPointsDistance) + xOffset;
    plot(xPoints(1:length(meanIntensityOnBlebs_mean)), meanIntensityOnBlebs_mean, 'LineWidth', 2, 'Color', [0.5, 0, 0.5]); hold on
    plot(xPoints(1:length(meanIntensityOffBlebs_mean)), meanIntensityOffBlebs_mean, 'LineWidth', 2, 'Color', [0, 0, 0.6])
    xPointsPatch = [xPoints(1:length(meanIntensityOnBlebs_mean)), fliplr(xPoints(1:length(meanIntensityOnBlebs_mean))), xPoints(1)];
    patch = fill(xPointsPatch, [meanIntensityOnBlebs_mean + meanIntensityOnBlebs_ci, fliplr(meanIntensityOnBlebs_mean - meanIntensityOnBlebs_ci), meanIntensityOnBlebs_mean(1) + meanIntensityOnBlebs_ci(1)], [0.5, 0, 0.5]);
    %set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', 0);
    hold on;
    xPointsPatch = [xPoints(1:length(meanIntensityOffBlebs_mean)), fliplr(xPoints(1:length(meanIntensityOffBlebs_mean))), xPoints(1)];
    patch = fill(xPointsPatch, [meanIntensityOffBlebs_mean + meanIntensityOffBlebs_ci, fliplr(meanIntensityOffBlebs_mean - meanIntensityOffBlebs_ci), meanIntensityOffBlebs_mean(1) + meanIntensityOffBlebs_ci(1)], [0, 0, 0.6]);
    %set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', 0);
    plot(xPoints(1:length(meanIntensityOnBlebs_mean)), meanIntensityOnBlebs_mean, 'LineWidth', 2, 'Color', [0.5, 0, 0.5])
    plot(xPoints(1:length(meanIntensityOffBlebs_mean)), meanIntensityOffBlebs_mean, 'LineWidth', 2, 'Color', [0, 0, 0.6])
    xlabel('Distance from bleb edge'); ylabel('Mean Intensity')
    title('Mean Intensity as a Function of Distance from Bleb Edge on and off Blebs');
    legend('On blebs', 'Off blebs');
    saveName = fullfile(p.savePath, 'intensityFunctionOfDistanceBlebsFaces');
    saveas(f, saveName, 'epsc'); savefig(f, saveName);
    saveName = fullfile(p.savePath, 'meanIntensityOnBlebs_mean');
    save([saveName '.txt'],'meanIntensityOnBlebs_mean','-ascii');
    saveName = fullfile(p.savePath, 'meanIntensityOffBlebs_mean');
    save([saveName '.txt'],'meanIntensityOffBlebs_mean','-ascii');
    saveName = fullfile(p.savePath, 'meanIntensityOnBlebs_ci');
    save([saveName '.txt'],'meanIntensityOnBlebs_ci','-ascii');
    saveName = fullfile(p.savePath, 'meanIntensityOffBlebs_ci');
    save([saveName '.txt'],'meanIntensityOffBlebs_ci','-ascii');
    forXaxis = xPoints(1:length(meanIntensityOffBlebs_mean));
    saveName = fullfile(p.savePath, 'forXaxis');
    save([saveName '.txt'],'forXaxis','-ascii');
    1;
 
    % this may not work!!!!!!!
    % plot intensity as a function of curvature on and off blebs with confidence intervals
    numCells = max(cellIndexFaces(:));
    numGridPointsCurvature = 64;
    curvatureFacesForIndex = curvatureFaces;
    curvatureFacesForIndex(curvatureFacesForIndex<-3) = -3;
    curvatureFacesForIndex(curvatureFacesForIndex>3) = 3;
    curvatureFacesForIndex = curvatureFacesForIndex + 3;
    curvatureIndex = 1 + round((numGridPointsCurvature-1)*(curvatureFacesForIndex/6));
    curvatureIndex(~isfinite(curvatureIndex)) = NaN;
    meanIntensityOnBlebs = nan(numCells,numGridPointsCurvature);
    meanIntensityOffBlebs = nan(numCells,numGridPointsCurvature);
    for i = 1:numGridPointsCurvature
        for nC = 1:numCells
            meanIntensityOnBlebs(nC,i) = nanmean(intensityFaces(curvatureIndex==i & isProtrusionFaces==1 & cellIndexFaces==nC));
            meanIntensityOffBlebs(nC,i) = nanmean(intensityFaces(curvatureIndex==i & isProtrusionFaces==0 & cellIndexFaces==nC));
        end
    end
    meanIntensityOnBlebs_mean = nanmean(meanIntensityOnBlebs, 1); meanIntensityOnBlebs_mean = meanIntensityOnBlebs_mean(isfinite(meanIntensityOnBlebs_mean));
    meanIntensityOffBlebs_mean = nanmean(meanIntensityOffBlebs, 1); meanIntensityOffBlebs_mean = meanIntensityOffBlebs_mean(isfinite(meanIntensityOffBlebs_mean));
    meanIntensityOnBlebs_ci = 1.96*nanstd(meanIntensityOnBlebs, 1)./sqrt(sum(~isnan(meanIntensityOnBlebs))); meanIntensityOnBlebs_ci = meanIntensityOnBlebs_ci(isfinite(meanIntensityOnBlebs_ci));
    meanIntensityOffBlebs_ci = 1.96*nanstd(meanIntensityOffBlebs, 1)./sqrt(sum(~isnan(meanIntensityOffBlebs))); meanIntensityOffBlebs_ci = meanIntensityOffBlebs_ci(isfinite(meanIntensityOffBlebs_ci));
    f = figure;
    xOffset = (1/(2*numGridPointsCurvature) )*(max(curvatureFacesForIndex(isfinite(curvatureFacesForIndex))))-3;
    xPoints = linspace(0, max(curvatureFacesForIndex(isfinite(curvatureFacesForIndex))), numGridPointsCurvature) + xOffset;
    plot(xPoints(1:length(meanIntensityOnBlebs_mean)), meanIntensityOnBlebs_mean, 'LineWidth', 2, 'Color', [0.5, 0, 0.5]); hold on
    plot(xPoints(1:length(meanIntensityOffBlebs_mean)), meanIntensityOffBlebs_mean, 'LineWidth', 2, 'Color', [0, 0, 0.6])
    xPointsPatch = [xPoints(1:length(meanIntensityOnBlebs_mean)), fliplr(xPoints(1:length(meanIntensityOnBlebs_mean))), xPoints(1)];
    patch = fill(xPointsPatch, [meanIntensityOnBlebs_mean + meanIntensityOnBlebs_ci, fliplr(meanIntensityOnBlebs_mean - meanIntensityOnBlebs_ci), meanIntensityOnBlebs_mean(1) + meanIntensityOnBlebs_ci(1)], [0.5, 0, 0.5]);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', 0.5);
    hold on;
    xPointsPatch = [xPoints(1:length(meanIntensityOffBlebs_mean)), fliplr(xPoints(1:length(meanIntensityOffBlebs_mean))), xPoints(1)];
    patch = fill(xPointsPatch, [meanIntensityOffBlebs_mean + meanIntensityOffBlebs_ci, fliplr(meanIntensityOffBlebs_mean - meanIntensityOffBlebs_ci), meanIntensityOffBlebs_mean(1) + meanIntensityOffBlebs_ci(1)], [0, 0, 0.6]);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', 0.5);
    plot(xPoints(1:length(meanIntensityOnBlebs_mean)), meanIntensityOnBlebs_mean, 'LineWidth', 2, 'Color', [0.5, 0, 0.5])
    plot(xPoints(1:length(meanIntensityOffBlebs_mean)), meanIntensityOffBlebs_mean, 'LineWidth', 2, 'Color', [0, 0, 0.6])
    xlabel('Mean Curvature'); ylabel('Mean Intensity')
    title('Mean Intensity as a Function of Mean Curvature on and off Blebs');
    legend('On blebs', 'Off blebs');
    saveName = fullfile(p.savePath, 'intensityFunctionOfCurvatureBlebsFaces');
    saveas(f, saveName, 'epsc'); savefig(f, saveName);
    saveName = fullfile(p.savePath, 'meanIntensityOnBlebsCurv_mean');
    save([saveName '.txt'],'meanIntensityOnBlebs_mean','-ascii');
    saveName = fullfile(p.savePath, 'meanIntensityOffBlebsCurv_mean');
    save([saveName '.txt'],'meanIntensityOffBlebs_mean','-ascii');
    saveName = fullfile(p.savePath, 'meanIntensityOnBlebsCurv_ci');
    save([saveName '.txt'],'meanIntensityOnBlebs_ci','-ascii');
    saveName = fullfile(p.savePath, 'meanIntensityOffBlebsCurv_ci');
    save([saveName '.txt'],'meanIntensityOffBlebs_ci','-ascii');
    forXaxis = xPoints(1:length(meanIntensityOffBlebs_mean));
    saveName = fullfile(p.savePath, 'forXaxisCurvature');
    save([saveName '.txt'],'forXaxis','-ascii');
    
   % plot curvature distances on and off blebs with confidence intervals
    numCells = max(cellIndexFaces(:));
    numGridPointsDistance = 64;
    distanceIndex = 1 + round((distanceFaces) ./ (max(distanceFaces(isfinite(distanceFaces))))*(numGridPointsDistance-1));
    distanceIndex(~isfinite(distanceIndex)) = NaN;
    meanCurvatureOnBlebs = nan(numCells,numGridPointsDistance);
    meanCurvatureOffBlebs = nan(numCells,numGridPointsDistance);
    %volumeBlebs = volumePatches(logical(isProtrusionPatches));
    for i = 1:numGridPointsDistance
        for nC = 1:numCells
            meanCurvatureOnBlebs(nC,i) = nanmean(curvatureFaces(distanceIndex==i & isProtrusionFaces==1 & cellIndexFaces==nC));
            meanCurvatureOffBlebs(nC,i) = nanmean(curvatureFaces(distanceIndex==i & isProtrusionFaces==0 & cellIndexFaces==nC));
        end
    end
    meanCurvatureOnBlebs_mean = nanmean(meanCurvatureOnBlebs, 1); %meanCurvatureOnBlebs_mean = meanCurvatureOnBlebs_mean(isfinite(meanCurvatureOnBlebs_mean));
    meanCurvatureOffBlebs_mean = nanmean(meanCurvatureOffBlebs, 1); %meanCurvatureOffBlebs_mean = meanCurvatureOffBlebs_mean(isfinite(meanCurvatureOffBlebs_mean));
    meanCurvatureOnBlebs_ci = 1.96*nanstd(meanCurvatureOnBlebs, 1)./sqrt(sum(~isnan(meanCurvatureOnBlebs))); %meanCurvatureOnBlebs_ci = meanCurvatureOnBlebs_ci(isfinite(meanCurvatureOnBlebs_ci));
    meanCurvatureOffBlebs_ci = 1.96*nanstd(meanCurvatureOffBlebs, 1)./sqrt(sum(~isnan(meanCurvatureOffBlebs))); %meanCurvatureOffBlebs_ci = meanCurvatureOffBlebs_ci(isfinite(meanCurvatureOffBlebs_ci));
    f = figure;
    xOffset = (1/(2*numGridPointsDistance) )*(max(distanceFaces(isfinite(distanceFaces))));
    xPoints = linspace(0, max(distanceFaces(isfinite(distanceFaces))), numGridPointsDistance) + xOffset;
    plot(xPoints(1:length(meanCurvatureOnBlebs_mean)), meanCurvatureOnBlebs_mean, 'LineWidth', 2, 'Color', [0.5, 0, 0.5]); hold on
    plot(xPoints(1:length(meanCurvatureOffBlebs_mean)), meanCurvatureOffBlebs_mean, 'LineWidth', 2, 'Color', [0, 0, 0.6])
    xPointsPatch = [xPoints(1:length(meanCurvatureOnBlebs_mean)), fliplr(xPoints(1:length(meanCurvatureOnBlebs_mean))), xPoints(1)];
    patch = fill(xPointsPatch, [meanCurvatureOnBlebs_mean + meanCurvatureOnBlebs_ci, fliplr(meanCurvatureOnBlebs_mean - meanCurvatureOnBlebs_ci), meanCurvatureOnBlebs_mean(1) + meanCurvatureOnBlebs_ci(1)], [0.5, 0, 0.5]);
    %set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', 0);
    hold on;
    xPointsPatch = [xPoints(1:length(meanCurvatureOffBlebs_mean)), fliplr(xPoints(1:length(meanCurvatureOffBlebs_mean))), xPoints(1)];
    patch = fill(xPointsPatch, [meanCurvatureOffBlebs_mean + meanCurvatureOffBlebs_ci, fliplr(meanCurvatureOffBlebs_mean - meanCurvatureOffBlebs_ci), meanCurvatureOffBlebs_mean(1) + meanCurvatureOffBlebs_ci(1)], [0, 0, 0.6]);
    %set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', 0);
    plot(xPoints(1:length(meanCurvatureOnBlebs_mean)), meanCurvatureOnBlebs_mean, 'LineWidth', 2, 'Color', [0.5, 0, 0.5])
    plot(xPoints(1:length(meanCurvatureOffBlebs_mean)), meanCurvatureOffBlebs_mean, 'LineWidth', 2, 'Color', [0, 0, 0.6])
    xlabel('Distance from bleb edge'); ylabel('Mean Curvature')
    title('Mean Curvature as a Function of Distance from Bleb Edge on and off Blebs');
    legend('On blebs', 'Off blebs');
    saveName = fullfile(p.savePath, 'curvatureFunctionOfDistanceBlebsFaces');
    saveas(f, saveName, 'epsc'); savefig(f, saveName);
    saveName = fullfile(p.savePath, 'meanCurvatureOnBlebs_mean');
    save([saveName '.txt'],'meanCurvatureOnBlebs_mean','-ascii');
    saveName = fullfile(p.savePath, 'meanCurvatureOffBlebs_mean');
    save([saveName '.txt'],'meanCurvatureOffBlebs_mean','-ascii');
    saveName = fullfile(p.savePath, 'meanCurvatureOnBlebs_ci');
    save([saveName '.txt'],'meanCurvatureOnBlebs_ci','-ascii');
    saveName = fullfile(p.savePath, 'meanCurvatureOffBlebs_ci');
    save([saveName '.txt'],'meanCurvatureOffBlebs_ci','-ascii');
    forXaxis = xPoints(1:length(meanCurvatureOffBlebs_mean));
    saveName = fullfile(p.savePath, 'forXaxis');
    save([saveName '.txt'],'forXaxis','-ascii');
    
    % plot distances for high intensity populations
    numGridPointsDistance = 64;
    stdCutoff = 0;
    distanceIndex = 1 + round((distanceFaces) ./ (max(distanceFaces(isfinite(distanceFaces))))*(numGridPointsDistance-1));
    distanceIndex(~isfinite(distanceIndex)) = NaN;
    countAllOn = nan(1,numGridPointsDistance); countAllOff = nan(1,numGridPointsDistance);
    countHighIntensityOn = nan(1,numGridPointsDistance); countHighIntensityOff = nan(1,numGridPointsDistance);
    meanIntensity = mean(intensityFaces); stdIntensity = std(intensityFaces);
    for i = 1:numGridPointsDistance
        countAllOn(1,i) = sum(distanceIndex==i & isProtrusionFaces==0);
        countAllOff(1,i) = sum(distanceIndex==i & isProtrusionFaces==0);
        countHighIntensityOn(1,i) = sum(distanceIndex==i & intensityFaces > (meanIntensity+stdCutoff*stdIntensity) & isProtrusionFaces==1);
        countHighIntensityOff(1,i) = sum(distanceIndex==i & intensityFaces > (meanIntensity+stdCutoff*stdIntensity) & isProtrusionFaces==0);
    end
    countHighIntensityOn = countHighIntensityOn./countAllOn;
    countHighIntensityOff = countHighIntensityOff./countAllOff;
    %countHighIntensity = countHighIntensity./sum(countHighIntensity);
    %countLowIntensity = countLowIntensity./sum(countLowIntensity);
    f = figure;
    xOffset = (1/(2*numGridPointsDistance) )*(max(distanceFaces(isfinite(distanceFaces))));
    plot(linspace(0, max(distanceFaces(isfinite(distanceFaces))), numGridPointsDistance) + xOffset, countHighIntensityOn, 'LineWidth', 2)
    hold on
    plot(linspace(0, max(distanceFaces(isfinite(distanceFaces))), numGridPointsDistance) + xOffset, countHighIntensityOff, 'LineWidth', 2)
    xlabel('Distance from bleb edge'); ylabel('Mean Intensity')
    title('High intensity: Normalized Count as a Function of Distance from Bleb Edge');
    legend('on bleb', 'off bleb');
    %axis([0 5 0 Inf]);
    saveName = fullfile(p.savePath, 'countFunctionOfDistanceHighFaces');
    saveas(f, saveName, 'epsc'); savefig(f, saveName);
    

    
    % plot curvature as a function of distances on and off blebs
    numGridPointsDistance = 64;
    distanceIndex = 1 + round((distanceFaces) ./ (max(distanceFaces(isfinite(distanceFaces))))*(numGridPointsDistance-1));
    distanceIndex(~isfinite(distanceIndex)) = NaN;
    meanCurvatureOnBlebs = nan(1,numGridPointsDistance);
    meanCurvatureOffBlebs = nan(1,numGridPointsDistance);
    for i = 1:numGridPointsDistance
        meanCurvatureOnBlebs(1,i) = nanmean(curvatureFaces(distanceIndex==i & isProtrusionFaces==1));
        meanCurvatureOffBlebs(1,i) = nanmean(curvatureFaces(distanceIndex==i & isProtrusionFaces==0));
    end
    f = figure;
    xOffset = (1/(2*numGridPointsDistance) )*(max(distanceFaces(isfinite(distanceFaces))));
    plot(linspace(0, max(distanceFaces(isfinite(distanceFaces))), numGridPointsDistance) + xOffset, meanCurvatureOnBlebs, 'LineWidth', 2)
    hold on
    plot(linspace(0, max(distanceFaces(isfinite(distanceFaces))), numGridPointsDistance) + xOffset, meanCurvatureOffBlebs, 'LineWidth', 2)
    xlabel('Distance from bleb edge'); ylabel('Mean Curvature')
    title('Mean Curvature as a Function of Distance from Bleb Edge on and off Blebs');
    legend('On blebs', 'Off blebs');
    saveName = fullfile(p.savePath, 'curvatureFunctionOfDistanceBlebsFaces');
    saveas(f, saveName, 'epsc'); savefig(f, saveName);
    1;

end

%% Plot polarization statistics
if p.analyzeVonMises == 1
    
    [numTotalRand,~] = size(vonMisesBlebsRand);
    [numTotalFrames,~] = size(vonMisesBlebs);
    
        % plot as cumsum line plots rather than line plots
    f = figure;
    binWidth = 0.02; binMax = 1;
    edges = -binMax:binWidth:binMax;
    edgeCenters = (edges(1:end-1)+edges(2:end))/2;
    randCor = 2*rand(2000,1)-1;
    motionIntensityDirection = histcounts(sum(vonMisesMotion(:,1:3).*vonMisesIntensity(:,1:3), 2), edges, 'Normalization', 'probability');
    %surfaceIntensityDirection = histcounts(sum(vonMisesSurface(:,1:3).*vonMisesIntensity(:,1:3), 2), edges, 'Normalization', 'probability');
    blebIntensityDirection = histcounts(sum(vonMisesBlebs(:,1:3).*vonMisesIntensity(:,1:3), 2), edges, 'Normalization', 'probability');
    blebRandIntensityDirection = histcounts(sum(vonMisesBlebsRand(:,1:3).*repelem(vonMisesIntensity(:,1:3), numTotalRand/numTotalFrames, 1),2), edges, 'Normalization', 'probability');
    randCorDirection = histcounts(randCor, edges, 'Normalization', 'probability');
    plot(edgeCenters, cumsum(motionIntensityDirection), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 12)
    hold on
    %plot(edgeCenters, cumsum(surfaceIntensityDirection), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 12)
    plot(edgeCenters, cumsum(blebIntensityDirection), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 12)
    plot(edgeCenters, cumsum(blebRandIntensityDirection), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 12)
    plot(edgeCenters, cumsum(randCorDirection), 'LineWidth', 2)
    legend({'motion-intensity', 'bleb-intensity', 'blebPerm-intensity', 'no cor'}); colormap(jet);
    xlabel('Directional Correlation'); ylabel('Frequency Sum');
    title('Cummulative Directional Correlation');
    %axis([-5 5 -inf inf])
    saveName = fullfile(p.savePath, 'blebsIntensityDirectionDotVonMisesCumLines');
    saveas(f, saveName, 'epsc'); savefig(f, saveName);
    blebby = sum(vonMisesBlebs(:,1:3).*vonMisesIntensity(:,1:3), 2);
    blebbyRand = sum(vonMisesBlebsRand(:,1:3).*repelem(vonMisesIntensity(:,1:3), numTotalRand/numTotalFrames, 1),2);
    [h,pValue]=kstest2(blebby, blebbyRand);
end