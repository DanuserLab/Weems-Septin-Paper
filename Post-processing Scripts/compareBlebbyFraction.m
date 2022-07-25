function compareCellStatsRevise(p)

% compareCellStats - compares distributions derived from two sets of cells


% accumulate data on sets of cells
stats1 = collectDataMovie(p.cellsList1, p.mainDirectory1);
stats2 = collectDataMovie(p.cellsList2, p.mainDirectory2);
% display motif percentage statistics
disp(['Motif percentage for ' p.name1 ': '])
disp(['   mean:' num2str(mean(stats1.percentMotif))])
disp(['   std: ' num2str(std(stats1.percentMotif))])

disp(['Motif percentage for ' p.name2 ': '])
disp(['   mean:' num2str(mean(stats2.percentMotif))])
disp(['   std: ' num2str(std(stats2.percentMotif))])

disp('One-sided ttest on the mean motif percentage');
[~,pVal,ci,~] = ttest2(stats1.percentMotif, stats2.percentMotif,'Tail','left');
disp(['   p value: ' num2str(pVal)]);
disp(['   confidence interval: ' num2str(ci(1)) ' to ' num2str(ci(2))]);
disp('kstest on the mean motif percentage');
[~,pVal] = kstest2(stats1.percentMotif, stats2.percentMotif);
disp(['   p value: ' num2str(pVal)]);
disp([' '])
%disp(['   effect size: ' num2str(d)]);

% make a percentMotif plot
figure
firstColor = [0,0,0.75];
secondColor = [0.75,0,0.75];
thirdColor = [0.75,0,0];
markerSize = 250;
alpha = 0.7;
lineWidth = 0.001;
scatter(ones(1,length(stats1.percentMotif)), stats1.percentMotif, markerSize,  '.', 'MarkerEdgeColor', firstColor, 'MarkerFaceColor', firstColor, 'LineWidth', lineWidth, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha);
hold on
scatter(2*ones(1,length(stats2.percentMotif)), stats2.percentMotif, markerSize,  '.', 'MarkerEdgeColor', secondColor, 'MarkerFaceColor', secondColor, 'LineWidth', lineWidth, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha);
axis([0 6 0 1])
ylabel('Percentage surface that is blebby')
title([p.name1 ' (left) to ' p.name2 ])



function stats = collectDataMovie(cellsList, mainDirectory)

% initialize variables
stats.numMotifs = nan(1,length(cellsList));
stats.percentMotif = nan(1,length(cellsList));
stats.meanVolume = nan(1,length(cellsList));
stats.curvature = [];

% iterate over the list of cells
for c = 1:length(cellsList)
    
    % load motif and intensity statistics
    analysisPath = fullfile(mainDirectory, cellsList{c}, 'Morphology', 'Analysis', 'IntensityBlebCompare');
    statsStruct = load(fullfile(analysisPath, 'stats.mat'));
    comparePatches = statsStruct.comparePatches;
    compareFaces = statsStruct.compareFaces;
    convert = statsStruct.convert;
    
    % calculate statistics
    stats.numMotifs(1,c) = sum(comparePatches.isProtrusion);
    stats.percentMotif(1,c) = sum(compareFaces.isProtrusion)/length(compareFaces.isProtrusion);
    stats.meanVolume(1,c) = median(convert.volume*comparePatches.volume(logical(comparePatches.isProtrusion)));
    stats.curvature = [stats.curvature, convert.curvature*compareFaces.curvature];
    
end