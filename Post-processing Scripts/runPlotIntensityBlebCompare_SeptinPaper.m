function runPlotIntensityBlebCompareSeptins2020()

% septins 2020
p.mainDirectory = '/archive/bioinformatics/Danuser_lab/melanoma/analysis/Weems/Menagerie/MV3/Sept6-GFP/Basal/2022reduxWeiner/';
p.savePath = '/archive/bioinformatics/Danuser_lab/melanoma/analysis/Weems/Menagerie/MV3/Sept6-GFP/Basal/2022reduxWeiner/figures5_23_22';




%% p.cellsList{1} = '181204/Cell3';
p.cellsList{1} = '181204/Cell4';
p.cellsList{2} = '181105/Cell1';
p.cellsList{3} = '181105/Cell7';
p.cellsList{4} = '181105/Cell8';
p.cellsList{5} = '200221/Cell2';
p.cellsList{6} = '200221/Cell7';
p.cellsList{7} = '200221/Cell8';
p.cellsList{8} = '200221/Cell9';
p.cellsList{9} = '171201/Cell5';
p.cellsList{10} = 'newCells/200218/Cell5';
p.cellsList{11} = 'newCells/200218/Cell6';
p.cellsList{12} = 'newCells/200218/Cell8';
p.cellsList{13} = 'newCells/200218/Cell9';
p.cellsList{14} = 'newCells/200218/Cell10';
p.cellsList{15} = 'newCells/190221/Cell8';
p.cellsList{16} = 'newCells/190221/Cell9';
p.cellsList{17} = '170922/Cell1';
p.analyzeDiffusion = 0;
p.analyzeForwardsMotion = 0;
p.analyzeVonMises = 1;
p.analyzeDistance = 1;

% make the plots
plotIntensityBlebCompare_SeptinPaper(p);
