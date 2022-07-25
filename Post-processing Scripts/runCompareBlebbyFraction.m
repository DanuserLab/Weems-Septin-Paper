function runCompareCellStatsAndrewRevise()

% runCompareCellStats - runs compareCellStats.m

% set the save directory
p.savePath = '/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/mdriscoll/Septin/revise/figures';


% Negative control, AktPH-GFP
p.name = 'AktPH-GFP';
p.mainDirectory = '/project/bioinformatics/Danuser_lab/melanoma/analysis/Weems/hiRes3D/ASLM/MV3, EtOH ON/AktPH-GFP-GFP/200305/';
p.cellsList{1} = '/Cell1';
p.cellsList{2} = '/Cell3';
p.cellsList{3} = '/Cell4';
p.cellsList{4} = '/Cell6';
p.cellsList{5} = '/Cell7';
p.cellsList{6} = '/Cell10';
p.cellsList{7} = '/Cell12';
p.cellsList{8} = '/Cell13';
p.cellsList{9} = '/Cell14';
p.cellsList{10} = '/Cell16';
p.cellsList{11} = '/Cell17';
p.cellsList{12} = '/Cell18';
p.cellsList1 = p.cellsList;
p.mainDirectory1 = p.mainDirectory;
p.name1 = p.name;

% FCF AktPH-GFP
p.name = 'FCF AktPH-GFP';
p.mainDirectory = '/project/bioinformatics/Danuser_lab/melanoma/analysis/Weems/hiRes3D/ASLM/MV3, FCF ON/AktPH-GFP-GFP/200305/';
p.cellsList{1} = '/Cell1';
p.cellsList{2} = '/Cell5';
p.cellsList{3} = '/Cell6';
p.cellsList{4} = '/Cell7';
p.cellsList{5} = '/Cell8';
p.cellsList{6} = '/Cell9';
p.cellsList{7} = '/Cell10';
p.cellsList{8} = '/Cell12';
p.cellsList{9} = '/Cell13';
p.cellsList{10} = '/Cell14';
p.cellsList{11} = '/Cell15';
p.cellsList{12} = '/Cell17';
p.cellsList{13} = '/Cell20';
p.cellsList{14} = '/Cell23';
p.cellsList{15} = '/Cell24';
p.cellsList2 = p.cellsList;
p.mainDirectory2 = p.mainDirectory;
p.name2 = p.name;


compareBlebbyFraction(p);