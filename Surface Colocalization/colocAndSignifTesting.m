clear
allPvals = [];
allSpearmanCorr = [];

imageList = [1,2,3,4,6,7,8,9,10]; %%%Sept6 vs AktPH (aligned)%%%
    for z = 1:10

for a = 1:length(imageList)
      
     %%%%% Dataset: Sept6 vs AktPH %%%%%%
    surfacePath = ['/archive/bioinformatics/Danuser_lab/melanoma/analysis/Weems/hiRes3D/ASLM/MV3/Sept6-GFP aktPH-SNAP-TMRStar/190822/Cell' num2str(imageList(a)) ' (Sept6 1p5)/Morphology/Analysis/Mesh/surface_1_1.mat'];
    refPath = ['/archive/bioinformatics/Danuser_lab/melanoma/analysis/Weems/hiRes3D/ASLM/MV3/Sept6-GFP aktPH-SNAP-TMRStar/190822/Cell' num2str(imageList(a)) ' (Sept6 1p5)/Morphology/Analysis/Intensity/intensity_1_1.mat'];
    expPath = ['/archive/bioinformatics/Danuser_lab/melanoma/analysis/Weems/hiRes3D/ASLM/MV3/Sept6-GFP aktPH-SNAP-TMRStar/190822/Cell' num2str(imageList(a)) '/Morphology/Analysis/Intensity/intensity_1_1.mat'];
  
    
    refInten = load(refPath);
    expInten = load(expPath);
    surf = load(surfacePath);
    
                
            %%%%%% The following code downsamples the intensity data by selecting a given number of faces equidistant from each other on the mesh %%%%%%
            
                %This script generates the faces index of a mesh, randomly distributed on
                %the surface. numDots is a parameter to give the numebr of seeds
                %(downsampling parameter)

                %load the surface 
                if ~exist('surf','var')
                error('load a surface');
                end 
                %choose the number of seeds for downsampling data on a mesh
                numDots=200; 

                %find the position of each face (if intensity is face-based)
                numFaces=size(surf.surface.faces,1); 
                facePositions = zeros(numFaces,3);
                for f = 1:numFaces
                    verticesFace = surf.surface.faces(f,:);
                    facePositions(f,:) = (surf.surface.vertices(verticesFace(1),:) + surf.surface.vertices(verticesFace(2),:) + surf.surface.vertices(verticesFace(3),:))/3;
                end

                %find the center of seeds randomly distributed on a surface
                polkaDotSeedsFaces = fps_euclidean(facePositions,numDots); 
    
                refSelectedFaces = refInten.faceIntensities.mean(polkaDotSeedsFaces,:);
                expSelectedFaces = expInten.faceIntensities.mean(polkaDotSeedsFaces,:);
   
    testCorrResults = [];
    
    parfor c = 1:1000
        %%%%%% Shuffles experimental intensity dataset %%%%%%
        rows = size(expSelectedFaces,1);
        perm = randperm(rows);
        shuffInten = expSelectedFaces(perm,:);

        %%%%%% Determines correlation coefficient b/w reference dataset and shuffled experimental dataset %%%%%%
        testCorrImmediate = corr(refSelectedFaces, shuffInten,'Type','Spearman');
        testCorrResults = [testCorrResults, testCorrImmediate,];
    end
    
%%%%%% Determines P-value (not to be confused with p-value) by measuring the fraction of testCorrResults above the correlation b/w reference dataset and UNSHUFFLED experimental dataset %%%%%%    
    realCorr = corr(refSelectedFaces, expSelectedFaces,'Type','Spearman');
    Pval = sum(testCorrResults>realCorr)/1000;
    allPvals = [allPvals,Pval];
%%%%%% Determines the Spearman Correlation b/w the 2 observed distributions %%%%%%  
    spearmanCorr = corr(refInten.faceIntensities.mean, expInten.faceIntensities.mean,'Type','Spearman');
    allSpearmanCorr = [allSpearmanCorr, spearmanCorr];
end

end

clearvars -except allPvals allSpearmanCorr;