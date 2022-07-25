imageList = [1,2,4:9,11:16];
allSept = [0];
meanCyto = [176,184,246,207,291,332,316,298,290,219,319,324,242,188]; %mean cytoplasmic intensities, determined by using FIJI to measure mean intensities of several large non-nuclear ROIs for all cells
for i = 1:length(imageList)
    disp(['-------Cell ' num2str(imageList(i)) '--------'])
    
    imagePath = ['/project/bioinformatics/Danuser_lab/melanoma/analysis/Weems/Menagerie/MV3/Sept6-GFP/Vitro_Gel/191119_' num2str(imageList(i)) '/Morphology/Analysis/Mesh/'];
    
%%Load Edge Image%%
    tiff_infoOUT = imfinfo([imagePath 'edgeImage_1_1.tif']); % return tiff structure, one element per image
    tiff_stackOUT = imread([imagePath 'edgeImage_1_1.tif'], 1) ; % read in first image
%%concatenate each successive tiff to tiff_stack
    for ii = 2 : size(tiff_infoOUT, 1)
    temp_tiffOUT = imread([imagePath 'edgeImage_1_1.tif'], ii);
    tiff_stackOUT = cat(3 , tiff_stackOUT, temp_tiffOUT);
    end
%%Remove Zeroes
    nonZ_tiffOUT = tiff_stackOUT(tiff_stackOUT~=0);

%%%Measure occupancy of voxels higher than meanCyto%%%
    abovecyto_count = nonZ_tiffOUT>meanCyto(i);
    abovecyto_tiffOUT = sum(abovecyto_count);
    Sept = abovecyto_tiffOUT/numel(nonZ_tiffOUT); %%percentage of voxels in edge above mean cytoplasmic value%%

%%%add to dataset of all metrics%%%
    temp_Sept = Sept;
    allSept = cat(1, allSept, temp_Sept);

end
