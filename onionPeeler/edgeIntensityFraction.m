all_edge_fractionCOUNT= [0];
all_edge_fractionINT = [0];
all_percent_changeCOUNTtoINT = [0];

imageList = [1,3:12,14:21,23];
for i = 1:length(imageList)
    disp(['-------Cell ' num2str(imageList(i)) '--------'])
    
    imagePath = ['/project/bioinformatics/Danuser_lab/melanoma/analysis/Weems/hiRes3D(cheap knockoff)/ASLM/MV3/NRAS-GFP/201009_GoodSurfaces_edge8/Cell' num2str(imageList(i)) '/Morphology/Analysis/Mesh/'];

    
%%%Load Interior Image%%%

    tiff_infoIN = imfinfo([imagePath 'interiorImage_1_1.tif']); % return tiff structure, one element per image
    tiff_stackIN = imread([imagePath 'interiorImage_1_1.tif'], 1) ; % read in first image
    
    %%%concatenate each successive tiff to tiff_stack
    for ii = 2 : size(tiff_infoIN, 1)
    temp_tiffIN = imread([imagePath 'interiorImage_1_1.tif'], ii);
    tiff_stackIN = cat(3 , tiff_stackIN, temp_tiffIN);
    end
    
    %%%Remove Zeroes
    nonZ_tiffIN = tiff_stackIN(tiff_stackIN~=0);
    
    %%%Sort and count voxels
    sorted_tiffIN = sort(nonZ_tiffIN, 'descend'); %sort descending
    count_IN = numel(sorted_tiffIN); %count rows 
   
    
%%%Load Edge Image%%%
    
    tiff_infoOUT = imfinfo([imagePath 'edgeImage_1_1.tif']); % return tiff structure, one element per image
    tiff_stackOUT = imread([imagePath 'edgeImage_1_1.tif'], 1) ; % read in first image
    
    %%%concatenate each successive tiff to tiff_stack
    for ii = 2 : size(tiff_infoOUT, 1)
    temp_tiffOUT = imread([imagePath 'edgeImage_1_1.tif'], ii);
    tiff_stackOUT = cat(3 , tiff_stackOUT, temp_tiffOUT);
    end
    
    %%%Remove Zeroes
    nonZ_tiffOUT = tiff_stackOUT(tiff_stackOUT~=0);
    
    %%%Sort and count voxels
    sorted_tiffOUT = sort(nonZ_tiffOUT, 'descend'); %sort descending
    count_OUT = numel(sorted_tiffOUT); %count rows 

    
%%%Analyze Images%%%
    
    %%%Calculate fraction of total voxels near the edge
    edge_fractionCOUNT = (count_OUT/(count_IN + count_OUT));
    
    %%%Calculate fraction of total intensity near the edge
    sum_nonZ_tiffOUT = sum(nonZ_tiffOUT);
    sum_nonZ_tiffIN = sum(nonZ_tiffIN);
    edge_fractionINT = (sum_nonZ_tiffOUT/(sum_nonZ_tiffIN + sum_nonZ_tiffOUT));
    
    
%%%Record Results%%%

    temp_all_edge_fractionINT = edge_fractionINT;
    all_edge_fractionINT = cat(1, all_edge_fractionINT, temp_all_edge_fractionINT);

    temp_all_edge_fractionCOUNT = edge_fractionCOUNT;
    all_edge_fractionCOUNT = cat(1, all_edge_fractionCOUNT, temp_all_edge_fractionCOUNT);
    
    temp_all_percent_changeCOUNTtoINT = ((edge_fractionINT-edge_fractionCOUNT)/edge_fractionCOUNT*100);
    all_percent_changeCOUNTtoINT = cat(1, all_percent_changeCOUNTtoINT, temp_all_percent_changeCOUNTtoINT);
end
