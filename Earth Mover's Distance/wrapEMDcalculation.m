function wrapEMDcalculation
fileDir='/project/bioinformatics/Danuser_lab/melanoma/analysis/Hanieh/AndrewData/demo';
sheet='EMD-FCFOn';
%load the cytosol background noise - it is measured by Andrew manually for
%each cell and saved in an excel file
[NumData, StrData]=xlsread(fullfile(['/project/bioinformatics/Danuser_lab/melanoma/analysis/Hanieh/U01Grant/EMDPlot'], 'EMD-MV3Cells.xlsx'),sheet);
CellList=NumData(:,1);
MeanCytoplasmic=NumData(:,6)/(2^16); % 2^16 to convert to bit

%calculate the EMD for a single cell
c=1; % index of a cell in the list
cellID=CellList(c);

%load the mesh for a chosen cell
filename='surface_1_1.mat';
meshPath = [fileDir '/Cell' num2str(cellID) '/Morphology/Analysis/Mesh'];
load(fullfile(meshPath,filename))
%load the intensity
filename='intensity_1_1.mat';
intensityPath = [fileDir '/Cell' num2str(cellID) '/Morphology/Analysis/Intensity'];
load(fullfile(intensityPath,filename))

%measure the EMD
% Compute LB eigenstuff using the FEM
nEigs = 100;
FEM = firstOrderFEM(surface.vertices,surface.faces); % set up FEM matrices
[evecs,evals] = eigs(FEM.laplacian,FEM.vtxInnerProds,nEigs,'sm'); %original test was 300
evals = diag(evals);

%1st distribution: normalized molecular distribution
rho0=faceIntensities.mean/MeanCytoplasmic(c); % normalized for the ratio
surfaceBackground=median(faceIntensities.mean);
ThresholdRatio=surfaceBackground/MeanCytoplasmic(c);
%normalize the intensity on the surface
rho0=rho0 - ThresholdRatio; %subtract surface bg
Ind=    find(rho0 < 0);%subtract surface bg
rho0(Ind)=0;%subtract surface bg

%2nd distribution: uniform distribution of molecular distribution
rho1=sum(rho0)/length(rho0)*ones(size(rho0));

% Compute EMD
structure = precomputeEarthMoversADMM(surface.vertices, surface.faces, evecs(:,2:end));
[EMD_dist] = earthMoversADMM(surface.vertices, surface.faces, rho0, rho1, structure);

% save EMD in the cellPath
savePath = [fileDir '/Cell' num2str(cellID) '/Morphology/Analysis/EMD']
if ~isdir(savePath) mkdir(savePath); end
filename='earthMoverDist_1_1.mat';
save(fullfile(savePath,filename),'EMD_dist')

end