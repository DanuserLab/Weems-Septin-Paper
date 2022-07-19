%This script generates the faces index of a mesh, randomly distributed on
%the surface. numDots is a parameter to give the numebr of seeds
%(downsampling parameter)

%load the surface 
if ~exist('surface','var')
error('load a surface')
end 
%choose the number of seeds for downsampling data on a mesh
numDots=10; 

%find the position of each face (if intensity is face-based)
numFaces=size(surface.faces,1); 
facePositions = zeros(numFaces,3);
for f = 1:numFaces
    verticesFace = surface.faces(f,:);
    facePositions(f,:) = (surface.vertices(verticesFace(1),:) + surface.vertices(verticesFace(2),:) + surface.vertices(verticesFace(3),:))/3;
end

%find the center of seeds randomly distributed on a surface
polkaDotSeeds_faces = fps_euclidean(facePositions,numDots); 
