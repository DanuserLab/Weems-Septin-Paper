function structure = precomputeEarthMoversADMM(X, T, curlFunctionBasis, harmonicFields)
% Precomputes structures that can be used for multiple computations of EMD
% on a mesh independently of the distribution.  Takes three arguments and
% one fourth optional argument.  The first two (X,T) determine a triangle
% mesh.  The third -- curlFunctionBasis -- is a set of per-vertex functions
% on the mesh whose gradients rotated by 90 degrees in the tangent plane
% will determine the curl component of J.  The optional argument is a basis
% of harmonic vector fields (per-face).

% rho1 and rho2 should be per-vertex and should integrate to 1.
% Can use two different functions for curl and div part -- this avoids
% boundary condition issues.

nf = size(T,1);
nv = size(X,1);

structure.FEM = firstOrderFEM(X, T);
structure.curlFunctionBasis = curlFunctionBasis;

structure.nHarmonic = 0;
if nargin == 4
    structure.harmonicFields = harmonicFields;
    structure.nHarmonic = length(harmonicFields);
end

ncurl = size(curlFunctionBasis,2);

if ncurl < 500
    structure.B = zeros(nf*3,ncurl); % dense faster for small matrices
else
    structure.B = sparse(nf*3,ncurl);
end

for i=1:3
    structure.B(i:3:end,:) = structure.FEM.faceInnerProds*structure.FEM.grad{i}*curlFunctionBasis;
end

if nargin == 4
    fprintf('Adding harmonic vector field basis...\n');
    for i=1:length(harmonicFields)
        rotBasis = cross(mesh.faceNormals, harmonicFields{i});
        newIdx = size(B,2)+1;
        for j=1:3
            innerProdMtx = structure.FEM.faceInnerProds*rotBasis(:,j);
            structure.B(j:3:end,newIdx) = innerProdMtx;
        end
    end
end

fprintf('Computing pseudoinverse...\n');
structure.leastSquaresMatrix = pinv(structure.B'*structure.B);
fprintf('Done computing pseudoinverse.\n');

%% Prefactor Laplacian via sparse cholesky

W = structure.FEM.laplacian; % check sign
AM = sum(structure.FEM.vtxInnerProds,2); % area weights
nopts = nv;

W = -W;

W(1, :) = zeros(1,nopts);
W(:, 1) = zeros(nopts,1);
W(1, 1) = 1;

s = symamd(W);
W=W(s,s);
R = chol(W);

structure.W = W;
structure.R = R;
structure.s = s;

%% Compute face normals

normalf = cross( X(T(:,2),:)'-X(T(:,1),:)', ...
                 X(T(:,3),:)'-X(T(:,1),:)' );
d = sqrt( sum(normalf.^2,1) ); d(d<eps)=1;
structure.faceNormals = (normalf ./ repmat( d, 3,1 ))';