function [distance,J,gradPart,curlPart] = earthMoversADMM(X, T, rho1, rho2, structure)
% rho1 and rho2 should be per-vertex and should sum to 1 (assume they're 
% already scaled by area weights).

n = size(X,1);
nf = size(T,1);

FEM = structure.FEM;

% We want to solve
%    min \int_M ||J||_2
%    s.t. div (J) = rho_1 - rho_0
% To do so, write J = sum_i [ a_i grad phi_i + b_i Rgrad phi_i ].  Then, we
% know that div(J) = laplacian(a), so we want laplacian(a)=rho_1-rho_0.
% This gives a in closed form:

RHS = (rho2-rho1);

batchsize = 500;
soln = zeros(n,size(rho1,2));
for i=1:batchsize:size(rho1,2)
    fprintf('Solve %d of %d...\n',i,size(rho1,2));
    idx2 = min(i+batchsize-1,size(rho1,2));
    
    curRhs = RHS(structure.s,i:idx2);
    curRhs(structure.s(1),:) = 0;
    
    soln(:,i:idx2) = -structure.R\(structure.R'\curRhs);
end
soln(structure.s,:) = soln;

integral = sum(FEM.vtxInnerProds*soln);
totalArea = sum(sum(FEM.vtxInnerProds));
a = bsxfun(@minus,soln,integral/totalArea);

% So now we just have to find b, so we'll set the grad part in stone.
nDists = size(rho1,2);
gradPart = cell(nDists,1);
for j=1:nDists
    for i=1:3
        gradPart{j}(:,i) = FEM.grad{i} * a(:,j);
    end
end

%% Find b

% Rather than minimize ||grad(a) + R grad(b)||_2 per face, we're going to
% apply R to grad(a) instead to get an alternative optimization.  This is
% computationally a bit easier because rotations are annoying to compute as
% matrix operators.

for j=1:nDists
    % Rotate grad part in the tangent plane
    rotGradPart = cross(structure.faceNormals, gradPart{j}); % check sign
    bb = rotGradPart';
    bb = bsxfun(@times,FEM.faceAreas',bb);
    b(:,j) = bb(:); % [r1x;r1y;r1z;r2x;r2y;r2z;...]
end

ncurl = size(structure.curlFunctionBasis,2);
B = structure.B;

%% Do ADMM

nharmonic = structure.nHarmonic;
c = zeros(ncurl+nharmonic,nDists);
y = zeros(3*nf,nDists);

% Fancier implementations would autotune the ADMM "rubber band" for
% enforcing constraints.  See Boyd's survey paper for some simple
% strategies for tuning this value.  ADMM converges regardless.
rho = 10000;

w = b; % rotated grad part

Bc = zeros(size(B,1),nDists); % since we initialized c=0
pp = structure.leastSquaresMatrix*B';

% See the supplemental documentation from the SIGGRAPH paper for derivation
% of the iterative ADMM formulas below. 

for i=1:100 % usually stops far before this point
    fprintf(1,'ADMM iteration %d...\n',i);
    
    % J step
    z = Bc + w - y/rho;
    rhoModZ = rho*sqrt(z(1:3:end,:).^2+z(2:3:end,:).^2+z(3:3:end,:).^2);
    aa = 1-4./rhoModZ;
    aa(rhoModZ <= 4) = 0;
    aaa = zeros(3*size(aa,1),size(aa,2));
    aaa(1:3:end,:) = aa;
    aaa(2:3:end,:) = aa;
    aaa(3:3:end,:) = aa;
    J = aaa.*z;
    
    % c step
    oldC = c;
    c = pp*(y/rho + J - w);
    
    % Convergence criterion
    if sum((c(:)-oldC(:)).^2/nDists) < 1e-5 && i > 3
        break
    end
    
    Bc = B*c;
    
    % y step
    y = y + rho*(J-Bc-w); 
end

% Post-ADMM endgame -- compute J and read off EMD

curX = c;
curlPart = cell(nDists,1);
J = cell(nDists,1);
for j=1:nDists
    for i=1:3
        curlPart{j}(:,i) = FEM.grad{i} * (structure.curlFunctionBasis * curX(1:ncurl,j)); % sign?
    end
    curlPart{j} = cross(curlPart{j},structure.faceNormals);
    J{j} = curlPart{j} + gradPart{j};
end

resid = B*curX + b;
residNorms = sqrt(resid(1:3:end,:).^2 + resid(2:3:end,:).^2 + resid(3:3:end,:).^2);
distance = sum(residNorms,1)';