%function poly = GetPolynomials(AF, TF, tilesPerRow)
numTiles = size(AF,3);
lin  = reshape( 1:numTiles, tilesPerRow, numTiles/tilesPerRow );
for j = 1:2
lin1 = lin(:,j:2:end);
lin1 = lin1(:);
A = AF(:,:,lin1);
T = TF(:,:,lin1);
tic
if size(A,3) ~= size(T,3)
    error('Needs same number of slices in labels and gray images');
end
%Get the indices and values from each tile
numCurrentTiles = length(lin1);
countFl    = zeros(numCurrentTiles,1);
countAF    = zeros(numCurrentTiles,1);
AFIndValsC = cell(numCurrentTiles,1);
FlIndValsC = cell(numCurrentTiles,1);
for i=1:numCurrentTiles
    curIm         = A(:,:,lin1(i));
    curLab        = T(:,:,lin1(i));
    currentAF     = curIm(curLab==1);
    [xAF,yAF]     = find(curLab==1);
    AFIndValsC{i} = [ currentAF xAF yAF ];
    countAF(i)    = length(yAF);
    currentFl     = curIm(curLab==2);
    [xFl,yFl]     = find(curLab==2);
    FlIndValsC{i} = [ currentFl xFl yFl ];
    countFl(i)    = length(yFl);
    clear curIm curLab currentAF currentFl xFl yFl xAF yAF
end
toc
tic
%Put the above into arrays
%Using third order exponential
AFInd  = zeros( sum(countAF), 9 );
AFVal  = zeros( sum(countAF), 1 );
FlInd  = zeros( sum(countFl), 9);
FlVal  = zeros( sum(countFl), 1 );
sumFl = 1;
sumAF = 1;
for i=1:numCurrentTiles
  if countFl(i)~=0
    x = FlIndValsC{i}(:,2); y = FlIndValsC{i}(:,3);
    FlInd(sumFl:sumFl+countFl(i)-1,:) = [ x.^3 y.^3 x.*y.^2 x.^2.*y ... 3rd
                                          x.^2 y.^2 x.*y ...            2nd
                                          x y ];      %ones(size(x,1),1)1st
    FlVal(sumFl:sumFl+countFl(i)-1,:) = log(cast(FlIndValsC{i}(:,1),...
                                                                'double'));
  end
  if countAF(i)~=0
    x = AFIndValsC{i}(:,2); y = AFIndValsC{i}(:,3);
    AFInd(sumAF:sumAF+countAF(i)-1,:) = [ x.^3 y.^3 x.*y.^2 x.^2.*y ... 3rd
                                          x.^2 y.^2 x.*y ...            2nd
                                          x y ];      %ones(size(x,1),1)1st
    AFVal(sumAF:sumAF+countAF(i)-1,:) = log(cast(AFIndValsC{i}(:,1),...
                                                                'double'));
  end
  sumFl = sumFl+countFl(i);
  sumAF = sumAF+countAF(i);
end
clear count sumFl sumAF countFl countAF FlIndValsC AFIndValsC x y i 
toc
param.numThreads=-1; % all cores (-1 by default)
param.verbose=true;   % verbosity, false by default
param.lambda=0.05; % regularization parameter
param.it0=10;      % frequency for duality gap computations
param.max_it=200; % maximum number of iterations
param.L0=0.1;
param.tol=1e-3;
param.intercept=true;
param.pos=false;
param.compute_gram=true;
fprintf('\nFISTA + Group Lasso\n');
param.loss='square';
param.regul='group-lasso-l2';
param.size_group=5;
W0=zeros(size(FlInd,2),size(FlVal,2));
%Normalize inputs
FlInd1=FlInd-repmat(mean(FlInd),[size(FlInd,1) 1]);
FlInd1=mexNormalize(FlInd1);
FlVal1=FlVal-repmat(mean(FlVal),[size(FlVal,1) 1]);
FlVal1=mexNormalize(FlVal1);
tic
[WFl optim_infoFl]=mexFistaFlat(FlVal1,FlInd1,W0,param);
t=toc;
%For Autoflour
AFInd=AFInd-repmat(mean(AFInd),[size(AFInd,1) 1]);
% AFInd=mexNormalize(AFInd);
AFVal=AFVal-repmat(mean(AFVal),[size(AFVal,1) 1]);
AFVal=mexNormalize(AFVal);
tic
[WAF optim_infoAF]=mexFistaFlat(AFVal,AFInd,W0,param);
t=toc;


end
clear lin1 A T
%
%end