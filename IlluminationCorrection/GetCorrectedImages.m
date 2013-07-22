function OutIm = GetCorrectedImages( A, T, poly, FlAv, AFAv, BGAv )

%Visualize the meshes
[ X,Y ] = meshgrid(1:1:2048,1:1:2048);
Z = [ X(:).^4 Y(:).^4 X(:).*Y(:).^3 X(:).^3.*Y(:) X(:).^2.*Y(:).^2  ... 4th
      X(:).^3 Y(:).^3 X(:).*Y(:).^2 X(:).^2.*Y(:)                   ... 3rd
      X(:).^2 Y(:).^2 X(:).*Y(:)                                    ... 2nd
      X(:) Y(:) ];                                                  ... 1st
ZAF = Z-repmat(poly(2).meanInd,[size(Z,1) 1]);
ZAF = ZAF./repmat(poly(2).normInd,[size(ZAF,1) 1]);
% ZAF = exp((ZAF*WAF).*repmat(normAFVal,[length(ZAF) 1]));
ZAF = ZAF*poly(2).coeffs';%.*normAFVal+repmat(meanAFVal,[length(ZAF) 1]);
ZAF = reshape(ZAF,size(X,1),size(X,2));

ZFl = Z-repmat(poly(1).meanInd,[size(Z,1) 1]);
ZFl = ZFl./repmat(poly(1).normInd,[size(ZFl,1) 1]);
% ZFl = exp((ZFl*WFl).*repmat(normFlVal,[length(ZFl) 1]));
ZFl = ZFl*poly(1).coeffs';%.*normFlVal+repmat(meanFlVal,[length(ZFl) 1]);
ZFl = reshape(ZFl,size(X,1),size(X,2));
ZFl = (ZFl+ZAF)./2;

ZBG = Z-repmat(poly(3).meanInd,[size(Z,1) 1]);
ZBG = ZBG./repmat(poly(3).normInd,[size(ZBG,1) 1]);
% ZBG = exp((ZBG*WBG).*repmat(normBGVal,[length(ZBG) 1]));
ZBG = ZBG*poly(3).coeffs';%.*normBGVal+repmat(meanBGVal,[length(ZBG) 1]);
ZBG = reshape(ZBG,size(X,1),size(X,2));

%Figure out the places in the image where the max and the min for the polys
%are located and get corresponding image values
FlAv = log(FlAv);
AFAv = log(AFAv);
BGAv = log(BGAv);

[rowAFMax,colAFMax] = GetDialatedMinMaxIndices( ZAF, 1 );
AFMax = GetTenPCMinMaxVals( rowAFMax, colAFMax, AFAv, 1 );
[rowAFMin,colAFMin] = GetDialatedMinMaxIndices( ZAF, 0 );
AFMin = GetTenPCMinMaxVals( rowAFMin, colAFMin, AFAv, 0 );

[rowFlMax,colFlMax] = GetDialatedMinMaxIndices( ZFl, 1 );
FlMax = GetTenPCMinMaxVals( rowFlMax, colFlMax, FlAv, 1 );
[rowFlMin,colFlMin] = GetDialatedMinMaxIndices( ZFl, 0 );
FlMin = GetTenPCMinMaxVals( rowFlMin, colFlMin, FlAv, 0 );

[rowBGMax,colBGMax] = GetDialatedMinMaxIndices( ZBG, 1 );
BGMax = GetTenPCMinMaxVals( rowBGMax, colBGMax, BGAv, 1 );
% [rowBGMin,colBGMin] = GetDialatedMinMaxIndices( ZBG, 0 );
% BGMin = GetTenPCMinMaxVals( rowBGMin, colBGMin, BGAv, 0 )
%Very noisy use AF range
BGMin = BGMax-(AFMax-AFMin);

FlSubSurf = GetRescaledSurface( ZFl, FlMin, FlMax );
AFSubSurf = GetRescaledSurface( ZAF, AFMin, AFMax );
BGSubSurf = GetRescaledSurface( ZBG, BGMin, BGMax );

figure,imagesc(AFSubSurf);colorbar
figure,imagesc(FlSubSurf);colorbar
figure,imagesc(BGSubSurf);colorbar

FlSubSurf = repmat(FlSubSurf,[1 1 size(A,3)]);
AFSubSurf = repmat(AFSubSurf,[1 1 size(A,3)]);
BGSubSurf = repmat(BGSubSurf,[1 1 size(A,3)]);
FlSubSurf(T~=2) = 0;
AFSubSurf(T~=1) = 0;
BGSubSurf(T~=0) = 0;
FlSubSurf = AFSubSurf+FlSubSurf+BGSubSurf;
clear AFSubSurf BGSubSurf

OutIm = exp(log(cast(A,'double'))-FlSubSurf);
OutIm = cast(OutIm,'uint16');

end