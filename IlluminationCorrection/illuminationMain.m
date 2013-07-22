tilesPerRow = 42;
reverseEvenRowsForWriting = 0;
%Add path to compiled SPAMS library
SpamsFolder = 'E:/SharedFolder/ARBs_MMI-31_28R_p4067/spams-matlab/build/';
addpath(SpamsFolder);
StatsFolder = 'E:/SharedFolder/ARBs_MMI-31_28R_p4067/stats/';
addpath(StatsFolder);
clear SpamsFolder

fileFolder = 'E:/SharedFolder/ARBs_MMI-31_28R_p4060/C2/';
fileNames1  = 'ARBs_MMI-31_28R_p4060C2*.tif';
fileName2  = 'out1.tif';
AFF = ReadImageSeries( fileFolder, fileNames1 );
TFF = ReadImageStack( fileFolder, fileName2 );
for i=1:size(AFF,3)
    AFF(:,:,i) = medfilt2( AFF(:,:,i), [5 5] );
end

tic
numTiles = size(AFF,3);
lin  = reshape( 1:numTiles, tilesPerRow, numTiles/tilesPerRow );
FlIms = zeros(size(AFF,1),size(AFF,2),2);
AFIms = zeros(size(AFF,1),size(AFF,2),2);
BGIms = zeros(size(AFF,1),size(AFF,2),2);
for j = 1:2
    lin1 = lin(:,j:2:end);
    lin1 = lin1(:);
    A = AFF(:,:,lin1);
    T = TFF(:,:,lin1);
    [ FlIms(:,:,j), AFIms(:,:,j), BGIms(:,:,j) ] = GetMeanImages( A, T );
end
clear A T
toc

tic
Poly = [];
for j=1:2
    poly = GetPolynomials( FlIms(:,:,j), AFIms(:,:,j), BGIms(:,:,j) );
    Poly = [Poly; poly];
end
toc

tic
OutIm = zeros(size(AFF),'uint16');
for j = 1:2
    lin1 = lin(:,j:2:end);
    lin1 = lin1(:);
    A = AFF(:,:,lin1);
    T = TFF(:,:,lin1);
    OutIm(:,:,lin1) = GetCorrectedImages( A, T, Poly(j,:), FlIms(:,:,j),...
                                              AFIms(:,:,j), BGIms(:,:,j) );
end
toc

if reverseEvenRowsForWriting
    OutImTmp = zeros(size(AFF),'uint16');
    lin1 = lin(:,1:2:end);
    lin1 = lin1(:);
    OutImTmp(:,:,lin1) = OutIm(:,:,lin1);
    lin1 = lin(:,2:2:end);
    lin2 = flipud(lin1);
    lin1 = lin1(:); lin2 = lin2(:);
    OutImTmp(:,:,lin2) = OutIm(:,:,lin1);
    OutIm = OutImTmp;
    clear OutImTmp lin2
end

tic
clear A T AFF TFF lin lin1 meana numTiles poly BGIms...
      tilesPerRow AFIms FlIms Poly StatsFolder fileName2 fileNames1 i j
WriteImageSeries( OutIm, fileFolder, 'corrected', 'tif' );
toc

% % Hold license for breaks aka I hate matlab
% i=1
% while i~=10
%     meana = rand(100,100); meana  = medfilt2(meana,[3 3]);
% %     meana = nanmean(meana);
% end