% function fltrdSurf = GetFilteredSurface(A, T)
tic
numTiles = size(A,3);
fltrdSurf = zeros(size(A));
parfor i=1:numTiles
	fltrdSurf(:,:,i) = cast(medfilt2(A(:,:,i),[5 5]),'double');
    %fltrdSurf(:,:,i) = mat2gray(medfilt2(A(:,:,i),[5 5]));
end
sortedNoThresh = sort(fltrdSurf,3,'descend');
%Take the top 10%
sortedNoThresh = cast(sortedNoThresh(:,:,1:ceil(numTiles*0.1)),'double');
toc
tic
sortedNoThresh = medfilt3( sortedNoThresh, [ 5 5 3 ] );
toc
meanInNoThresh = mean(sortedNoThresh,3);
stdInNoThresh  = std(sortedNoThresh,0,3);
figure,imagesc(stdInNoThresh); colorbar;
figure,imagesc(meanInNoThresh); colorbar;
figure,imagesc(log(stdInNoThresh)); colorbar;
figure,imagesc(log(meanInNoThresh)); colorbar;
for i=1:ceil(numTiles*0.1)
    figure(10),imagesc(log(sortedNoThresh(:,:,i))); colorbar;
    pause(0.25)
end

%Estimating the FG even in the top 10%
threshSortedThresh = graythresh(mat2gray(sortedNoThresh));
threshSortedThresh = threshSortedThresh*max(sortedNoThresh(:));
sortedThresh = sortedNoThresh;
sortedThresh(sortedThresh<threshSortedThresh) = 0;
count = sum(sortedThresh~=0,3);
meanSortedThresh = sum(sortedThresh,3);
meanSortedThresh = meanSortedThresh./count;
figure,imagesc(log(meanSortedThresh)); colorbar;



fltrdSurfThresh = fltrdSurf;
fltrdSurfThresh(T==0) = 0;
sortedThresh = sort(fltrdSurfThresh,3,'descend');
%Take the top 10%
sortedThresh = cast(sortedThresh(:,:,1:ceil(numTiles*0.1)),'double');
threshSortedThresh = graythresh(mat2gray(sortedThresh));
threshSortedThresh = threshSortedThresh*max(sortedThresh(:));
sortedThresh(sortedThresh<threshSortedThresh) = 0;
count = sum(sortedThresh~=0,3);
meanSortedThresh = sum(sortedThresh,3);
meanSortedThresh = meanSortedThresh./count;
% meanInThresh = mean(sortedThresh,3);
% stdInThresh  = std(sortedThresh,0,3);
% meannThresh  = medfilt2(meanInA,[20 20]);
% myfilter = fspecial('gaussian',[3 3], 6);
% myfilteredimage = imfilter(meannA, myfilter, 'replicate');
% figure,imagesc(stdInThresh)
figure,imagesc(meanSortedThresh);colorbar;
figure,imagesc(log(medfilt2(meanSortedThresh,[11 11])));colorbar;
% figure,imagesc(log(stdInThresh))
for i=1:ceil(numTiles*0.1)
    figure(11),imshow(mat2gray(sortedThresh(:,:,i)));% colorbar;
    pause(0.25)
end
medfilt3
%%%%%%%Old Code%%%%%%%%%%%
% % % % % medKernSz = 9;
% % % % % windowSz  = 128;
% % % % % 
% % % % % numTiles = size(A,3);
% % % % % %Median filter to remove thermal noise
% % % % % fltrdSurf = zeros(size(A));
% % % % % parfor i=1:numTiles
% % % % %     fltrdSurf(:,:,i) = mat2gray(medfilt2(A(:,:,i),[3 3]));
% % % % % end
% % % % % 
% % % % % %Extract contours of the images
% % % % % contourIm = zeros(size(A));
% % % % % parfor i=1:numTiles
% % % % %     %Median or mean filter with a large kernel
% % % % %     contourIm(:,:,i) = medfilt2(fltrdSurf(:,:,i),[medKernSz medKernSz]);
% % % % %     contourIm(:,:,i) = mat2gray(contourIm(:,:,i)-fltrdSurf(:,:,i));
% % % % %     contourIm(:,:,i) = im2bw(contourIm(:,:,i));
% % % % % end
% % % % % 
% % % % % %Compute greythresh levels for each tile and then compute lower bound
% % % % % greyThreshLevels = zeros( 1, numTiles );
% % % % % parfor i=1:numTiles
% % % % %     greyThreshLevels(i) = graythresh( fltrdSurf(:,:,i) );
% % % % % end
% % % % % 
% % % % % greyThreshLevelAll = graythresh( greyThreshLevels );
% % % % % greyThreshLevelLb  = greyThreshLevelAll - std( greyThreshLevels )
% % % % % 
% % % % % %Windowed greythresh to get the foreground regions
% % % % % startRow = 1:windowSz:size(A,1);
% % % % % endRow   = windowSz:windowSz:size(A,1);
% % % % % if endRow(end)~=size(A,1)
% % % % %     endRow = [endRow size(A,1)];
% % % % % end
% % % % % startCol = 1:windowSz:size(A,2);
% % % % % endCol   = windowSz:windowSz:size(A,2);
% % % % % if endCol(end)~=size(A,2)
% % % % %     endCol = [endCol size(A,2)];
% % % % % end
% % % % % 
% % % % % parfor n=1:numTiles
% % % % %     currentTile = fltrdSurf(:,:,n);
% % % % %     AA = currentTile(startRow:endRow,startCol:endCol);
% % % % %     curThresh = graythresh( AA );
% % % % %     if curThresh < greyThreshLevelLb
% % % % %         curThresh = greyThreshLevelAll;
% % % % %     end
% % % % %     AA(AA<curThresh) = 0;
% % % % % 	currentTile(startRow:endRow,startCol:endCol) = AA;
% % % % %     fltrdSurf(:,:,n) = currentTile;
% % % % % end
% % % % % 
% % % % % %Show thresholded images
% % % % % for n=1:numTiles:10
% % % % %     figure(10)
% % % % %     subplot(1,2,1),imshow(fltrdSurf(:,:,n));
% % % % %     subplot(1,2,2),imshow(mat2gray(A(:,:,n)));
% % % % %     pause(0.25)
% % % % % end
% % % % % for n = 1:numTiles
% % % % %     
% % % % % end
% end
