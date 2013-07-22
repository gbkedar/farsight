function [ FlValsIm, AFValsIm, BGValsIm ] = GetMeanImages( A, T )

if size(A,3) ~= size(T,3)
    error('Needs same number of slices in labels and gray images');
end

%Separate the fluorescence and the auto-fluorescence
FlVals       = zeros(size(A));
AFVals       = zeros(size(A));
FlVals(T==2) = A(T==2);
AFVals(T==1) = A(T==1);

A = sort(A,3,'ascend');
BGValsIm5 = zeros( size(A,1), size(A,2), 10 );
for i=1:10
    BGValsIm5(:,:,i) = medfilt2( A(:,:,i), [10 10] );
end
BGValsIm = mean(BGValsIm5,3);
BGValsIm(BGValsIm<1) = 1;
% figure,imagesc(BGValsIm); colorbar;

FlValsSum = sum(FlVals~=0,3);
FlValsIm  = sum(FlVals,3)./FlValsSum;
AFValsSum = sum(AFVals~=0,3);
AFValsIm  = sum(AFVals,3)./AFValsSum;

end