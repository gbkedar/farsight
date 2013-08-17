h = fspecial('gaussian', [4096 4096],30000);
h1 = h(1500:(1500+2047),1500:(1500+2047));
h1 = h1 - repmat(min(min(h1)),size(h1,1),size(h1,1));
h1 = h1 ./ repmat(max(max(h1)),size(h1,1),size(h1,1));
h1 = h1 .* repmat(15,size(h1,1),size(h1,1));
A = zeros( 2048, 2048, 100, 'uint8' );
for i=1:100
    h2 = h1 + repmat((100+round(rand(1)*15)),size(h1,1),size(h1,1));
    h2 = cast(h2,'uint8');
    h2 = imnoise(h2,'poisson');
    if i>1
        imwrite(h2, 'simulatedShading.tif','tif', 'Compression', 'none', 'WriteMode', 'append');
    else
        imwrite(h2, 'simulatedShading.tif','tif', 'Compression', 'none');
    end
end
figure,imshow(h2); colorbar