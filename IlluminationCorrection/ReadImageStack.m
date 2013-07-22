function A = ReadImageStack( fileFolder, fileNames )

filename = strcat(fileFolder,fileNames);
info = imfinfo(filename);
numImages = numel(info);

%Allocate space
A = imread(filename, 1, 'Info', info);
A = zeros( size(A,1), size(A,2), numImages, 'uint8' );
for k = 1:numImages
    A(:,:,k) = imread(filename, k, 'Info', info);
end

end