function A = ReadImageSeries( fileFolder, fileNames )

allFiles   = dir(fullfile(fileFolder,fileNames));
fileNames  = {allFiles.name}';
numImages  = numel(fileNames);

%Allocate space
A = imread(strcat(fileFolder,fileNames{1}));
A = zeros( size(A,1), size(A,2), numImages, 'uint16' );
for k = 1:numImages
    A(:,:,k) = imread(strcat(fileFolder,fileNames{k}));
end

end