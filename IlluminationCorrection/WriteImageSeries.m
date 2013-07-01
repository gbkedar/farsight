function WriteImageSeries(imageStack, fileFolder, fileNameTemp, extension)
for k = 1:size(imageStack,3)
    A = imageStack(:,:,k);
    imwrite(A,strcat(fileFolder,fileNameTemp,num2str(k-1,3),extension));
end
end