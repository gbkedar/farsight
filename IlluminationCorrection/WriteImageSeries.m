function WriteImageSeries(imageStack, fileFolder, fileNameTemp, ext)
for k = 1:size(imageStack,3)
    kk = cast(k,'uint16');
    A = imageStack(:,:,kk);
    fileStr=strcat(fileFolder,fileNameTemp,num2str(kk-1,'%3.3d'),'.',ext);
    imwrite(A,fileStr,ext);
end
end