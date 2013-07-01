%Add path to compiled SPAMS library
SpamsFolder = 'E:/SharedFolder/ARBs_MMI-31_28R_p4067/spams-matlab/build/';
addpath(SpamsFolder);
tilesPerRow = 42;
clear SpamsFolder

fileFolder = 'E:/SharedFolder/ARBs_MMI-31_28R_p4060/C2/';
fileNames1  = 'ARBs_MMI-31_28R_p4060C2*.tif';
fileName2  = 'out.tif';
AF = ReadImageSeries( fileFolder, fileNames1 );
TF = ReadImageStack( fileFolder, fileName2 );

% Poly = GetPolynomials( AF, TF, tilesPerRow );