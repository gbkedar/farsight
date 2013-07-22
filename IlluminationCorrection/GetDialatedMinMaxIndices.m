function [row,col] = GetDialatedMinMaxIndices( InpIm, Flag)

minLength = 10.^2;

if Flag
    [row,col] = find(InpIm == max(InpIm(:)));
else
    [row,col] = find(InpIm == min(InpIm(:)));
end

%If the search region is too small then
newrow = row; newcol = col;
while length(newrow)<minLength
    row = newrow; col = newcol;
    tempIm = logical(size(InpIm));
    tempIm(row,col) = 1;
    se = strel('disk',2);
    tempIm = imdilate(tempIm,se);
    [newrow,newcol] = find(tempIm);
end

end