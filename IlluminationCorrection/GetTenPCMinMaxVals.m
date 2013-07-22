function OutVal = GetTenPCMinMaxVals( rows, cols, AvImg, Flag )

vals = AvImg(rows, cols);

if Flag
    vals = vals(:);
    vals(isnan(vals)) = 0;
    vals = sort(vals(:),'descend');
else
    vals = sort(vals(:),'ascend');
end

OutVal = mean(vals(1:10));

end