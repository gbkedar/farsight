function SubSurf = GetRescaledSurface( InpSurf, Min, Max )

diff = Max-Min;
Min = min(InpSurf(:));
SubSurf = InpSurf-repmat(Min,[size(InpSurf,1) size(InpSurf,2)]);
Max = max(SubSurf(:));
ratio = diff/Max;
SubSurf = SubSurf.*repmat(ratio,[size(InpSurf,1) size(InpSurf,2)]);

end