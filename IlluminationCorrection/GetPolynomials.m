function poly = GetPolynomials(FlValsIm, AFValsIm, BGValsIm)

%Store the indices and the values
Fl = log(FlValsIm(~isnan(FlValsIm)));
AF = log(AFValsIm(~isnan(AFValsIm)));
BG = log(BGValsIm(~isnan(BGValsIm)));
[xFl,yFl] = find(~isnan(FlValsIm));
[xAF,yAF] = find(~isnan(AFValsIm));
[xBG,yBG] = find(~isnan(BGValsIm));

%Put the above into arrays
%Using third order exponential
FlInd = [ xFl(:).^4 yFl(:).^4 xFl(:).*yFl(:).^3      ... 4th
             xFl(:).^3.*yFl(:) xFl(:).^2.*yFl(:).^2  ... 4th continued
          xFl.^3 yFl.^3 xFl.*yFl.^2 xFl.^2.*yFl      ... 3rd
          xFl.^2 yFl.^2 xFl.*yFl                     ... 2nd
          xFl yFl ];            %ones(size(xFl,1),1)     1st
AFInd = [ xAF(:).^4 yAF(:).^4 xAF(:).*yAF(:).^3      ... 4th
             xAF(:).^3.*yAF(:) xAF(:).^2.*yAF(:).^2  ... 4th continued
          xAF.^3 yAF.^3 xAF.*yAF.^2 xAF.^2.*yAF      ... 3rd
          xAF.^2 yAF.^2 xAF.*yAF                     ... 2nd
          xAF yAF ];            %ones(size(xAF,1),1)1st
BGInd = [ xBG(:).^4 yBG(:).^4 xBG(:).*yBG(:).^3      ... 4th
             xBG(:).^3.*yBG(:) xBG(:).^2.*yBG(:).^2  ... 4th continued
          xBG.^3 yBG.^3 xBG.*yBG.^2 xBG.^2.*yBG      ... 3rd
          xBG.^2 yBG.^2 xBG.*yBG                     ... 2nd
          xBG yBG ];            %ones(size(xBG,1),1)     1st

%Normalize inputs for Autofluorescence
meanAFInd = mean(AFInd);            meanAFVal = mean(AF);
normAFInd = sqrt(sum(AFInd.^2,1));  normAFVal = sqrt(sum(AF.^2,1));
AFInd=AFInd-repmat(meanAFInd,[size(AFInd,1) 1]);
AFInd=mexNormalize(AFInd);
AFVal=AF-repmat(meanAFVal,[size(AF,1) 1]);
AFVal=mexNormalize(AFVal);

WAF = RegressionWithReguAndFTest( AFInd, AFVal );
%WAF = WAF.*normAFInd';

%Normalize inputs for Fluorescence
meanFlInd = mean(FlInd);            meanFlVal = mean(Fl);
normFlInd = sqrt(sum(FlInd.^2,1));  normFlVal = sqrt(sum(Fl.^2,1));
FlInd=FlInd-repmat(meanFlInd,[size(FlInd,1) 1]);
FlInd=mexNormalize(FlInd);
FlVal=Fl-repmat(meanFlVal,[size(Fl,1) 1]);
FlVal=mexNormalize(FlVal);

WFl = RegressionWithReguAndFTest( FlInd, FlVal );
% WFl = WFl.*normFlInd';

%Normalize inputs for the background
meanBGInd = mean(BGInd);            meanBGVal = mean(BG);
normBGInd = sqrt(sum(BGInd.^2,1));  normBGVal = sqrt(sum(BG.^2,1));
BGInd=BGInd-repmat(meanBGInd,[size(BGInd,1) 1]);
BGInd=mexNormalize(BGInd);
BGVal=BG-repmat(meanBGVal,[size(BG,1) 1]);
BGVal=mexNormalize(BGVal);

WBG = RegressionWithReguAndFTest( BGInd, BGVal );
% WBG = (BGInd'*BGInd)\(BGInd'*BGVal);
% WBG = WBG.*normBGInd';

poly = struct;
poly(1).coeffs  = WFl';
poly(1).meanInd = meanFlInd;
poly(1).normInd = normFlInd;
poly(1).meanVal = meanFlVal;
poly(1).normVal = normFlVal;
poly(2).coeffs  = WAF';
poly(2).meanInd = meanAFInd;
poly(2).normInd = normAFInd;
poly(2).meanVal = meanAFVal;
poly(2).normVal = normAFVal;
poly(3).coeffs  = WBG';
poly(3).meanInd = meanBGInd;
poly(3).normInd = normBGInd;
poly(3).meanVal = meanBGVal;
poly(3).normVal = normBGVal;

end