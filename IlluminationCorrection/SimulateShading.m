Im = zeros( 2048,2048,100,'uint16' );
lamdas = [ 500, 2000, 10000 ];
R = rand( 20*20*100, 1 ).*3+0.5.*ones(20*20*100,1);
R = round(R); R = cast(R,'int32');
count = 1;

%Create shading surfaces
[X,Y] = meshgrid(-1023:1024, -1023:1024);
X = abs(X)./1024; Y = abs(Y)./1024;
c = rand(14,1)-0.5; c = c./sum(abs(c));
mul = [ 1, 1, -1]; add = [ 0, 5, 5 ]; scale = [100,500,1000];
surfs = zeros( 2048,2048,3);
for i=1:3
    surfs(:,:,i) = c(1).*X.^4 + c(2).*X.^3.*Y + c(3).*X.^2.*Y.^2 + ...
                   c(4).*X.*Y.^3 + c(5).*Y.^4 +  ...
         c(6).*X.^3 + c(7).*X.^2.*Y + c(8).*X.*Y.^2 + c(9).*Y.^3 + ...
         c(10).*X.^2 + c(11).*X.*Y + c(12).*Y.^2 + ...
         mul(i).*(add(i)+c(13)).*X + mul(i).*(add(i)+c(14)).*Y;
    surfsmin = min(min(surfs(:,:,i)));
    surfs(:,:,i) = surfs(:,:,i) - surfsmin.*ones(size(X,1),size(X,2));
    surfs(:,:,i) = surfs(:,:,i)./max(max(surfs(:,:,i))).*scale(i);
end
surfs = cast(surfs,'uint16');

%Create image
for i=1:size(Im,3)
    for j=1:128:2048
        for k=1:128:2048
            RIm = cast( poissrnd( R(count), 128,128 ), 'uint16' );
            Im(j:j+127,k:k+127,i) = RIm + surfs(j:j+127,k:k+127,R(count));
            count = count+1;
        end
    end
end
WriteImageSeries(Im, 'E:/Kedar/', 'ShadingSimulated', 'tif');
WriteImageSeries(surfs, 'E:/Kedar/', 'ShadingSurfaces', 'tif');