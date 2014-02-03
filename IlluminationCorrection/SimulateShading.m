Im = zeros( 2048,2048,100,'uint16' );
Lab = zeros( 2048,2048,100,'uint8' );
CountIm = zeros( 2048,2048,3,'uint8' );
lamdas = [ 50, 200, 1000 ];
R = rand( 20*20*100, 1 ).*3+0.5.*ones(20*20*100,1);
R = round(R); R = cast(R,'int32');

%Create shading surfaces
%c = rand(14,1)-0.5; c = c./sum(abs(c));
[X,Y] = meshgrid(-1023:1024, -1023:1024);
X = X./1024; Y = abs(Y)./1024;
mul = [ 0.95, 0.93, 0.91]; add = [ 0, 0, 0 ]; scale = [200,400,800];
add = add./1024;
surfs = zeros( 2048,2048,3);
for i=1:3
    surfs(:,:,i) = c(1).*(add(i)+X).^4 + c(2).*(add(i)+X).^3.*Y +...
       c(3).*(add(i)+X).^2.*Y.^2 + c(4).*(add(i)+X).*Y.^3 + c(5).*Y.^4 +...
       c(6).*(add(i)+X).^3 + c(7).*(add(i)+X).^2.*Y +...
       c(8).*(add(i)+X).*Y.^2 + c(9).*Y.^3 + c(10).*(add(i)+X).^2 +...
       c(11).*(add(i)+X).*Y + c(12).*Y.^2 + c(13).*(add(i)+X) + c(14).*Y;
    surfsmin = min(min(surfs(:,:,i)));
    surfs(:,:,i) = surfs(:,:,i) - surfsmin.*ones(size(X,1),size(X,2));
    surfsmax = max(max(surfs(:,:,i)));
    surfs(:,:,i) = surfs(:,:,i)./surfsmax;
    figure,imshow(surfs(:,:,i));
    surfs(:,:,i) = surfs(:,:,i).*(1-mul(i))+ones(size(X,1),size(X,2)).*mul(i);
end

%WriteImageSeries(surfs, '/data/kedar/sim/ht/', 'ShSurf', 'tif');
tic
count = 1;
%Create image
for i=1:size(Im,3)
    for j=1:128:2048
        for k=1:128:2048
            RIm = poissrnd( lamdas(R(count)), 128,128 );
            Im(j:j+127,k:k+127,i) = ...
                cast( RIm.*surfs(j:j+127,k:k+127,R(count)), 'uint16' );
            Lab(j:j+127,k:k+127,i) = ones(128,128,'uint8').*...
                                     cast( (R(count)-1),'uint8');
            CountIm(j:j+127,k:k+127,R(count)) = ...
                CountIm(j:j+127,k:k+127,R(count))+ones(128,128,'uint8');
            count = count+1;
        end
    end
end
WriteImageSeries(Im, '/data/kedar/sim/ht/sim/', 'ShSim', 'tif');
WriteImageSeries(Lab, '/data/kedar/sim/ht/sim_label/', 'ShSimlab', 'tif');
% WriteImageSeries(CountIm, '/data/kedar/sim/ht/', 'Counts', 'tif');
toc