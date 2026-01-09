function Y = SavaMedian3D(X)
%-------------------------------
%3-D MEDIAN FILTERING
%--------------------------------
%Y-3D output image (x,y,z)
%X--input image(noisy image)
Y=X;
[Ni,Nj,Nk]=size(X);

% do 2D median filtering on 3D matrix faces
Ayz = squeeze(X(1,:,:));
Y(1,1:Nj,1:Nk) = medfilt2(Ayz, [3 3]);
Ayz = squeeze(X(Ni,:,:));
Y(Ni,1:Nj,1:Nk) = medfilt2(Ayz, [3 3]);

Axz = squeeze(X(:,1,:));
Y(1:Ni,1,1:Nk) = medfilt2(Axz, [3 3]);
Axz = squeeze(X(:,Nj,:));
Y(1:Ni,Nj,1:Nk) = medfilt2(Axz, [3 3]);

Axy = squeeze(X(:,:,1));
Y(1:Ni,1:Nj,1) = medfilt2(Axy, [3 3]);
Axy = squeeze(X(:,:,Nk));
Y(1:Ni,1:Nj,Nk) = medfilt2(Axy, [3 3]);

% do 3D median filter inside the matrix
for i=2:Ni-1,
    i
    for j=2:Nj-1,
        for k=2:Nk-1,
            aa=X(i-1:i+1,j-1:j+1,k-1:k+1);
            Y(i,j,k)=median(aa(:));
        end;
    end;
end;
