function Convert3DmatToTIFF(X, filename)

X = X./max(X(:)).*65535;

[Nz, Nx, Ny] = size(X);

for i=1:Nz,
    fr = uint16(squeeze(X(i,:,:)));
    imwrite(fr,filename,'tiff','Compression','none','WriteMode','append');
end;
