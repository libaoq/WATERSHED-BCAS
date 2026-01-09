function adjustTIFFhistogram
%%
% [fname pname] = uigetfile('*.*','Select 3D Tiff Stack');
% %%
% hfig = imfinfo([pname fname]);
% Nz = size(hfig,1);
% Nx = hfig(1).Height;
% Ny = hfig(1).Width;
% data = uint16(zeros(Nz,Nx,Ny));
% %%
% for ii=1:Nz,
%     data(ii,:,:) = imread([pname fname],'tiff',ii);
% end;

% dataMed3d = data;
load('HandlesTotalStack.mat')
dataMed3d = TotalStack1;
[Nz,Nx,Ny] = size(dataMed3d);
%%
% dataMed3d = medianImageStack(double(data));
%%
data_noBackground = dataMed3d;
for ii=1:Nz,
    ii
    foo = squeeze(dataMed3d(ii,:,:));
    foo1 = uint16(filter2(fspecial('average',50),foo));
    %foo1 = imfilter(foo,fspecial('average',50),'symmetric');
    foo = foo-foo1;
    foo(find(foo(:)<0))=0;
    data_noBackground(ii,:,:) = foo;
end;
%% don't use this section for now
data_final = double(data_noBackground);
for ii=1:Nz,
    foo = double(squeeze(data_noBackground(ii,:,:)));
    foo = foo./max(foo(:));
    poo = sort(foo(:));
    poo(find(poo==0))=[];
    Im = poo(round(length(poo)*0.999));
    boo = imadjust(foo,[0 Im],[0 1]);
    data_final(ii,:,:) = boo;
end;
save data_final.mat data_final
% ----------up here
%%
data_final(1:10,:,:)=[];
figure; imshow(squeeze(data_final(10,:,:)),[]);
% save data
% fname1 = [fname(1:end-4) '_adjusted.mat'];
% save([pname fname1],'data_final','-mat');

save data_final.mat data_final

%% display
figure('color','w');
rangei = round(linspace(1,231,17));
j=1;
%for i = 2:length(rangei)
for i = 14:17
    subplot(1,4,j);
    j=j+1;
    foo = squeeze(max(data_final(rangei(i-1):rangei(i),:,:),[],1));
    imshow(foo,[]);
    axis off
    title(['MIP: ', num2str(rangei(i-1)*2-2), '-', num2str(rangei(i)*2), 'um']);
end

