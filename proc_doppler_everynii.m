close all; clear; clc

%% load data
[filename,datapath]  = uigetfile('*.NII','Add data','MultiSelect','on');
filelist = filename;

bsurf = 140;
bdep = 128;

for file_i = 1:numel(filelist)

    % read data
    NII = ReadNII3([datapath,filelist{file_i}]);  
    RR = NII2RR(NII);
    [nz,nx,ny] = size(RR);
    clear NII

    file_str = filelist{file_i};
    file_str = file_str(15:17);
    
    % calculate doppler planes
    ilst = 35:nz; % z indices for motion correction
    indices = bsurf:bsurf+bdep-1; % these are the Z indices to calculate doppler on

    % widths = [];
    d0 = [];
    
    hwait = waitbar(0,sprintf('Calculating slice 0 of %d',ny));
    for ii=1:size(RR,3)
        waitbar(ii/size(RR,3),hwait,sprintf('Calculating slice %d of %d',ii,size(RR,3)));
        [a, d0(:,:,ii)] = octopusDoppler_step1(RR(:,:,ii), ilst, indices); 
        % [a, d0(:,:,ii)] = getDoppler('getPhase', RR(:,:,ii), ilst, indices, 1, 46780, 1, 3);
    end
    close(hwait);

    % save d0
	eval(['dangles_',file_str,' = d0;'])
    eval(['save ','dangles_', file_str,'.mat', ' dangles_', file_str, ' -v7.3'])
    
    % clear d0
    eval(['clear ', 'dangles_', file_str])

    % break;
end

%% display doppler planes
cm1a = [32:-1:1]'*[0 1 0]/32; cm1a(:,3) = 1;
cm1b = [32:-1:1]'*[0 0 1]/32;
cm2 = hot(64);

dangles = dangles_001;

figure;

colormap([cm1a; cm1b; cm2])
for ii=1:size(dangles,1)

    subplot(1,1, 1)
    foo = squeeze(dangles(ii,:,:));
    % foo = medfilt2(angles(:,:,ii), [3 8]);
    % imagesc(foo,[-1 1])
    imagesc(foo)
    title(num2str(ii))
    % title( sprintf('slice %d',idxY(ii)) )
    
    %    ha(2)=subplot(2,2,4);
    %    imagesc([1:size(G,1)],[1:8:512*8],squeeze(G(:,:,idxY(ii))),[0 20])
    
    pause(1)
end
