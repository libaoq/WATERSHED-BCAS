function CorrectZinZstacks

% cretae slash/backslash depending on PC/UNIX 
if ispc,
    bsl = '\';
else
    bsl = '/';
end;

% first select root directory where Z stacks are stored
topDirName = uigetdir(pwd, 'Select directory name where individual Z Stacks are stored');
% then select multiple directories with XY frames
d=dir(topDirName);
strng = {d.name};
% manual selection of frames
[s,v] = listdlg('ListString',strng,'PromptString','Select Z-stacks:','SelectionMode','multiple','ListSize',[250 300]);

%%
for ii = 1:length(s),
    load([topDirName bsl strng{s(ii)}]);
    Z{ii} = 25-squeeze(ZframeStruct.XYZposition(:,3));
end;
[Nz Nx Ny] = size(ZframeStruct.Data);
%%

nZ = length(s);
for ii = 1:length(s),
    foo = Z{ii};
    z1(ii) = foo(1);
end;
[z1ascend, idxAsc]=sort(z1);
for ii = 1:length(s),
    strngAsc{ii} = strng{s(idxAsc(ii))};
end;
% above for loop creates strngAsc(1:nZ) with proper order of stack names
% and z1ascend has correct starting Z positions

%%
% now pick up data and save files
z1ascend_um = round(z1ascend.*1000);
for ii = 1:length(s),
    load([topDirName bsl strngAsc{ii}]);
    ChNumber = ZframeStruct.ChannelNumber;
    data = ZframeStruct.Data;
    [Nz,Nx,Ny] = size(data); 
    Zi = 25000 - round(1000.*squeeze(ZframeStruct.XYZposition(:,3))) - z1ascend_um(1)+1;
    zstep = round(Zi(2)-Zi(1));
    if zstep > 1, % Z steps are greater than 1 um
        Zi_frames = round((Zi+1)/2);
    else
        Zi_frames = Zi;
    end;
        
    clear 'ZframeStruct';
    
    % give the name to the chunk and save it...
    zChunkName = ['Ch' num2str(ChNumber) '_Z_' sprintf('%05.0f',min(Zi_frames)) '_' sprintf('%05.0f',max(Zi_frames)) '_' ...
    sprintf('%05.0f',min(Zi)) '_' sprintf('%05.0f',max(Zi)) '.mat'];      % Ch3_Z_00859_00969_00101_00211.mat
    eval(['Ch' num2str(ChNumber) ' = data;']);
    dName =['Ch' num2str(ChNumber)];
    save([topDirName bsl zChunkName], dName, '-mat');
    clear dName;
end;
    




    
