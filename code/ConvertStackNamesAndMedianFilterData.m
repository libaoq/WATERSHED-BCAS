function ConvertStackNamesAndMedianFilterData

% Load data files
[FileName,PathName] = uigetfile('*.mat','Select Matlab File Names With Data to filter', 'MultiSelect','on','Location',[500 500]);
if iscell(FileName),
    NumFiles = size(FileName,2);
else
    NumFiles = 1;
end;

%%
% find minimum and maximum Z position
for i=1:NumFiles,
    if iscell(FileName),
        fname = FileName{i};
    else
        fname = FileName;
    end;
    Zmin(i) = str2num(fname(7:11));
    Zmax(i) = str2num(fname(13:17));
end;
ZrelMin = Zmin-min(Zmin)+1;
ZrelMax = Zmax-min(Zmin)+1;
%%

% Do you want to median 3D filter files?
button = questdlg('Do you want 3D median filter?','Median 3D (Yes/No)','Yes','No','Yes');

% loop through the files and rename them / filter them...
for i=1:NumFiles,
    %rename old file
    if iscell(FileName),
        fname = FileName{i};
    else
        fname = FileName;
    end;
    disp(['Processing file ' fname ' ...']);
    newname = [fname(1:end-4) '_' sprintf('%05d',ZrelMin(i)) '_' sprintf('%05d',ZrelMax(i)) '.mat'];
    if iscell(FileName),
        NewFileName{i} = newname;
    else
        NewFileName = newname;
    end;
    movefile([PathName fname],[PathName newname]);
    
    % filter file...
    if strcmp(button,'Yes'),
        infname = [PathName newname];
        x=whos('-file',infname);
        load(infname);
        eval(['Data = ' x.name ';' 'clear ' x.name ';']);
        DataMed3D = medianImageStack(Data);
        clear Data;
        fname3d = [x.name '_Med3D'];
        eval([fname3d ' = DataMed3D; clear DataMed3D; ']);
        tosave = [infname(1:end-4) '_Med3D.mat'];
        save(tosave,fname3d);
    end;
    disp('... done!');
end;



