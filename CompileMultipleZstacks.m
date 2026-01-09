function CompileMultipleZstacks

% cretae slash/backslash depending on PC/UNIX 
if ispc,
    bsl = '\';
else
    bsl = '/';
end;

% first select root directory where Z stacks are stored
grandTopDirName = uigetdir(pwd, 'Select root directory name where individual Z Stacks are stored');
% then select multiple directories with XY frames
gd=dir(grandTopDirName);
gstrng = {gd.name};
% manual selection of frames
[gs,gv] = listdlg('ListString',gstrng,'PromptString','Select Z-stack directories:','SelectionMode','multiple','ListSize',[250 300]);
numZstacks = length(gs);
for ii = 1:numZstacks,
    workDirName = [grandTopDirName bsl gstrng{gs(ii)}];
    disp(['working in ' gstrng{gs(ii)} ' directory...']);
    CompileZstackFromSingleFrames(workDirName);
end;
    
    

end % end function

function ZframeStruct = CompileZstackFromSingleFrames(topDirName)

% select directories with data to process...
% select multuiple directories, each containing single XY frame
% data. This procedure will go in each of them, pick up individual
% frames into matrices, and save matlab structure with data matrix and
% frame details



% cretae slash/backslash depending on PC/UNIX 
if ispc,
    bsl = '\';
else
    bsl = '/';
end;

% first select root directory where Z and T series directories are stored
%topDirName = uigetdir(pwd, 'Select root directory name where individual Z frames are stored');

% then select multiple directories with XY frames
d=dir(topDirName);
strng = {d.name};
% manual selection of frames
%[s,v] = listdlg('ListString',strng,'PromptString','Select directories with data:','SelectionMode','multiple','ListSize',[250 300]);
% automatic selection of frames
i=0;
for idx = 1:length(strng),    
    if (~strcmp(strng{idx},'.') && ~strcmp(strng{idx},'..') ), 
        i=i+1;
        s(i) = idx;
    end;
end;

% get some common characteristics of the files (#pixels, #channels, #z
% frames...
pathName = [topDirName bsl strng{s(1)}];
disp(['Getting Frame Parameters from directory ' strng{s(1)} ' ...']);
xmlFname = [strng{s(1)} '.xml'];
fid = fopen([pathName bsl xmlFname]);
C = textscan(fid, '%s');
fclose(fid);
StringCounter = 0;
ReachedEndOfFile = 0;

while ReachedEndOfFile == 0,
    
    StringCounter = StringCounter + 1;
    ch=char(C{1,1}(StringCounter));
    
    if strcmp(ch,'</LVData>'), % if end of file, then exit
        ReachedEndOfFile = 1;
    else
        if strncmp(ch,'<Name>Required',14) & strncmp(char(C{1,1}(StringCounter+3)),'B-scans</Name>',14), % if string is on Numer of B scans
            ch = char(C{1,1}(StringCounter+4));
            Marks1 = strfind(ch,'>'); % <Val>512</Val>
            Marks2 = strfind(ch,'<');
            ZframeStruct.NumBscans = str2num(ch(Marks1(1)+1:Marks2(2)-1));  % get number of B scans per frame
        end;
        if strncmp(ch,'<Name>Required',14) & strncmp(char(C{1,1}(StringCounter+3)),'A-scans/B-scan</Name>',21), % <Name>Required No. of A-scans/B-scan</Name>
            ch = char(C{1,1}(StringCounter+4));
            Marks1 = strfind(ch,'>'); % <Val>512</Val>
            Marks2 = strfind(ch,'<');
            ZframeStruct.NumAscans = str2num(ch(Marks1(1)+1:Marks2(2)-1));  % get number of B scans per frame
        end;
        if strncmp(ch,'<Name>Number',12) & strncmp(char(C{1,1}(StringCounter+2)),'Active',6) & strncmp(char(C{1,1}(StringCounter+3)),'Channels</Name>',15), % <Name>Number of Active Channels</Name>
            ch = char(C{1,1}(StringCounter+4));
            Marks1 = strfind(ch,'>'); % <Val>512</Val>
            Marks2 = strfind(ch,'<');
            ZframeStruct.NumActiveChannels = str2num(ch(Marks1(1)+1:Marks2(2)-1));  % get number of B scans per frame
        end;
        if strncmp(ch,'<Name>Active',12) & strncmp(char(C{1,1}(StringCounter+1)),'Channel',7) & strncmp(char(C{1,1}(StringCounter+2)),'Array</Name>',12), % <Name>Active Channel Array</Name>
            StringCounter = StringCounter+3;
            for chIdx = 1:ZframeStruct.NumActiveChannels,
                StringCounter = StringCounter+3;
                ch = char(C{1,1}(StringCounter));
                StringCounter = StringCounter+1;
                Marks1 = strfind(ch,'>'); % <Val>512</Val>
                Marks2 = strfind(ch,'<');
                ZframeStruct.ChannelNumber(chIdx) = str2num(ch(Marks1(1)+1:Marks2(2)-1));  % get number of B scans per frame
            end;
        end;
    end;
end;

ZframeStruct.NumZframes = length(s);
if ZframeStruct.NumActiveChannels > 1,
    ZframeStruct.Data = zeros(ZframeStruct.NumZframes,ZframeStruct.NumActiveChannels,ZframeStruct.NumAscans,ZframeStruct.NumBscans);
else
    ZframeStruct.Data = zeros(ZframeStruct.NumZframes,ZframeStruct.NumAscans,ZframeStruct.NumBscans);
end;
    

        

            

% loop through the list of selected directories, and inside each of them
% take out individual frames


for i=1:length(s),
    pathName = [topDirName bsl strng{s(i)}];
    disp(['Processing directory ' strng{s(i)} ' ...']);
    fName = [strng{s(i)} '.bin'];
    
    fid = fopen([pathName bsl fName]);
     % import binary file into Data matrix...
    if ZframeStruct.NumActiveChannels > 1,
        xx=fread(fid,[ZframeStruct.NumActiveChannels,ZframeStruct.NumAscans,ZframeStruct.NumBscans],'int32');
        for idxChan = 1:ZframeStruct.NumActiveChannels,
            xxy = squeeze(xx(idxChan,:,:));
            ZframeStruct.Data(i,idxChan,:,:) = flipud(rot90(xxy));
        end;
    else
        xx=fread(fid,[ZframeStruct.NumAscans,ZframeStruct.NumBscans],'int32');
        ZframeStruct.Data(i,:,:) = flipud(rot90(xx));
    end;
    fclose(fid);  
    
    % get frame time from frame name
    timeString = strng{s(i)};
    timeString=timeString(end-9:end);
    ZframeStruct.frmTime_msec(i) = str2num(timeString(5:10)) + str2num(timeString(3:4))*60 + str2num(timeString(1:2))*3600;
    
    % import PMT gain, EOM voltage, Z position, frame time, etc.
    fName = [strng{s(i)} '.xml'];
    fid = fopen([pathName bsl fName]);
    C = textscan(fid, '%s');
    fclose(fid);
    StringCounter = 0;
    ReachedEndOfFile = 0;

    while ReachedEndOfFile == 0,
    
        StringCounter = StringCounter + 1;
        ch=char(C{1,1}(StringCounter));
    
        if strcmp(ch,'</LVData>'), % if end of file, then exit
            ReachedEndOfFile = 1;
        else
            if strncmp(ch,'<Name>EOM',9) & strncmp(char(C{1,1}(StringCounter+1)),'data</Name>',11), % if string is on EOM opening...
                ch = char(C{1,1}(StringCounter+2));
                Marks1 = strfind(ch,'>'); % <Val>512</Val>
                Marks2 = strfind(ch,'<');
                ZframeStruct.EOMvoltage(i) = str2num(ch(Marks1(1)+1:Marks2(2)-1));  % get number of B scans per frame
            end;
            %%if strncmp(ch,'<Name>Position</Name>',21), % XYZ position
            if strncmp(ch,'<Name>Value</Name>',18), % XYZ position
                StringCounter = StringCounter+2;
                for positionIdx = 1:3,
                    StringCounter = StringCounter+2;
                    ch = char(C{1,1}(StringCounter));
                    StringCounter = StringCounter+2;
                    Marks1 = strfind(ch,'>'); 
                    Marks2 = strfind(ch,'<');
                    ZframeStruct.XYZposition(i,positionIdx) = str2num(ch(Marks1(1)+1:Marks2(2)-1));  % get xyz position
                end;
            end;
            if strncmp(ch,'<Name>PMT',9) & strncmp(char(C{1,1}(StringCounter+1)),'true',4) & strncmp(char(C{1,1}(StringCounter+2)),'Voltages</Name>',15), % <Name>PMT true Voltages</Name>
                StringCounter = StringCounter+4;
                for chIdx = 1:4,
                    StringCounter = StringCounter+2;
                    ch = char(C{1,1}(StringCounter));
                    StringCounter = StringCounter+2;
                    Marks1 = strfind(ch,'>'); 
                    Marks2 = strfind(ch,'<');
                    ZframeStruct.PMTvolts(i,chIdx) = str2num(ch(Marks1(1)+1:Marks2(2)-1));  % get PMT voltages
                end;
            end;    
        end;
    end;   
    
    disp(['... done!']);
end;

%[fileOut, pathOut] = uiputfile( '*.mat','Save Collected Structures');
pathOut = [topDirName bsl];

[foo1, fileOut, foo2] = fileparts(pathOut(1:end-1));

fileOut = [fileOut '.mat'];

toSaveName = [pathOut fileOut];
save(toSaveName,'ZframeStruct','-mat');
end