function CombineRepetativeScans

% cretae slash/backslash depending on PC/UNIX 
if ispc,
    bsl = '\';
else
    bsl = '/';
end;

% first select root directory where individual scans are stored
topDirName = uigetdir(pwd, 'Select root directory name where individual point scans are stored');

% then select multiple directories with XY frames
d=dir(topDirName);
strng = {d.name};
[s,v] = listdlg('ListString',strng,'PromptString','Select directories with single scans:','SelectionMode','multiple','ListSize',[250 300]);

numDirectories = length(s);

% go to the top directory (below are individual point scans)
pathComb = [topDirName bsl 'Combined'];
currDir = pwd;
%cd(topDirName);
%[fileComb, pathComb] = uiputfile( '*.mat','Create Combined Directory (type some file name too...)');

for idxPtScans = 1:numDirectories,
    tmpDirName = [topDirName bsl strng{s(idxPtScans)}];
    cd(tmpDirName);
    load 'binned_data.txt' -ascii;
    load 'image_parameters.txt' -ascii;
    if idxPtScans == 1,
        % get/create bigger matrix for combined scans
        [nScans nTimePts] = size(binned_data);
        totNumScans = nScans-1;
        binned_data_combined = zeros(numDirectories*totNumScans+1,nTimePts);
        binned_data_combined(1,:)=binned_data(1,:);
        % get parameter files
        load 'image_parameters.txt' -ascii;
        interExcTime = image_parameters(4);
        load 'selected_points.txt' -ascii;
        load 'survey_scan_image.txt' -ascii;
        nPts = size(selected_points,1);
        nReps = totNumScans/nPts;
        image_parameters(5) = numDirectories*nReps;
        survey_scan_image = survey_scan_image .* numDirectories;
        %tmpData = zeros(nPts,nReps*numDirectories,nTimePts);   
    end;
    % if different trials have diferent interexcitation times, take a
    % minimal time in all
    nTimePts = min( nTimePts, size(binned_data,2));
    bined_data = binned_data(:,1:nTimePts);
    binned_data_combined = binned_data_combined(:,1:nTimePts);
    interExcTime = min( interExcTime, image_parameters(4));
    image_parameters(4) = interExcTime;
    image_parameters(5) = numDirectories*nReps;
    
    for ii = 1:nPts,
        %tmpData(ii,:,:) = binned_data((ii-1)*nReps+2:ii*nReps+1,:);
        binned_data_combined(  (ii-1)*nReps*numDirectories + (idxPtScans-1)*nReps + 2 : (ii-1)*nReps*numDirectories + idxPtScans*nReps + 1, 1:nTimePts ) = binned_data( (ii-1)*nReps+2 : ii*nReps+1 ,1:nTimePts);
    end;
    %binned_data_combined((idxPtScans-1)*(nScans-1)+2:idxPtScans*(nScans-1)+1,:)=binned_data(2:end,:);
end;

cd(pathComb);
save('binned_data.mat','binned_data_combined','-mat','-v7.3');
save('image_parameters.txt','image_parameters','-ascii');
save('selected_points.txt','selected_points','-ascii');
cd(currDir);


