function ExtractLifetimeDataIntoDATfile
%%

%%
% cretae slash/backslash depending on PC/UNIX 
if ispc,
    bsl = '\';
else
    bsl = '/';
end;
% select dataStruct_delay15.mat file with individual decays...
[fname pname] = uigetfile('*.mat','Select DataStruct mat file with decays to extract.','DataStruct_delay15.mat');

load([pname fname]);

t=dataStruct.fittedDecayTimePoints_us;
t=t';
data=squeeze(dataStruct.summed_lifetime_data(:,dataStruct.fittedDecayTimePoints_indices));
data = data';
toExport = [t data];

newFname = [fname(1:end-4) '_OnlyDecays.dat'];
save([pname newFname],'toExport','-ascii');

 