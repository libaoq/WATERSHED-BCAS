function CoregisterStackAndPO2points_SingleImages_SelectedPoints
%%

% cretae slash/backslash depending on PC/UNIX 
if ispc,
    bsl = '\';
else
    bsl = '/';
end;

% first select root directory where individual planes are stored
planeDirName = uigetdir(pwd, 'Select directory name where data from one plane is stored');

cd(planeDirName);
pause(0);
[SurveyFileName, SurveySlicePathName] = uigetfile('.txt', 'Select Survey Image');
survFrm = load([SurveySlicePathName SurveyFileName],'-ascii');
load('image_parameters.txt','-ascii');

[SurveyFn1, SurveyPn1] = uigetfile('.jpg', 'Select the Same JPEG Survey Image');
survJPG = imread([SurveyPn1 SurveyFn1],'jpg');
image_parameters(1)=size(survJPG,1);
image_parameters(2)=size(survJPG,2);


survFrm = flipud(rot90(reshape(survFrm,image_parameters(1),image_parameters(2))));
%survFrm = reshape(survFrm,[image_parameters(1) image_parameters(2)]);
%survFrm = survFrm';
survFrm = medfilt2(survFrm, [2 2],'symmetric');

selected_points_indices = load([planeDirName bsl 'selected_points.txt']);
[dStr_fname dStr_pname] = uigetfile('*.mat','Select dataStruct.mat file','dataStruct_delay15.mat');
load([dStr_pname dStr_fname]); %load dataStruct structure...
% extract relative errors...
%%
pO2_relErr = dataStruct.pO2_rel_err;
%%

%%

[fitFileName fitDirName] = uigetfile('*.txt','Select file with Sergeis fitted pO2 values. Delay is ');
NewFitPO2 = load([fitDirName fitFileName],'-ascii');
dataStruct.pO2 = squeeze(NewFitPO2(:,2));

% plot images

fitcFig = figure(102);

bds = dataStruct.surveyImageGrayBounds;
survFrm(find(survFrm(:)<bds(1)))=bds(1);
survFrm(find(survFrm(:)>bds(2)))=bds(2);

hFITC = imagesc(survFrm);  colormap gray;cbFITC = colorbar;
set(gca,'XTick',[]); set(gca,'YTick',[]);
%set(get(gca,'YLabel'),'FontSize',15);set(get(gca,'YLabel'),'FontWeight','Bold'); set(get(gca,'XLabel'),'FontSize',15); set(get(gca,'XLabel'),'FontWeight','Bold');
set(gca,'FontSize',35);set(gca,'FontWeight','Bold');
%title('pO_2 (mmHg) Points Over FITC MIP');
set(fitcFig,'Position',[150 150 600 500]);
set(gca,'Position',[0.07 0.1 0.8 0.8]); %set(gca,'nextplot','replacechildren');
daspect([1 1 1]);
  

% change colormap...

colormap([gray(64);jet(64)]);

% set PO2 image
%C_PO2 = get(hPO2,'CData');
ncmp=64;

maxFTC = max(survFrm(:));
minFTC = min(survFrm(:));
c_FITC = min(ncmp,  round(   (ncmp-1).*(survFrm-minFTC)./(maxFTC-minFTC)  )+1  );
c_FITC = c_FITC;

set(hFITC,'CDataMapping','Direct');
set(hFITC,'CData',c_FITC);
pause(0);
figure(fitcFig);
caxis([min(c_FITC(:)) max(c_FITC(:))]);
pause(0);
set(cbFITC,'YLim',[ncmp+1 2*ncmp]);

%%
% set colormap values for PO2 points
% first set selected points only
% remove points with relErr > 0.1
pO2_relErr = abs(pO2_relErr); % just in case that it was not already done
selIdx = find(pO2_relErr < 0.15);
%selIdx = 1:length(dataStruct.pO2); %[4 9 12 13 14 16 17 20 21 22 23 24 25];
%selIdx = [4 9 12 13 14 16 17 19 20 22 23 24 25];
%selIdx([1 2 21 26 27])=[];
%%
dataStruct.pO2sel = dataStruct.pO2(selIdx);

% then create proper colormap
cmap1=jet(64);
CmapN = size(cmap1,1); % length of colormap[N,3] matrix is N
minPO2 = min(dataStruct.pO2sel);
maxPO2 = max(dataStruct.pO2sel);

%minPO2 = 5; 
%maxPO2 = 50;
dataStruct.pO2sel( find(dataStruct.pO2sel > maxPO2) )=maxPO2;
dataStruct.pO2sel( find(dataStruct.pO2sel < minPO2) )=minPO2;

cmapPO2idx = ones(length(dataStruct.pO2sel),1); % vector of colormap indexes for pO2
if maxPO2 == minPO2,
    cmapPO2idx = cmapPO2idx.*round(CmapN/2);
else
    cmapPO2idx = round(   (  CmapN.*(dataStruct.pO2sel-minPO2)+maxPO2-dataStruct.pO2sel  )    ./(maxPO2-minPO2) );
    cmapPO2idx(find(cmapPO2idx <= 0)) = 1;
    cmapPO2idx(find(cmapPO2idx > CmapN)) = CmapN;
end;


for i = 1:length(dataStruct.pO2),
    if ismember(i,selIdx),
        colMap = squeeze(cmap1(cmapPO2idx(find(selIdx==i)),:));
        rectangle('Position',[(selected_points_indices(i,3) - 2),(selected_points_indices(i,2) - 2), 4, 4],'Curvature',[0,0],'FaceColor',colMap,'LineStyle','none');
        daspect ([1,1,1])
        %text((new_selected_points_indices(i,3) - 5),(new_selected_points_indices(i,2) - 5),[num2str(dataStruct.pO2(i),3) ],'FontSize',10,'FontWeight', 'Bold','Color',dotColor);
    end;
end;

step = 10;         
newTickLabels = step*ceil(minPO2/step)  : step  : step*floor(maxPO2/step);
newTicks =( newTickLabels - (minPO2*ncmp-maxPO2)/(ncmp-1) ) .* (ncmp-1) ./ (maxPO2-minPO2);

%newTicks = 10:10:ncmp;
%newTickLabels = newTicks.*(maxPO2-minPO2)./(ncmp-1)+(minPO2*ncmp-maxPO2)/(ncmp-1);
newTickLabelsStr=num2str(newTickLabels(1),2);
for idxLab = 2:length(newTickLabels),
    newTickLabelsStr = char(newTickLabelsStr,num2str(newTickLabels(idxLab),2));
end;
pause(0);
set(cbFITC,'YTick',newTicks+ncmp);
set(cbFITC,'YTickLabel',newTickLabelsStr);
 

%%

%imageOutName = [planeDirName bsl 'pO2_pointsOverMIP' fnm(11:end-4) '.jpg'];
%print(tmpFig,'-djpeg',imageOutName,'-r600');
imageOutName = [planeDirName bsl 'pO2_pointsOverSurvey_Delay15.ps'];
print(fitcFig,'-dpsc2',imageOutName);

%%

% save pO2 values, errors, and x-y positions in txt file...
toSave_pO2 = dataStruct.pO2(selIdx);
toSave_XY  = selected_points_indices(selIdx,1:3);
toSave_err = pO2_relErr(selIdx);
toSave = [toSave_XY toSave_pO2 toSave_err];
save 'forImView_idx_XY_pO2_relErr.dat' toSave -ascii;
%%

decayDelay = dataStruct.delayDecayProcessing_us;
transformStruct.transformedSelectedPointsIndices = new_selected_points_indices;
transformStruct.FITCmipFrame = bareFITCframe;
transformStruct.newSurveyImage = newSurveyImage;
transformStruct.xData = xdata_po2;
transformStruct.yData = ydata_po2;
dataStruct.transformStruct=transformStruct;
save([planeDirName bsl fnm],'dataStruct','-mat');



% do image transformation on frm
% ask for the file name






                
                
