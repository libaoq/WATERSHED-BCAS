function varargout = GUI_ConnectStacks_Test(varargin)
% GUI_CONNECTSTACKS_TEST M-file for GUI_ConnectStacks_Test.fig
%      GUI_CONNECTSTACKS_TEST, by itself, creates a new GUI_CONNECTSTACKS_TEST or raises the existing
%      singleton*.
%
%      H = GUI_CONNECTSTACKS_TEST returns the handle to a new GUI_CONNECTSTACKS_TEST or the handle to
%      the existing singleton*.
%
%      GUI_CONNECTSTACKS_TEST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_CONNECTSTACKS_TEST.M with the given input arguments.
%
%      GUI_CONNECTSTACKS_TEST('Property','Value',...) creates a new GUI_CONNECTSTACKS_TEST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_ConnectStacks_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_ConnectStacks_Test_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_ConnectStacks_Test

% Last Modified by GUIDE v2.5 22-Mar-2010 11:32:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_ConnectStacks_Test_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_ConnectStacks_Test_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before GUI_ConnectStacks_Test is made visible.
function GUI_ConnectStacks_Test_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_ConnectStacks_Test (see VARARGIN)

% Choose default command line output for GUI_ConnectStacks_Test
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using GUI_ConnectStacks_Test.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end

% UIWAIT makes GUI_ConnectStacks_Test wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_ConnectStacks_Test_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)




% Load Individual Stacks
% --- Executes on button press in LoadStacksButton.
function LoadIndividualStacks(hObject, eventdata, handles)
% hObject    handle to LoadStacksButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%gets input file(s) from user
[input_file,pathname] = uigetfile({'*.mat','MAT-files (*.mat)'},'Select all stacks in correct order','MultiSelect', 'on');
%if file selection is cancelled, pathname should be zero and nothing should happen
if pathname == 0,
    return;
end;
%if they only select one file, then the data will not be a cell
%if more than one file selected at once,
%then the data is stored inside a cell
if iscell(input_file),
    handles.NumberOfStacks = length(input_file);
else
    handles.NumberOfStacks = 1;
end;
handles.PathName = pathname;
handles.LocalInputFileNames = input_file;

for n = 1:handles.NumberOfStacks,
    if handles.NumberOfStacks > 1,
        name = input_file{n};
    else
        name = input_file;
    end;
    name = name(1:end-4);
    tmp = load([pathname name '.mat']);
    StackName = fieldnames(tmp);
    eval(['stk = tmp.' StackName{1} ';']);
    stk = double(stk);
    clear tmp;
    handles.Stack{n} = stk;
    handles.StackDim(n,1:3) = size(stk);
    clear stk;
    handles.StackZplace(n,1) = str2num(name(7:11));
    if n == 1,
        FirstZ = handles.StackZplace(n,1);
    end;
    handles.StackZplace(n,2) = str2num(name(13:17));
    handles.StackZplace(n,1) = handles.StackZplace(n,1) - FirstZ +1;
    handles.StackZplace(n,2) = handles.StackZplace(n,2) - FirstZ +1;
    PopUpMenu2{n} = ['Z = ' num2str(handles.StackZplace(n,1)) ' -> ' num2str(handles.StackZplace(n,2))];
end;
set(handles.popupmenu2,'String',PopUpMenu2);
set(handles.popupmenu4,'String',PopUpMenu2);
handles.dN = round(double(handles.StackDim(1,2))/25);
set(handles.slider1,'Min',-handles.dN);
set(handles.slider1,'Max',handles.dN);
set(handles.slider2,'Min',-handles.dN);
set(handles.slider2,'Max',handles.dN);
guidata(hObject, handles);



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over LoadStacksButton.
function LoadStacksButton_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to LoadStacksButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
contents = get(hObject,'Value');
stk = handles.Stack{contents}; % stk is Z-X-Y matrix
contents1 = get(handles.popupmenu3,'Value');
if contents1 == 1,
    Img = squeeze(max(stk,[],1)); % project on X-Y plane
elseif contents1 == 2,
    Img = squeeze(max(stk,[],3)); % project on X-Z plane
else
    Img = squeeze(max(stk,[],2)); % project on Y-Z plane
end;
set(handles.edit3,'String',num2str(  min(min(min(Img))) ));
set(handles.edit4,'String',num2str(  max(max(max(Img))) ));
guidata(hObject,handles);
axes(handles.axes1);
imagesc(Img); 
colormap gray;
colorbar;


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4
contents = get(hObject,'Value');
n=contents;

ans=0;
if isfield(handles,'AdjustHoriz'),
    if n <= length(handles.AdjustHoriz),
        if ~isempty(handles.AdjustHoriz(n).Top),
            ans=1;
        end;
    end;
end;
if ans,
    axes(handles.axes1);
    imagesc(squeeze(max(handles.AdjustHoriz(n).Joint,[],1))); 
    colormap gray;
    colorbar;
    set(handles.edit1,'String',num2str(handles.AdjustHoriz(n).dX));
    set(handles.slider1,'Value',handles.AdjustHoriz(n).dX);
    set(handles.edit2,'String',num2str(handles.AdjustHoriz(n).dY));
    set(handles.slider2,'Value',handles.AdjustHoriz(n).dY);
else
    stk1 = handles.Stack{contents-1}; % stk is Z-X-Y matrix
    stk2 = handles.Stack{contents};
    Ztp1 = handles.StackZplace(contents-1,1);
    Ztp2 = handles.StackZplace(contents-1,2);
    Zbt1 = handles.StackZplace(contents,1);
    Zbt2 = handles.StackZplace(contents,2);
            
    Xpos = 0;
    handles.AdjustHoriz(n).dX = 0;
    Ypos = 0;
    handles.AdjustHoriz(n).dY = 0;
    set(handles.edit1,'String',num2str(handles.AdjustHoriz(n).dX));
    set(handles.edit2,'String',num2str(handles.AdjustHoriz(n).dY));
    set(handles.slider1,'Value',handles.AdjustHoriz(n).dX);
    set(handles.slider2,'Value',handles.AdjustHoriz(n).dY);
    handles.AdjustHoriz(n).Top = stk1(Zbt1-Ztp1+1:Ztp2-Ztp1+1,:,:);
    Top = handles.AdjustHoriz(n).Top;
    handles.AdjustHoriz(n).Btm = stk2(Zbt1-Zbt1+1:Ztp2-Zbt1+1,:,:);
    Btm = handles.AdjustHoriz(n).Btm;
    handles.AdjustHoriz(n).OverlapSize = Ztp2 - Zbt1 +1;
    Nz = handles.AdjustHoriz(n).OverlapSize;
    Nx = handles.StackDim(1,2);
    Ny = handles.StackDim(1,3);
    Joint = zeros(2*Nz, Nx+2*handles.dN, Ny+2*handles.dN);
    Joint(1:Nz,handles.dN+1:handles.dN+1+Nx-1,handles.dN+1:handles.dN+1+Ny-1) = Top(:,:,:);
    Joint(Nz+1:2*Nz,handles.dN+1+Xpos:handles.dN+1+Xpos+Nx-1,handles.dN+1+Ypos:handles.dN+1+Ypos+Ny-1) = Btm(:,:,:);
    handles.AdjustHoriz(n).Joint=Joint;
    axes(handles.axes1);
    imagesc(squeeze(max(handles.AdjustHoriz(n).Joint,[],1))); 
    colormap gray;
    colorbar;
end;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

n = get(handles.popupmenu4,'Value');
Xpos = get(hObject,'Value');
handles.AdjustHoriz(n).dX = Xpos;
set(handles.edit1,'String',num2str(Xpos));

Ypos = get(handles.slider2,'Value');
Top = handles.AdjustHoriz(n).Top;
Btm = handles.AdjustHoriz(n).Btm;
Nz = handles.AdjustHoriz(n).OverlapSize;
Nx = handles.StackDim(1,2);
Ny = handles.StackDim(1,3);
Joint = zeros(2*Nz, Nx+2*handles.dN, Ny+2*handles.dN);
Joint(1:Nz,handles.dN+1:handles.dN+1+Nx-1,handles.dN+1:handles.dN+1+Ny-1) = Top(:,:,:);
Joint(Nz+1:2*Nz,handles.dN+1+Xpos:handles.dN+1+Xpos+Nx-1,handles.dN+1+Ypos:handles.dN+1+Ypos+Ny-1) = Btm(:,:,:);
handles.AdjustHoriz(n).Joint(:,:,:) = Joint(:,:,:);
axes(handles.axes1);
imagesc(squeeze(max(Joint,[],1))); 
colormap gray;
colorbar;
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
n = get(handles.popupmenu4,'Value');
Ypos = get(hObject,'Value');
handles.AdjustHoriz(n).dY = Ypos;
set(handles.edit2,'String',num2str(Ypos));

Xpos = get(handles.slider1,'Value');
Top = handles.AdjustHoriz(n).Top;
Btm = handles.AdjustHoriz(n).Btm;
Nz = handles.AdjustHoriz(n).OverlapSize;
Nx = handles.StackDim(1,2);
Ny = handles.StackDim(1,3);
Joint = zeros(2*Nz, Nx+2*handles.dN, Ny+2*handles.dN);
Joint(1:Nz,handles.dN+1:handles.dN+1+Nx-1,handles.dN+1:handles.dN+1+Ny-1) = Top(:,:,:);
Joint(Nz+1:2*Nz,handles.dN+1+Xpos:handles.dN+1+Xpos+Nx-1,handles.dN+1+Ypos:handles.dN+1+Ypos+Ny-1) = Btm(:,:,:);
handles.AdjustHoriz(n).Joint(:,:,:) = Joint(:,:,:);
axes(handles.axes1);
imagesc(squeeze(max(Joint,[],1))); 
colormap gray;
colorbar;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hfigure = imcontrast(handles.axes1);




function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text9,'BackgroundColor',[1 0 0]);
drawnow;

NewMin = str2double(get(handles.edit3,'String'));
NewMax = str2double(get(handles.edit4,'String'));
n = get(handles.popupmenu2,'Value');
stk = double(handles.Stack{n});
for i=1:handles.StackDim(n,1),
    fr = double(squeeze(stk(i,:,:)));
    fr1 = imadjust(uint16(fr),[NewMin/65535, NewMax/65535],[0 0.5]);
    stk(i,:,:) = fr1(:,:);
end;

contents1 = get(handles.popupmenu3,'Value');
if contents1 == 1,
    Img = squeeze(max(stk,[],1)); % project on X-Y plane
elseif contents1 == 2,
    Img = squeeze(max(stk,[],3)); % project on X-Z plane
else
    Img = squeeze(max(stk,[],2)); % project on Y-Z plane
end;

axes(handles.axes1);
imagesc(Img); 
colormap gray;
colorbar;
NewMin = min(min(Img));
NewMax = max(max(Img));
set(handles.edit3,'String',num2str(NewMin));
set(handles.edit4,'String',num2str(NewMax));

handles.Stack{n} = stk;
guidata(hObject,handles);
set(handles.text9,'BackgroundColor',[0 1 0]);


% Join stacks
% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Nx = handles.StackDim(1,2);
Ny = handles.StackDim(1,3);

N = length(handles.AdjustHoriz);
Tot_dX = 0;
Tot_dY = 0;

for n = 1:N, % adjust 
    stk = handles.Stack{n};
    handles.Stack{n}=[]; %clear memory 
    Nz = handles.StackDim(n,1);
    %EmptyFrame = zeros(Nx+2*handles.dN,Ny+2*handles.dN);
    BigStk = zeros(Nz,Nx+2*handles.dN,Ny+2*handles.dN);
    if n == 1,
        Tot_dX = 0;
        Tot_dY = 0;
    else
        Tot_dX = Tot_dX + handles.AdjustHoriz(n).dX;
        Tot_dY = Tot_dY + handles.AdjustHoriz(n).dY;
    end;
    
    BigStk(1:Nz,handles.dN+1+Tot_dX:handles.dN+Tot_dX+Nx,handles.dN+1+Tot_dY:handles.dN+Tot_dY+Ny) = stk(:,:,:);
    clear stk;
    handles.BigStack{n} = BigStk;
    clear BigStck;
end;

TotNz = handles.StackZplace(N,2)-handles.StackZplace(1,1)+1;
Stack = zeros(TotNz,Nx+2*handles.dN,Ny+2*handles.dN);
BigStk = handles.BigStack{1};
handles.BigStack{1} = []; %clear from memory
Stack(1:handles.StackDim(1,1),:,:) = BigStk(:,:,:);

for n=2:N,
    clear BigStk;
    BigStk = double(handles.BigStack{n});
    handles.BigStack{n} = [];
    Overlap =  handles.AdjustHoriz(n).OverlapSize;
    Stack(handles.StackZplace(n,1):handles.StackZplace(n-1,2),:,:) = (Stack(handles.StackZplace(n,1):handles.StackZplace(n-1,2),:,:) + BigStk(1:Overlap,:,:) )./2;
    Stack(handles.StackZplace(n-1,2)+1:handles.StackZplace(n,2),:,:) = BigStk(Overlap+1:handles.StackDim(n,1),:,:);
end;
    
TotalStack = uint16(Stack);
clear Stack;
TotalStack1 = TotalStack(1:TotNz,handles.dN+1:handles.dN+Nx,handles.dN+1:handles.dN+Ny);
clear TotalStack;
%handles.TotalStack = TotalStack1;
save 'HandlesTotalStack.mat' TotalStack1;
Convert3DmatToTIFF(TotalStack1, 'HandlesTotalStack.tiff');

axes(handles.axes1);
imagesc(squeeze(max(TotalStack1,[],1))); 
colormap gray;
colorbar;

guidata(hObject,handles);



% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
save 'GUI_Handle' handles -V7.3;




% --- Executes on button press in SaveCurrentStackPushButton.
function SaveCurrentStackPushButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveCurrentStackPushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n = get(handles.popupmenu2,'Value');
stk = double(handles.Stack{n});
save 'HandlesTotalSingleStack.mat' stk;
Convert3DmatToTIFF(stk, 'HandlesTotalSingleStack.tiff');
