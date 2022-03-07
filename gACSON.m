%Copyright (c) 2020 Andrea Behanova, Ali Abdollahzadeh, Alejandra Sierra, Jussi Tohka
%A.I. Virtanen Institute for Molecular Sciences, University of Eastern
%Finland, Finland

%version 1.0  11.02.2020

%Permission is hereby granted, free of charge, to any person obtaining a copy 
%of this software and associated documentation files (the "Software"), to deal
%in the Software without restriction, including without limitation the rights
%to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%copies of the Software, and to permit persons to whom the Software is
%furnished to do so, subject to the following conditions:

%The above copyright notice and this permission notice shall be included in all
%copies or substantial portions of the Software.

%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%SOFTWARE.

%This software uses BM4D Matlab software, which can be found in http://www.cs.tut.fi/~foi/GCF-BM3D/
%This software uses Bio-Formats 5.9.2 package, which can be found in https://www.openmicroscopy.org/bio-formats/downloads/

function varargout = gACSON(varargin)
% GACSON MATLAB code for gACSON.fig
%      GACSON, by itself, creates a new GACSON or raises the existing
%      singleton*.
%
%      H = GACSON returns the handle to a new GACSON or the handle to
%      the existing singleton*.
%
%      GACSON('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GACSON.M with the given input arguments.
%
%      GACSON('Property','Value',...) creates a new GACSON or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gACSON_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gACSON_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gACSON

% Last Modified by GUIDE v2.5 11-Feb-2022 16:43:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gACSON_OpeningFcn, ...
                   'gui_OutputFcn',  @gACSON_OutputFcn, ...
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


% --- Executes just before gACSON is made visible.
function gACSON_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gACSON (see VARARGIN)

% Choose default command line output for gACSON
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gACSON wait for user response (see UIRESUME)
% uiwait(handles.figure1);
handles.image = [];
handles.size = [];
handles.vox_size = [];
handles.dispChoice = [];
handles.overChoice = [];
handles.indx = [];
handles.indxLBL = [];
handles.last_path = [];

set(handles.listFiles,'String','File list:');
set(handles.DisplayPopUp,'String','No Display');
set(handles.DisplayPopUp,'Value',1);
set(handles.OverlayPopUp,'String','No Overlay');
set(handles.readImage,'String','Select image...');
set(handles.morph_IAS_lbl,'String','Select IAS label...');
set(handles.Morph_myelin_lbl,'String','Select myelin label...');
set(handles.SlidesText, 'String', 'Slices:')

%morphometry
handles.morph.im = [];
handles.morph.myelin = [];
handles.morph.grid = 15;

%slices
handles.S = 1;
set(handles.sliderSlices, 'Value', 1)

%alpha
handles.alpha = 0.3;
set(handles.alphafactor,'Value',handles.alpha);

%buttons
set(handles.Erase_button, 'Value', 0)

handles.opt = [];
handles.opt.filt_im = 1;
handles.opt.read_mask = {};

txt = get(handles.myelin_indx,'String');
txt = {txt;'Threshold selection'; 'Manual selection'; 'Machine learning'};
set(handles.myelin_indx,'String',txt);
handles.opt.myelin_inx = [];
handles.opt.read_im = [];

txt = get(handles.seed_location,'String');
txt = {txt; 'Watershed centroids'; 'Manual selection'};
set(handles.seed_location,'String',txt);
% handles.opt.seed_loc = 'Regional maxima';
handles.seed_loc = 'Regional maxima';

set (gcf, 'WindowScrollWheelFcn', {@mouseScroll, handles}); 

handles.tp = [];
handles.fp = [];
handles.fn = 0;
handles.tp_idx = [];
handles.fp_idx = [];
handles.fn_idx = [];

handles.sbplt = 0;

handles.res = [];

%add paths

currentFolder = fileparts(which('gACSON.m'));
addpath(genpath(currentFolder));

%path to the Bio Formats
addpath('..\bfmatlab')

%Machine Learning - myelin detection
handles.ML.myelin.pos = [];
handles.ML.bg.pos = [];
handles.ML.train = [];

%subplots
handles.subplots.slider1 = 1;
handles.subplots.slider2 = 1;
handles.subplots.slider3 = 1;

guidata(hObject, handles);


function dispSelection(hObject, eventdata, handles)
%selection change in display pop up menu
if isempty(handles.dispChoice)
    return
end
    
if handles.sbplt == 1
    errordlg('Change the view', 'One image should be displayed');        
    return
end

if strcmp(handles.dispChoice, 'No Display')
    handles.dispChoice = [];
    handles.indx = [];
    set(get(gca,'children'),'Visible','off')
else
    set(handles.subplot, 'Enable', 'on')
    indx = strcmp(handles.listFiles.String, handles.dispChoice);
    indx = indx(2:end);
    handles.indx = indx;
    image = handles.image(indx);
    inx = handles.size{handles.indx};
    if length(inx)<3
        inx = [inx, 1];
    end
    
    if handles.S>inx(3)
        handles.S=inx(3);
    end
    axes(handles.Picture)
    imshow(image{1,1}(:,:,handles.S),[]);
    hold on
    
    set(handles.Picture,'xlim',[0 inx(2)],'ylim',[0 inx(1)])
end

if isempty(handles.indx)
    return
end

inx = handles.size(handles.indx);
if length(inx{1,1})==2
    sno = 1;
else
    sno = inx{1,1}(3);
end
set(handles.SlidesText, 'String', sprintf('Slices: %d/%d', handles.S, sno))

guidata(hObject, handles);


function Disp(hObject, eventdata, handles)
% set(handles.readImage, 'Value', 1)

% if handles.sbplt == 1
%     return
% end



%set(get(gca,'children'),'Visible','off')

dispSelection(hObject, eventdata, handles)

if handles.sbplt == 1
    handles.subplots.slider1 = 1;
    handles.subplots.slider2 = 1;
    handles.subplots.slider3 = 1;
    indx = strcmp(handles.listFiles.String, handles.dispChoice);
    indx = indx(2:end);
    handles.indx = indx;
    subplots(hObject, eventdata, handles, handles.subplots.slider1, handles.subplots.slider2, handles.subplots.slider3)
end
colors = lines;

%Machine learning myelin
if ~isempty(handles.ML.myelin.pos)
    act_pos = handles.ML.myelin.pos(:,3)==handles.S;
    if sum(act_pos)~=0
        plot(handles.ML.myelin.pos(act_pos,1),handles.ML.myelin.pos(act_pos,2),'.','Color',colors(1,:));
    end
end

if ~isempty(handles.ML.bg.pos)
    act_pos = handles.ML.bg.pos(:,3)==handles.S;
    if sum(act_pos)~=0
        plot(handles.ML.bg.pos(act_pos,1),handles.ML.bg.pos(act_pos,2),'.','Color',colors(2,:));
    end
end

%selection change in overlay pop up menu
if strcmp(handles.overChoice, 'No Overlay')
    set(get(gca,'children'),'Visible','off')
    handles.overChoice = [];
    handles.indxLBL = [];
    dispSelection(hObject, eventdata, handles)
elseif isempty(handles.overChoice)
    return
else
    indx = strcmp(handles.listFiles.String, handles.overChoice);
    indx = indx(2:end);
    handles.indxLBL = indx;
    image = handles.image(indx);
    slice = handles.S;
    lbl = image{1,1}(:,:,slice);
    Lrgb = label2rgb(round(lbl),'jet', 'k', 'shuffle');
    %Lrgb = label2rgb(round(lbl),'lines', 'k');
    axes(handles.Picture)
    himage = imshow(Lrgb); himage.AlphaData = handles.alpha;
    handles.himage = himage;
    hold off
end

if ~isempty(handles.tp_idx)
    for i = 1:size(handles.tp_idx,1)
        if handles.tp_idx(i,3)==handles.S
            text(handles.tp_idx(i,2),handles.tp_idx(i,1),'TP','Color','green','FontSize',10);
        end
    end
end

if ~isempty(handles.fp_idx)
    for i = 1:size(handles.fp_idx,1)
        if handles.fp_idx(i,3)==handles.S
            text(handles.fp_idx(i,2),handles.fp_idx(i,1),'FP','Color','red','FontSize',10);
        end
    end
end

if ~isempty(handles.fn_idx)
    for i = 1:size(handles.fn_idx,1)
        if handles.fn_idx(i,3)==handles.S
            text(handles.fn_idx(i,2),handles.fn_idx(i,1),'FN','Color','yellow','FontSize',10);
        end
    end
end


    
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = gACSON_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in DisplayPopUp.
function DisplayPopUp_Callback(hObject, eventdata, handles)
% hObject    handle to DisplayPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns DisplayPopUp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from DisplayPopUp

% selection
handles.MLseg.Nlbl = 2;
handles.MLseg.lbl0 = [];
handles.MLseg.lbl1 = [];
handles.MLseg.lbl2 = [];
handles.MLseg.lbl3 = [];
handles.MLseg.lbl4 = [];
handles.MLseg.train = [];

handles.res = [];

contents = cellstr(get(hObject,'String'));
display_choice = contents(get(hObject, 'Value'));
handles.dispChoice = display_choice;
guidata(hObject, handles);
Disp(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function DisplayPopUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DisplayPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in OverlayPopUp.
function OverlayPopUp_Callback(hObject, eventdata, handles)
% hObject    handle to OverlayPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns OverlayPopUp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from OverlayPopUp

%selection
if handles.sbplt == 1
    errordlg('Change the view','One image should be displayed');
    set(handles.OverlayPopUp, 'Value', 1);
    return
end

contents = cellstr(get(hObject,'String'));
overlay_choice = contents(get(hObject, 'Value'));
handles.overChoice = overlay_choice;

handles.tp = [];
handles.fp = [];
handles.fn = 0;
handles.tp_idx = [];
handles.fp_idx = [];
handles.fn_idx = [];

set(handles.presition, 'String', 'Precision = 0');
set(handles.recall, 'String', 'Recall = 0');
set(handles.F1, 'String', 'F1 score = 0');
set(handles.TP_string, 'String', '= 0');
set(handles.FP_string, 'String', '= 0');
set(handles.FN_string, 'String', '= 0');

guidata(hObject, handles);
Disp(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function OverlayPopUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OverlayPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuOpen_Callback(hObject, eventdata, handles)
% hObject    handle to MenuOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.last_path)
    [Files, pathname] = uigetfile('*.*', 'Select files to load:','MultiSelect','on');
else
    path = handles.last_path;
    [Files, pathname] = uigetfile('*.*', 'Select files to load:',path,'MultiSelect','on');
end

if size(Files,2) == 1
    return
end

handles.last_path = pathname;
guidata(hObject, handles);

if iscell(Files)
    L = size(Files,2);
    f = waitbar(1/L,'Loading files');
else
    L = 1;
    filename = Files;
end

for i = 1:L
    if iscell(Files)
        filename = Files{1,i};
    end
    
    filepath = [pathname filename];

    [pth,name,ext] = fileparts(filepath);

    prompt = {'x:','y:','z:'};
    dlgtitle = ['Voxel size of ', filename];
    definput = {'1','1','1'};
    dims = [ones(size(definput')) ones(size(definput'))*50];
    answer = inputdlg(prompt,dlgtitle,dims,definput,'on');
    
    if isempty(answer)
        answer = {'1';'1';'1'};
    end
    
    voxel_size = [str2double(answer{1,1}), str2double(answer{2,1}), str2double(answer{3,1})];

    if strcmp(ext,'.mat')
%         f = waitbar(0.3,'Loading');
        image = importdata(fullfile(pathname, filename));
%         waitbar(0.99,f,'Loading');
%         close(f)
    elseif strcmp(ext,'.am')
%         f = waitbar(0.3,'Loading');
        [header,image] = LoadData_Amira(filepath);
%         waitbar(0.99,f,'Loading');
%         close(f)
    elseif strcmp(ext,'.h5')
%         f = waitbar(0.3,'Loading');
        info_h5 = h5info(filepath);
        dataset = {info_h5.Datasets.Name};
        
        [indx,tf] = listdlg('PromptString',{'Select dataset...',''},'ListString',dataset);
        
        image = {};
        for i = 1:length(indx)
            dtst = info_h5.Datasets(indx(i),1).Name;
            dataset = ['/', dtst];
            image{i,1} = h5read(filepath,dataset);
        end
        
%         waitbar(0.99,f,'Loading');
%         close(f)
    else
        data = bfopen([pathname filename]);
        im = data{1,1};
        im = im(:,1);
        if size(im,1)>1
            for i = 1:size(im,1)
                image(:,:,i)=im{i,1};
            end
        else
            image=im{1,1};
        end
    end
    
    if iscell(image)
        n = length(image);
    else
        n = 1;
    end
    
    set_image = image;
    fil_name = filename;
    for j = 1:n
        if iscell(set_image)
            image = set_image{j,1};
            filename = [fil_name, '_',info_h5.Datasets(indx(j),1).Name];
        end
        
        image = squeeze(image);
        sz = size(image);

        % Saving features to memory
        handles.image = [handles.image;{image}];
        handles.size = [handles.size;{sz}];
        handles.vox_size = [handles.vox_size; {voxel_size}];

        % Adding file to file list
        txt = get(handles.listFiles,'String');
        txt = [txt; {filename}];
        set(handles.listFiles,'String',txt);

        % Adding file to Pup Up menus
        poptxt = get(handles.DisplayPopUp,'String');
        poptxt = [poptxt; {filename}];
        set(handles.DisplayPopUp,'String',poptxt);

        poptxt2 = get(handles.OverlayPopUp,'String');
        poptxt2 = [poptxt2; {filename}];
        set(handles.OverlayPopUp,'String',poptxt2);

        %read file pop-up menu
        rdfile = get(handles.readImage,'String');
        rdfile = [rdfile; {filename}];
        set(handles.readImage,'String',rdfile);

        %morphometry
        morphfile1 = get(handles.morph_IAS_lbl,'String');
        morphfile1 = [morphfile1; {filename}];
        set(handles.morph_IAS_lbl,'String',morphfile1);

        morphfile2 = get(handles.Morph_myelin_lbl,'String');
        morphfile2 = [morphfile2; {filename}];
        set(handles.Morph_myelin_lbl,'String',morphfile2);
    end

    if iscell(Files)
        waitbar(i/L,f)
    end
    
end

if iscell(Files)
    close(f)
end

guidata(hObject, handles);

% --------------------------------------------------------------------
function import_workspace_Callback(hObject, eventdata, handles)
% hObject    handle to import_workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

list_ws = evalin('base', 'who');
list_ws = [{'Select variable...'}; list_ws];
[index,tf] = listdlg('PromptString',{'matlab workspace variables:',''},'SelectionMode','single','ListString',list_ws);

if isempty(index)
    return
end

v_name = list_ws{index};

if strcmp(v_name,'Select variable...')
    return
end

image = evalin('base', v_name);
filename = v_name;

prompt = {'x:','y:','z:'};
dlgtitle = ['Voxel size of ', filename];
definput = {'1','1','1'};
dims = [ones(size(definput')) ones(size(definput'))*50];
answer = inputdlg(prompt,dlgtitle,dims,definput,'on');
if isempty(answer)
    answer = {'1';'1';'1'};
end
voxel_size = [str2double(answer{1,1}), str2double(answer{2,1}), str2double(answer{3,1})];
handles.vox_size = [handles.vox_size; {voxel_size}];

image = squeeze(image);
sz = size(image);

% Saving features to memory
handles.image = [handles.image;{image}];
handles.size = [handles.size;{sz}];

% Adding file to file list
txt = get(handles.listFiles,'String');
txt = [txt; {filename}];
set(handles.listFiles,'String',txt);

% Adding file to Pup Up menus
poptxt = get(handles.DisplayPopUp,'String');
poptxt = [poptxt; {filename}];
set(handles.DisplayPopUp,'String',poptxt);

poptxt2 = get(handles.OverlayPopUp,'String');
poptxt2 = [poptxt2; {filename}];
set(handles.OverlayPopUp,'String',poptxt2);

%read file pop-up menu
rdfile = get(handles.readImage,'String');
rdfile = [rdfile; {filename}];
set(handles.readImage,'String',rdfile);

%morphometry
morphfile1 = get(handles.morph_IAS_lbl,'String');
morphfile1 = [morphfile1; {filename}];
set(handles.morph_IAS_lbl,'String',morphfile1);

morphfile2 = get(handles.Morph_myelin_lbl,'String');
morphfile2 = [morphfile2; {filename}];
set(handles.Morph_myelin_lbl,'String',morphfile2);

guidata(hObject, handles);


function mouseScroll(hObject, eventdata,handles)

handles = guidata(hObject);

if handles.sbplt == 1
    return
end

if isempty(handles.size) || (isempty(handles.indx) && isempty(handles.indxLBL))
    maxi = 5;
    set(handles.sliderSlices,'Value', 1)
else    
    if isempty(handles.indx)
        inx = handles.size(handles.indxLBL);
    else
        inx = handles.size(handles.indx);
    end
    
    if length(inx{1,1})==2
        sno = 1;
    else
        sno = inx{1,1}(3);
    end
    maxi = sno;
end

set(handles.sliderSlices, 'Max', maxi);
set(handles.sliderSlices, 'Min', 1);
set(handles.sliderSlices, 'SliderStep' , [1/maxi,0.2] );
sliderValue = get(handles.sliderSlices, 'Value');
handles.S = round(sliderValue);
set(handles.sliderSlices, 'Value', handles.S)

indx = strcmp(handles.listFiles.String, handles.dispChoice);
% if isempty(handles.dispChoice)
%     indx = strcmp(handles.listFiles.String, handles.opt.read_im);
% end
    
indx = indx(2:end);
handles.indx = indx;
img = handles.image(indx);

if isempty(img)
    indx = strcmp(handles.listFiles.String, handles.overChoice);
    indx = indx(2:end);
    handles.indx = indx;
    img = handles.image(indx);
    if isempty(img)
        return
    end
end

img = img{1,1};

S = round(get(handles.sliderSlices,'Value'));
sz = size(img);
if length(sz)==2
    return
end

UPDN = eventdata.VerticalScrollCount;
S = S - UPDN;
if (S < 1); S = 1; elseif (S > sz(3)); S = sz(3); end

set(handles.sliderSlices,'Value',S);
set(handles.SlidesText, 'String', sprintf('Slices: %d/%d', S, sz(3)))

handles.S = S;

Disp(hObject, eventdata, handles);


guidata(hObject, handles);


% --- Executes on slider movement.
function sliderSlices_Callback(hObject, eventdata, handles)
% hObject    handle to sliderSlices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
if isempty(handles.size) || (isempty(handles.indx) && isempty(handles.indxLBL))
    maxi = 5;
    set(handles.sliderSlices, 'Value', 1)
else
    if isempty(handles.indx)
        inx = handles.size(handles.indxLBL);
    else
        inx = handles.size(handles.indx);
    end
    
    if length(inx{1,1})==2
        sno = 1;
    else
        sno = inx{1,1}(3);
    end
    maxi = sno;
end

set(handles.sliderSlices, 'Max', maxi);
set(handles.sliderSlices, 'Min', 1);
set(handles.sliderSlices, 'SliderStep' , [0.01, 0.2] );
sliderValue = get(handles.sliderSlices, 'Value');
handles.S = round(sliderValue);
set(handles.sliderSlices, 'Value', handles.S)

set (gcf, 'WindowScrollWheelFcn', {@mouseScroll, handles}); 
Disp(hObject, eventdata, handles);
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function sliderSlices_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderSlices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in merge_button.
function merge_button_Callback(hObject, eventdata, handles)
% hObject    handle to merge_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.sbplt == 1
    errordlg('Change the view','One image should be displayed');
    return
end

if isempty(handles.overChoice)
    return
end

inx = handles.size(handles.indx);

if length(inx{1,1})==2
    sno = 1;
else
    sno = inx{1,1}(3);
end

r = inx{1,1}(1);
c = inx{1,1}(2);
S = handles.S;
lbl = handles.image(handles.indxLBL);

if isempty(lbl)
    f = errordlg('Load the labels','Labels not found');
    return
end

lbl = lbl{1,1};
handles.pixInx = label2idx(lbl);

[y,x] = getpts; y = round(y(:,1)); x = round(x(:,1));
lx = length(x);
linearInd = zeros(lx,1);

for i = 1:lx
    if (x(i)>0) && (x(i)<r) && (y(i)>0) && (y(i)<c)
        linearInd(i) = sub2ind([r,c,sno], x(i), y(i), S);
    end
end
linearInd = nonzeros(linearInd);
l = nonzeros(lbl(linearInd));
t_ref = l(1);
t = unique(l);

ref_inx = handles.pixInx{t_ref};
for j = 1:length(t)
    if t(j)~=t_ref
        ref_inx = [ref_inx;handles.pixInx{t(j)}];
        lbl(handles.pixInx{t(j)}) = t_ref;
        handles.pixInx{t(j)} = zeros(0,1);
    end
end
handles.pixInx(t_ref) = {ref_inx};

Lrgb = label2rgb(lbl(:,:,S),'jet', 'k', 'shuffle'); hold on
%set(handles.himage, 'cdata', Lrgb ,'AlphaData', handles.alpha); hold off

%saving handles
lbl2{1,1} = lbl;
handles.image(handles.indxLBL) = lbl2;
Disp(hObject, eventdata, handles);
guidata(hObject, handles);


% --- Executes on button press in Erase_button.
function Erase_button_Callback(hObject, eventdata, handles)
% hObject    handle to Erase_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.sbplt == 1
    errordlg('Change the view','One image should be displayed');
    set(handles.Erase_button, 'Value', 0)
    return
end

if isempty(handles.overChoice)
    set(handles.Erase_button, 'Value', 0)
    return
end

inx = handles.size(handles.indx);
if length(inx{1,1})==2
    sno = 1;
else
    sno = inx{1,1}(3);
end
r = inx{1,1}(1);
c = inx{1,1}(2);
S = handles.S;
lbl = handles.image(handles.indxLBL);

if isempty(lbl)
    f = errordlg('Load the labels','Labels not found');
    set(handles.Erase_button, 'Value', 0)
    return
end

lbl = lbl{1,1};

val = get(hObject,'Value');

if val == 0
    Disp(hObject, eventdata, handles);
end

while val
    g = get(gca,'children');
    h = imfreehand();
    if isempty(h)
        return
    end

    BW = createMask(h, g(1));
    delete(h)
    [x,y] = find(BW);
    linearInd = sub2ind([r,c,sno], x, y, repmat(S, length(x),1));
    lbl(linearInd) = 0;
    Lrgb = label2rgb(lbl(:,:,S),'jet', 'k', 'shuffle'); hold on
    %set(handles.himage, 'cdata', Lrgb, 'AlphaData', handles.alpha); hold off
    
    %saving handles
    lbl2{1,1} = lbl;
    handles.image(handles.indxLBL) = lbl2;
    Disp(hObject, eventdata, handles);
end
guidata(hObject, handles);


% --- Executes on button press in Fill_button.
function Fill_button_Callback(hObject, eventdata, handles)
% hObject    handle to Fill_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.sbplt == 1
    errordlg('Change the view','One image should be displayed');
    set(handles.Fill_button, 'Value', 0)
    return
end

if isempty(handles.overChoice)
    set(handles.Fill_button, 'Value', 0)
    return
end

inx = handles.size(handles.indxLBL);

if length(inx{1,1})==2
    sno = 1;
else
    sno = inx{1,1}(3);
end
r = inx{1,1}(1);
c = inx{1,1}(2);
S = handles.S;
lbl = handles.image(handles.indxLBL);

if isempty(lbl)
    f = errordlg('Load the labels','Labels not found');
    set(handles.Fill_button, 'Value', 0)
    return
end

lbl = lbl{1,1};

val = get(hObject,'Value');

if val == 0
    Disp(hObject, eventdata, handles);
end

while val
    g = get(gca,'children');
    h = imfreehand();
    if isempty(h)
        return
    end
    BW = createMask(h, g(1));
    delete(h)
    [x,y] = find(BW);
    linearInd = sub2ind([r,c,sno], x, y, repmat(S, length(x),1));
    lbl(linearInd) = max(lbl(:))+1;
    Lrgb = label2rgb(lbl(:,:,S),'jet', 'k', 'shuffle'); hold on
    %set(handles.himage, 'cdata', Lrgb, 'AlphaData', handles.alpha); hold off
    
    %saving handles
    lbl2{1,1} = lbl;
    handles.image(handles.indxLBL) = lbl2;
    Disp(hObject, eventdata, handles);
end
guidata(hObject, handles);

function mouseScrollAlpha(hObject, eventdata,handles)
handles = guidata(hObject);

S  = get(handles.alphafactor,'Value');

UPDN = eventdata.VerticalScrollCount;
S = S - UPDN*0.05;
if (S < 0.05); S = 0.05; elseif (S > 1); S = 1; end

set(handles.alphafactor,'Value',S);

handles.alpha = S;
guidata(hObject, handles);
Disp(hObject, eventdata, handles);

% --- Executes on slider movement.
function alphafactor_Callback(hObject, eventdata, handles)
% hObject    handle to alphafactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
if handles.sbplt == 1
    errordlg('Change the view','One image should be displayed');
    return
end

set(handles.alphafactor, 'Max',1);
set(handles.alphafactor, 'Min', 0);
set(handles.alphafactor, 'SliderStep' , [0.05,0.05] );
sliderValue = get(handles.alphafactor,'Value');
handles.alpha = sliderValue;
set(handles.alphafactor,'Value', handles.alpha)

set (gcf, 'WindowScrollWheelFcn', {@mouseScrollAlpha, handles});

guidata(hObject, handles);
Disp(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function alphafactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alphafactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set(gcf,'units','normalized','outerpos',[0 0 1 1]);
step = 1/15;

if handles.sbplt
    set(handles.Picture, 'Position', [0.185 0.67 0.38 0.3]);
else
    set(handles.Picture, 'Position', [0.165 0.315 0.83 0.705]);
end

set(handles.SlidesText, 'Position', [0.165 0.285 0.83 0.03]);
set(handles.sliderSlices, 'Position', [0.165 0.26 0.83 0.025]);

set(handles.listFiles, 'Position', [0.02 0.55 1/7 0.43]);
set(handles.uibuttongroup1, 'Position', [0.02 0.16+3*step 1/7 1/15]);
set(handles.uipanel1, 'Position', [0.02 0.10+5*step 1/7 1/9]);
set(handles.text4, 'Position', [0.01 0.05 0.3 0.3]);
set(handles.text3, 'Position', [0.01 0.5 0.3 0.3]);
set(handles.OverlayPopUp, 'Position', [0.33 0 0.6 0.4]);
set(handles.DisplayPopUp, 'Position', [0.33 0.45 0.6 0.4]);
set(handles.alphafactor, 'Position', [0.05 0.3 0.9 0.4]);

set(handles.manual_segm,'Units', 'normalized');
set(handles.manual_segm, 'Position', [0.02 0.16+2*step 1/7 1/15]);
set(handles.Erase_button, 'Position', [0.045 0.2 0.28 0.6]);
set(handles.Fill_button, 'Position', [0.365 0.2 1/4 0.6]);
set(handles.merge_button, 'Position', [0.655 0.2 0.3 0.6]);

set(handles.Iso_surf, 'Position', [0.02 0.26 0.035 0.03]);
set(handles.Delete,'Units', 'normalized');
set(handles.Delete, 'Position', [0.056 0.26 0.055 0.03]);
set(handles.save_label,'Units', 'normalized');
set(handles.save_label, 'Position', [0.112 0.26 0.05 0.03]);

set(handles.text25,'Units', 'normalized');
set(handles.text25,'Position', [0.19 0.01 0.8 0.03]);
set(handles.percentage,'Units', 'normalized');
set(handles.percentage,'Position', [0.19 0.01 0 0.03]);
set(handles.perc_string,'Units', 'normalized');
set(handles.perc_string,'Position', [0.55 0.033 0.05 0.03]);


%Pipeline part
set(handles.pipeline,'Units', 'normalized');
set(handles.pipeline,'Position', [0.02 0.02 0.17 0.235]);
set(handles.readImage,'Units', 'normalized');
set(handles.readImage,'Position', [0.01 0.89 0.5 0.1]);
set(handles.uipanel2,'Units', 'normalized');
set(handles.uipanel2,'Position', [0.01 0.61 0.5 0.25]);
set(handles.myelin_indx,'Units', 'normalized');
set(handles.myelin_indx,'Position', [0.05 0.5 0.9 0.4]);

set(handles.uipanel3,'Units', 'normalized');
set(handles.uipanel3,'Position', [0.01 0.32 0.5 0.3]);
set(handles.min_myelin_vol,'Units', 'normalized');
set(handles.min_myelin_vol,'Position', [0 0.65 0.5 0.33]);
set(handles.max_myelin_vol,'Units', 'normalized');
set(handles.max_myelin_vol,'Position', [0 0.2 0.5 0.33]);
set(handles.min_myel_vol,'Units', 'normalized');
set(handles.min_myel_vol,'Position', [0.55 0.64 0.3 0.36]);
set(handles.max_myel_vol,'Units', 'normalized');
set(handles.max_myel_vol,'Position', [0.55 0.16 0.3 0.38]);

set(handles.uipanel4,'Units', 'normalized');
set(handles.uipanel4,'Position', [0.01 0.01 0.5 0.32]);
set(handles.text12,'Units', 'normalized');
set(handles.text12,'Position', [0 0.65 0.5 0.33]);
set(handles.text13,'Units', 'normalized');
set(handles.text13,'Position', [0 0.2 0.5 0.33]);
set(handles.min_lbl_vol,'Units', 'normalized');
set(handles.min_lbl_vol,'Position', [0.55 0.64 0.3 0.36]);
set(handles.max_lbl_vol,'Units', 'normalized');
set(handles.max_lbl_vol,'Position', [0.55 0.16 0.3 0.38]);

set(handles.text8,'Units', 'normalized');
set(handles.text8,'Position', [0.51 0.81 0.3 0.17]);
set(handles.similarity_th,'Units', 'normalized');
set(handles.similarity_th,'Position', [0.78 0.86 0.13 0.1]);
set(handles.text18,'Units', 'normalized');
set(handles.text18,'Position', [0.51 0.73 0.5 0.07]);
set(handles.Small_vol_th,'Units', 'normalized');
set(handles.Small_vol_th,'FontSize', 9);
set(handles.Small_vol_th,'Position', [0.86 0.72 0.13 0.1]);

set(handles.mask,'Units', 'normalized');
set(handles.mask,'Position', [0.51 0.6 0.5 0.15]);
set(handles.Filter,'Units', 'normalized');
set(handles.Filter,'Position', [0.51 0.48 0.5 0.15]);
set(handles.seed_loc_text,'Units', 'normalized');
set(handles.seed_loc_text,'Position', [0.51 0.26 0.49 0.25]);
set(handles.seed_location,'Units', 'normalized');
set(handles.seed_location,'Position', [0.05 0.05 0.9 0.8]);


set(handles.buttonSegmentation,'Units', 'normalized');
set(handles.buttonSegmentation,'Position', [0.52 0.03 0.46 0.2]);
% set(handles.choice_segmen,'Units', 'normalized');
% set(handles.choice_segmen,'Position', [0.75 0.35 0.2 0.2]);
% set(handles.new_segm,'Units', 'normalized');
% set(handles.new_segm,'Position', [0.75 0.1 0.2 0.2]);


%Evaluation
set(handles.Categorization,'Units', 'normalized');
set(handles.Categorization,'Position', [0.9 0.02 0.09 0.235]);
set(handles.TP,'Units', 'normalized');
set(handles.TP,'Position', [0.05 0.85 0.25 0.12]);
set(handles.FP,'Units', 'normalized');
set(handles.FP,'Position', [0.05 0.7 0.25 0.12]);
set(handles.FN,'Units', 'normalized');
set(handles.FN,'Position', [0.05 0.55 0.25 0.12]);
set(handles.TP_string,'Units', 'normalized');
set(handles.TP_string,'Position', [0.30 0.83 0.4 0.12]);
set(handles.FP_string,'Units', 'normalized');
set(handles.FP_string,'Position', [0.3 0.68 0.4 0.12]);
set(handles.FN_string,'Units', 'normalized');
set(handles.FN_string,'Position', [0.3 0.53 0.4 0.12]);
set(handles.presition,'Units', 'normalized');
set(handles.presition,'Position', [0.05 0.4 0.9 0.12]);
set(handles.recall,'Units', 'normalized');
set(handles.recall,'Position', [0.05 0.3 0.9 0.12]);
set(handles.F1,'Units', 'normalized');
set(handles.F1,'Position', [0.05 0.2 0.9 0.12]);
set(handles.save,'Units', 'normalized');
set(handles.save,'Position', [0.35 0.05 0.4 0.15]);

%Morphometry
set(handles.morphometry_panel,'Units', 'normalized');
set(handles.morphometry_panel,'Position', [0.193 0.09 0.125 0.165]);
set(handles.morph_IAS_lbl,'Units', 'normalized');
set(handles.morph_IAS_lbl,'Position', [0.05 0.95 0.9 0.05]);
set(handles.Morph_myelin_lbl,'Units', 'normalized');
set(handles.Morph_myelin_lbl,'Position', [0.05 0.72 0.9 0.05]);
set(handles.grid_size_text,'Units', 'normalized');
set(handles.grid_size_text,'Position', [0.05 0.4 0.5 0.1]);
set(handles.grid_size,'Units', 'normalized');
set(handles.grid_size,'Position', [0.45 0.4 0.2 0.1]);
set(handles.morphometry_button,'Units', 'normalized');
set(handles.morphometry_button,'Position', [0.3 0.02 0.5 0.2]);

%%Machine learning
% colors = lines;
% set(handles.ML_panel,'Units', 'normalized');
% set(handles.ML_panel,'Position', [0.32 0.09 0.15 0.165]);
% set(handles.label0_txt,'Units', 'normalized');
% set(handles.label0_txt,'Position', [0.01 0.86 0.3 0.15]);
% set(handles.label0,'Units', 'normalized');
% set(handles.label0,'Position', [0.5 0.86 0.2 0.15]);
% set(handles.label0,'BackgroundColor',colors(1,:)); 
% 
% set(handles.label1_txt,'Units', 'normalized');
% set(handles.label1_txt,'Position', [0.01 0.7 0.3 0.15]);
% set(handles.label1,'Units', 'normalized');
% set(handles.label1,'Position', [0.5 0.7 0.2 0.15]);
% set(handles.label1,'BackgroundColor',colors(2,:)); 
% 
% set(handles.label2_txt,'Units', 'normalized');
% set(handles.label2_txt,'Position', [0.01 0.54 0.3 0.15]);
% set(handles.label2,'Units', 'normalized');
% set(handles.label2,'Position', [0.5 0.54 0.2 0.15]);
% set(handles.label2,'BackgroundColor',colors(3,:));
% set(handles.label2_txt,'Visible','off'); 
% set(handles.label2,'Visible','off'); 
% 
% set(handles.label3_txt,'Units', 'normalized');
% set(handles.label3_txt,'Position', [0.01 0.38 0.3 0.15]);
% set(handles.label3,'Units', 'normalized');
% set(handles.label3,'Position', [0.5 0.38 0.2 0.15]);
% set(handles.label3,'BackgroundColor',colors(4,:));
% set(handles.label3_txt,'Visible','off'); 
% set(handles.label3,'Visible','off');
% 
% set(handles.label4_txt,'Units', 'normalized');
% set(handles.label4_txt,'Position', [0.01 0.22 0.3 0.15]);
% set(handles.label4,'Units', 'normalized');
% set(handles.label4,'Position', [0.5 0.22 0.2 0.15]);
% set(handles.label4,'BackgroundColor',colors(5,:));
% set(handles.label4_txt,'Visible','off'); 
% set(handles.label4,'Visible','off');
% 
% set(handles.ML_addlabel,'Units', 'normalized');
% set(handles.ML_addlabel,'Position', [0.01 0.01 0.3 0.2]);
% set(handles.reset_segm,'Units', 'normalized');
% set(handles.reset_segm,'Position', [0.7 0.01 0.3 0.2]);
% set(handles.MLseg_done,'Units', 'normalized');
% set(handles.MLseg_done,'Position', [0.35 0.01 0.3 0.2]);


% --- Executes on button press in Iso_surf.
function Iso_surf_Callback(hObject, eventdata, handles)
% hObject    handle to Iso_surf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.sbplt == 1
    errordlg('Change the view','One image should be displayed');
    return
end

if isempty(handles.overChoice)
    return
end

inx = handles.size(handles.indx);

if length(inx{1,1})==2
    sno = 1;
else
    sno = inx{1,1}(3);
end
r = inx{1,1}(1);
c = inx{1,1}(2);
S = handles.S;
lbl = handles.image(handles.indxLBL);

if isempty(lbl)
    f = errordlg('Load the labels','Labels not found');
    return
end
[y,x] = ginputc(1); y = round(y(:,1)); x = round(x(:,1));

lbl = lbl{1,1};
handles.pixInx = label2idx(lbl);
stats = regionprops(lbl,'BoundingBox');

if (x>0) && (x<r) && (y>0) && (y<c)
    linearInd = sub2ind([r,c,sno], x, y, S);            
    l = nonzeros(lbl(linearInd));            
    bw = false(r,c,sno);
    bw(handles.pixInx{l}) = true;           
    bb = stats(l).BoundingBox;
    [bbw,~] = util_extract_bounded_obj(bw,bb,[0,0,0]);
    figure; isosurface(bbw,0.5); axis equal
end

% --- Executes on button press in Delete.
function Delete_Callback(hObject, eventdata, handles)
% hObject    handle to Delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.sbplt == 1
    errordlg('Change the view','One image should be displayed');
    return
end

if isempty(handles.overChoice)
    return
end

inx = handles.size(handles.indxLBL);

if length(inx{1,1})==2
    sno = 1;
else
    sno = inx{1,1}(3);
end
r = inx{1,1}(1);
c = inx{1,1}(2);
S = handles.S;
lbl = handles.image(handles.indxLBL);

if isempty(lbl)
    f = errordlg('Load the labels','Labels not found');
    return
end

lbl = lbl{1,1};

[y,x] = ginput(1); y = round(y(:,1)); x = round(x(:,1));
if (x>0) && (x<r) && (y>0) && (y<c)
    lbl_idx = lbl(x,y,handles.S);
    lbl(lbl==lbl_idx) = 0;
    handles.image{handles.indxLBL} = lbl;
end
Disp(hObject, eventdata, handles);
guidata(hObject, handles);

% --- Executes on button press in save_label.
function save_label_Callback(hObject, eventdata, handles)
% hObject    handle to save_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.overChoice)
    return
end

inx = handles.size(handles.indxLBL);

if length(inx{1,1})==2
    sno = 1;
else
    sno = inx{1,1}(3);
end
r = inx{1,1}(1);
c = inx{1,1}(2);
S = handles.S;
lbl = handles.image(handles.indxLBL);

if isempty(lbl)
    f = errordlg('Load the labels','Labels not found');
    return
end

lbl = lbl{1,1};

[filename, filepath] = uiputfile('*.mat', 'Save the label to a file:',handles.overChoice{1,1});
if logical(filename)
    FileName = fullfile(filepath, filename);
    save(FileName, 'lbl', '-v7.3');
end


% --- Executes on selection change in readImage.
function readImage_Callback(hObject, eventdata, handles)
% hObject    handle to readImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns readImage contents as cell array
%        contents{get(hObject,'Value')} returns selected item from readImage
if handles.sbplt == 1
    set(handles.subplot,'State','off')
    set(handles.SlidesText,'Visible','on')
    set(handles.sliderSlices,'Visible','on')
    set (gcf, 'WindowScrollWheelFcn', {@mouseScroll, handles}); 

    handles.sbplt = 0;

    set(handles.Picture, 'Position', [0.165 0.315 0.83 0.705]);
    set(handles.hs2,'Visible','off')
    set(handles.hs3,'Visible','off')
    set(handles.hs4,'Visible','off')
    set(handles.subplot4,'Visible','off')
    set(handles.slider_sbplt1,'Visible','off')
    set(handles.slider_sbplt2,'Visible','off')
    set(handles.slider_sbplt3,'Visible','off')
    set(get(gca,'children'),'Visible','off')
end


contents = cellstr(get(hObject,'String'));
readIM_choice = contents(get(hObject, 'Value'));
handles.opt.read_im = readIM_choice{1,1};

handles.res = [];

if strcmp(readIM_choice, 'Select image...')
    set(get(gca,'children'),'Visible','off')
else
    indx = strcmp(handles.listFiles.String, readIM_choice);
    indx = indx(2:end);
    handles.indx = indx;
    image = handles.image(indx);
    inx = handles.size{handles.indx};
    if length(inx)<3
        inx = [inx, 1];
    end
    
    if handles.S>inx(3)
        handles.S=inx(3);
    end
    
    axes(handles.Picture)
    imshow(image{1,1}(:,:,handles.S),[]);
    
    set(handles.Picture,'xlim',[0 inx(2)],'ylim',[0 inx(1)])
end

if isempty(handles.indx)
    return
end

inx = handles.size(handles.indx);
if length(inx{1,1})==2
    sno = 1;
else
    sno = inx{1,1}(3);
end
set(handles.SlidesText, 'String', sprintf('Slices: %d/%d', handles.S, sno))

handles.dispChoice = handles.opt.read_im;
handles.overChoice = [];

set(handles.DisplayPopUp, 'Value', 1)
set(handles.OverlayPopUp, 'Value', 1)

%Machine Learning myelin
handles.ML.myelin.pos = [];
handles.ML.bg.pos = [];
handles.ML.train = [];

%machine learning segmentation
handles.MLseg.lbl0 = [];
handles.MLseg.lbl1 = [];
handles.MLseg.lbl2 = [];
handles.MLseg.lbl3 = [];
handles.MLseg.lbl4 = [];
handles.MLseg.train = [];

% Disp(hObject, eventdata, handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function readImage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to readImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in mask.
function mask_Callback(hObject, eventdata, handles)
% hObject    handle to mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mask
mask = get(handles.mask,'value');

if mask == 0
    handles.opt.read_mask = {};
else
    handles.opt.read_mask = uigetfile;
end

guidata(hObject, handles);



function similarity_th_Callback(hObject, eventdata, handles)
% hObject    handle to similarity_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of similarity_th as text
%        str2double(get(hObject,'String')) returns contents of similarity_th as a double
sim_th = get(handles.similarity_th,'String');
handles.similarity_th.Value = str2double(sim_th);
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function similarity_th_CreateFcn(hObject, eventdata, handles)
% hObject    handle to similarity_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function min_myel_vol_Callback(hObject, eventdata, handles)
% hObject    handle to min_myel_vol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_myel_vol as text
%        str2double(get(hObject,'String')) returns contents of min_myel_vol as a double
min_myel = get(handles.min_myel_vol,'String');
handles.min_myel_vol.Value = str2double(min_myel);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function min_myel_vol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_myel_vol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max_myel_vol_Callback(hObject, eventdata, handles)
% hObject    handle to max_myel_vol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_myel_vol as text
%        str2double(get(hObject,'String')) returns contents of max_myel_vol as a double
max_myel = get(handles.max_myel_vol,'String');
handles.max_myel_vol.Value = str2double(max_myel);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function max_myel_vol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_myel_vol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in myelin_indx.
function myelin_indx_Callback(hObject, eventdata, handles)
% hObject    handle to myelin_indx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns myelin_indx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from myelin_indx


if isempty(handles.opt.read_im)
    f = errordlg('Load the image','Image not found');
    return
end

indx = strcmp(handles.listFiles.String, handles.opt.read_im);
indx = indx(2:end);
handles.indx = indx;
image = handles.image(indx);
image3D = image{1,1};
image = image{1,1}(:,:,handles.S);
set_low_intensity = [];
handles.opt.myelin_inx = [];

contents = cellstr(get(hObject,'String'));
choice = contents(get(hObject, 'Value'));
if strcmp(choice, 'Default - random')
    for i = 1:size(image,1)
        for j = 1:size(image,2)
            if image(i,j)<90
                set_low_intensity = [set_low_intensity; i, j];
            end
        end
    end 

    random_index = randsample(length(set_low_intensity),100);
    handles.opt.myelin_inx = [set_low_intensity(random_index,:), (handles.S+1)*ones(1,length(random_index))'];
elseif strcmp(choice, 'Threshold selection')
    pic = image;
    myelin_detection(pic,handles)
    handles.opt.myelin_inx = 0;
    
   
elseif strcmp(choice, 'Manual selection')
    [y,x] = getpts; y = round(y(:,1)); x = round(x(:,1));
    z = (handles.S+1)*ones(1,length(x))';
    handles.opt.myelin_inx = [x, y, z];

elseif strcmp(choice, 'Machine learning') 
    set(handles.mach_learn,'Visible','on')
    set(handles.mach_learn,'Units', 'normalized');
    set(handles.mach_learn, 'Position', [0.9 0.85 0.1 0.15]);
    set(handles.myelin_ML,'Visible','on')
    set(handles.myelin_ML,'Units', 'normalized');
    set(handles.myelin_ML, 'Position', [0 0.5 0.5 0.5]);
    set(handles.background_ML,'Visible','on')
    set(handles.background_ML,'Units', 'normalized');
    set(handles.background_ML, 'Position', [0.5 0.5 0.5 0.5]);
    set(handles.ML_done,'Visible','on')
    set(handles.ML_done,'Units', 'normalized');
    set(handles.ML_done, 'Position', [0 0 1 0.5]);
    
    handles.opt.myelin_inx = 1;
%     g = get(gca,'children');
%  h = imfreehand;
%     h = draw('LineWidth',5,'Closed',0,'HandleVisibility','off', 'FaceAlpha', 0,'Deletable',false, 'InteractionsAllowed', 'none');
%     [myobj,xs,ys] = freehanddraw(gca,'color','r','linewidth',5,'closed',0);
    
    
%     i=1;
% 
% 
% 
%     while i == 1
%         answer = questdlg('Do you want to add more?', ...
%         'More data', ...
%         'Continue','Finish','Continue');
%         switch answer
%             case 'Continue'
%                 draw(hObject,handles,image,image3D)
%             case 'Finish'
%                 i=2;
%         end
%     end

end
    
    
    
%     pos = round(h.Position);
%     pos2 = sub2ind(size(image),pos(:,2),pos(:,1));
%     intensities = image(pos2);
    
% Disp(hObject, eventdata, handles)

guidata(hObject, handles);


    



% --- Executes during object creation, after setting all properties.
function myelin_indx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to myelin_indx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function min_lbl_vol_Callback(hObject, eventdata, handles)
% hObject    handle to min_lbl_vol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_lbl_vol as text
%        str2double(get(hObject,'String')) returns contents of min_lbl_vol as a double
min_lbl = get(handles.min_lbl_vol,'String');
handles.min_lbl_vol.Value = str2double(min_lbl);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function min_lbl_vol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_lbl_vol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max_lbl_vol_Callback(hObject, eventdata, handles)
% hObject    handle to max_lbl_vol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_lbl_vol as text
%        str2double(get(hObject,'String')) returns contents of max_lbl_vol as a double
max_lbl = get(handles.max_lbl_vol,'String');
handles.max_lbl_vol.Value = str2double(max_lbl);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function max_lbl_vol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_lbl_vol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Small_vol_th_Callback(hObject, eventdata, handles)
% hObject    handle to Small_vol_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Small_vol_th as text
%        str2double(get(hObject,'String')) returns contents of Small_vol_th as a double
vol_th = get(handles.Small_vol_th,'String');
handles.Small_vol_th.Value = str2double(vol_th);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Small_vol_th_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Small_vol_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Filter.
function Filter_Callback(hObject, eventdata, handles)
% hObject    handle to Filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Filter
 filter = get(handles.Filter,'value');

if filter == 0
    handles.opt.filt_im = 1;
else
    handles.opt.filt_im = 0;
    
        answer = questdlg('Which profile do you want to use?', ...
	'Filter profile', ...
	'Low complexity profile','Normal profile','Modified profile','Low complexity profile');
    switch answer
        case 'Low complexity profile'
            handles.opt.filt_prof = 'lc';
        case 'Normal profile'
            handles.opt.filt_prof = 'np';
        case 'Modified profile'
            handles.opt.filt_prof = 'mp';
    end
end

guidata(hObject, handles);


% --- Executes on button press in buttonSegmentation.
function buttonSegmentation_Callback(hObject, eventdata, handles)
% hObject    handle to buttonSegmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.opt.read_im)
    f = errordlg('Load the image','Image not found');
    return
end

choice = handles.seed_loc;
if strcmp(choice, 'Regional maxima')
    handles.opt.similarity_th = get(handles.similarity_th,'value');
    handles.opt.min_max_myelin_volume = [get(handles.min_myel_vol,'value'), ...
        get(handles.max_myel_vol,'value');];
    handles.opt.min_max_lbl_volume = [get(handles.min_lbl_vol,'value'), ...
        get(handles.max_lbl_vol,'value');];
    handles.opt.closing_small_volume_threshold = get(handles.Small_vol_th,'value');
    handles.opt.foreground = 1;
    handles.opt.background = 0;



    if isempty(handles.opt.myelin_inx)
        indx = strcmp(handles.listFiles.String, handles.opt.read_im);
        indx = indx(2:end);
        handles.indx = indx;
        image = handles.image(indx);
        image = image{1,1}(:,:,handles.S);
        set_low_intensity = [];

        for i = 1:size(image,1)
            for j = 1:size(image,2)
                if image(i,j)<70
                    set_low_intensity = [set_low_intensity; i, j];
                end
            end
        end

        random_index = randsample(length(set_low_intensity),100);
        handles.opt.myelin_inx = [set_low_intensity(random_index,:), (handles.S+1)*ones(1,length(random_index))'];
    end

    opt = handles.opt;

    guidata(hObject, handles);

    util_HM_SBEM_segmentation_pipeline(opt,hObject,handles);

    

elseif strcmp(choice, 'Watershed centroids')
    handles.opt.similarity_th = get(handles.similarity_th,'value');
    handles.opt.min_max_myelin_volume = [get(handles.min_myel_vol,'value'), ...
        get(handles.max_myel_vol,'value');];
    handles.opt.min_max_lbl_volume = [get(handles.min_lbl_vol,'value'), ...
        get(handles.max_lbl_vol,'value');];
    handles.opt.closing_small_volume_threshold = get(handles.Small_vol_th,'value');
    handles.opt.foreground = 1;
    handles.opt.background = 0;

    if isempty(handles.opt.myelin_inx)
        indx = strcmp(handles.listFiles.String, handles.opt.read_im);
        indx = indx(2:end);
        handles.indx = indx;
        image = handles.image(indx);
        image = image{1,1}(:,:,handles.S);
        set_low_intensity = [];

        for i = 1:size(image,1)
            for j = 1:size(image,2)
                if image(i,j)<70
                    set_low_intensity = [set_low_intensity; i, j];
                end
            end
        end

        random_index = randsample(length(set_low_intensity),100);
        handles.opt.myelin_inx = [set_low_intensity(random_index,:), (handles.S+1)*ones(1,length(random_index))'];
    end

    opt = handles.opt;

    guidata(hObject, handles);

    new_segmentation_pipeline(opt,hObject,handles);
elseif strcmp(choice, 'Manual selection')
    handles.opt.similarity_th = get(handles.similarity_th,'value');
    handles.opt.min_max_myelin_volume = [get(handles.min_myel_vol,'value'), ...
        get(handles.max_myel_vol,'value');];
    handles.opt.min_max_lbl_volume = [get(handles.min_lbl_vol,'value'), ...
        get(handles.max_lbl_vol,'value');];
    handles.opt.closing_small_volume_threshold = get(handles.Small_vol_th,'value');
    handles.opt.foreground = 1;
    handles.opt.background = 0;

    if isempty(handles.opt.myelin_inx)
        indx = strcmp(handles.listFiles.String, handles.opt.read_im);
        indx = indx(2:end);
        handles.indx = indx;
        image = handles.image(indx);
        image = image{1,1}(:,:,handles.S);
        set_low_intensity = [];

        low_inten = find(image<70);
        [x,y] = ind2sub(size(image),low_inten);
        set_low_intensity = [x,y];

        random_index = randsample(length(set_low_intensity),100);
        handles.opt.myelin_inx = [set_low_intensity(random_index,:), (handles.S+1)*ones(1,length(random_index))'];
    end

    opt = handles.opt;

    guidata(hObject, handles);

    label = Choosed_seed_segmenttaion(opt,hObject,handles);
    n_labels=max(label(:));
    fff = msgbox(['Segmentation of manually selected seeds finished: ',num2str(n_labels),' labels detected']);
    
    image = squeeze(label);
    sz = size(image);
    voxel_size = handles.vox_size(indx);
    voxel_size = voxel_size{1,1};
    filename = ['seg_', opt.read_im];

    % Saving features to memory
    handles.image = [handles.image;{image}];
    handles.size = [handles.size;{sz}];
    handles.vox_size = [handles.vox_size; {voxel_size}];

    % Adding file to file list
    txt = get(handles.listFiles,'String');
    txt = [txt; {filename}];
    set(handles.listFiles,'String',txt);

    % Adding file to Pup Up menus
    poptxt = get(handles.DisplayPopUp,'String');
    poptxt = [poptxt; {filename}];
    set(handles.DisplayPopUp,'String',poptxt);

    poptxt2 = get(handles.OverlayPopUp,'String');
    poptxt2 = [poptxt2; {filename}];
    set(handles.OverlayPopUp,'String',poptxt2);

    %read file pop-up menu
    rdfile = get(handles.readImage,'String');
    rdfile = [rdfile; {filename}];
    set(handles.readImage,'String',rdfile);

    %morphometry
    morphfile1 = get(handles.morph_IAS_lbl,'String');
    morphfile1 = [morphfile1; {filename}];
    set(handles.morph_IAS_lbl,'String',morphfile1);

    morphfile2 = get(handles.Morph_myelin_lbl,'String');
    morphfile2 = [morphfile2; {filename}];
    set(handles.Morph_myelin_lbl,'String',morphfile2);
    
%     if length(size(label))==3
%         k = figure; isosurface(label,0.5); axis equal
%     else
%         k = figure; imshow(label,[]);
%     end
%     
%     answer = questdlg('Would you like to add this label to selected label in Overlay?', ...
% 	'Label', ...
% 	'Yes','No','Save this label','Yes');
% 
%     switch answer
%         case 'Yes'
%             
%             idx = find(label);
%             lbl = handles.image{handles.indxLBL};
%             max_lbl = max(lbl(:));
%             lbl(idx) = max_lbl+1;
%             close(k)
%             handles.image{handles.indxLBL} = lbl;
%             Disp(hObject, eventdata, handles);
%         case 'No'
%             close(k)
%         case 'Save this label'
%             [filename, filepath] = uiputfile('*.mat', 'Save the single label in a file:','Manually_segmented_lbl');
%             FileName = fullfile(filepath, filename);
%             CC = bwconncomp(label);
%             save(FileName, 'CC', '-v7.3');
%             close(k)
%     end
    
end

guidata(hObject, handles);


% --- Executes on button press in TP.
function TP_Callback(hObject, eventdata, handles)
% hObject    handle to TP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.sbplt == 1
    errordlg('Change the view','One image should be displayed');
    return
end

if isempty(handles.overChoice)
    return
end

inx = handles.size(handles.indxLBL);
if length(inx{1,1})==2
    sno = 1;
else
    sno = inx{1,1}(3);
end
r = inx{1,1}(1);
c = inx{1,1}(2);
S = handles.S;
lbl = handles.image(handles.indxLBL);

if isempty(lbl)
    f = errordlg('Load the labels','Labels not found');
    return
end
lbl = lbl{1,1};



[y,x] = getpts; y = round(y(:,1)); x = round(x(:,1));
lx = length(x);
linearInd = zeros(lx,1);
mark = [];

for i = 1:lx
    if (x(i)>0) && (x(i)<r) && (y(i)>0) && (y(i)<c)
        linearInd(i) = sub2ind([r,c,sno], x(i), y(i), S);
    end
end
linearInd = nonzeros(linearInd);
l = lbl(linearInd);

for j = 1:length(l)
    if sum(handles.fp==l(j))>0
        answer = questdlg('Do you want to rewrite FP with TP?', 'Overlapping TP and FP', 'Yes','No','No');
        switch answer
            case 'Yes'
                idex = handles.fp==l(j);
                handles.fp(idex) = 0;
                handles.fp_idx(idex,:) = [0 0 0];
            case 'No'
                l(j)=0;
                x(j)=0;
                y(j)=0;
        end
    end
end
l = nonzeros(l);
x = nonzeros(x);
y = nonzeros(y);
h = [];
for i = 1:size(handles.fp_idx,1)
    if sum(handles.fp_idx(i,:))~=0
        h = [h; handles.fp_idx(i,:)];
    end
end
handles.fp_idx = h;
        
handles.fp = nonzeros(handles.fp);


indexies = [handles.tp; l];
[C, ia, ic] = unique(indexies);
handles.tp = nonzeros(C);

slice = S*ones(length(x),1);
idx = [x,y,slice];
handles.tp_idx = [handles.tp_idx; idx];
handles.tp_idx = handles.tp_idx(ia,:);

set(handles.TP_string, 'String', sprintf('= %d', length(handles.tp)));

[P,R,F1] = Evaluation(length(handles.tp),length(handles.fp),handles.fn);
set(handles.presition, 'String', sprintf('Precision = %0.2f', P));
set(handles.recall, 'String', sprintf('Recall = %0.2f', R));
set(handles.F1, 'String', sprintf('F1 score = %0.2f', F1));


Disp(hObject, eventdata, handles);
guidata(hObject, handles);

assignin('base', 'tp_labels', handles.tp);


% --- Executes on button press in FP.
function FP_Callback(hObject, eventdata, handles)
% hObject    handle to FP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.sbplt == 1
    errordlg('Change the view','One image should be displayed');
    return
end

if isempty(handles.overChoice)
    return
end

inx = handles.size(handles.indxLBL);
if length(inx{1,1})==2
    sno = 1;
else
    sno = inx{1,1}(3);
end
r = inx{1,1}(1);
c = inx{1,1}(2);
S = handles.S;
lbl = handles.image(handles.indxLBL);

if isempty(lbl)
    f = errordlg('Load the labels','Labels not found');
    return
end
lbl = lbl{1,1};



[y,x] = getpts; y = round(y(:,1)); x = round(x(:,1));
lx = length(x);
linearInd = zeros(lx,1);

for i = 1:lx
    if (x(i)>0) && (x(i)<r) && (y(i)>0) && (y(i)<c)
        linearInd(i) = sub2ind([r,c,sno], x(i), y(i), S);
    end
end
linearInd = nonzeros(linearInd);
l = lbl(linearInd);

for j = 1:length(l)
    if sum(handles.tp==l(j))>0
        answer = questdlg('Do you want to rewrite TP with FP?', 'Overlapping TP and FP', 'Yes','No','No');
        switch answer
            case 'Yes'
                idex = handles.tp==l(j);
                handles.tp(idex) = 0;
                handles.tp_idx(idex,:) = [0 0 0];
            case 'No'
                l(j)=0;
                x(j)=0;
                y(j)=0;
        end
    end
end
l = nonzeros(l);
x = nonzeros(x);
y = nonzeros(y);
h = [];
for i = 1:size(handles.tp_idx,1)
    if sum(handles.tp_idx(i,:))~=0
        h = [h; handles.tp_idx(i,:)];
    end
end
handles.tp_idx = h;
        
handles.tp = nonzeros(handles.tp);


indexies = [handles.fp; l];
[C, ia, ic] = unique(indexies);
handles.fp = nonzeros(C);

slice = S*ones(length(x),1);
idx = [x,y,slice];
handles.fp_idx = [handles.fp_idx; idx];
handles.fp_idx = handles.fp_idx(ia,:);

set(handles.FP_string, 'String', sprintf('= %d', length(handles.fp)));

[P,R,F1] = Evaluation(length(handles.tp),length(handles.fp),handles.fn);
set(handles.presition, 'String', sprintf('Precision = %0.2f', P));
set(handles.recall, 'String', sprintf('Recall = %0.2f', R));
set(handles.F1, 'String', sprintf('F1 score = %0.2f', F1));

Disp(hObject, eventdata, handles);
guidata(hObject, handles);

assignin('base', 'fp_labels', handles.fp);


% --- Executes on button press in FN.
function FN_Callback(hObject, eventdata, handles)
% hObject    handle to FN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.sbplt == 1
    errordlg('Change the view','One image should be displayed');
    return
end

if isempty(handles.overChoice)
    return
end

[y,x] = getpts; x = round(x(:,1));
lx = length(x);

% for i = 1:lx
%     text(y(i),x(i),'FN','Color','yellow','FontSize',10)
% end

handles.fn = handles.fn + lx;

S = handles.S;
slice = S*ones(length(x),1);
idx = [x,y,slice];
handles.fn_idx = [handles.fn_idx; idx];

set(handles.FN_string, 'String', sprintf('= %d', handles.fn));

[P,R,F1] = Evaluation(length(handles.tp),length(handles.fp),handles.fn);
set(handles.presition, 'String', sprintf('Precision = %0.2f', P));
set(handles.recall, 'String', sprintf('Recall = %0.2f', R));
set(handles.F1, 'String', sprintf('F1 score = %0.2f', F1));

Disp(hObject, eventdata, handles);
guidata(hObject, handles);
assignin('base', 'fn_num', handles.fn);


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.overChoice)
    return
end

if isempty(handles.tp) && isempty(handles.fp) && handles.fn == 0
    f = errordlg('Make a categorization','Categorization');
    return
end

if isempty(handles.tp) 
    handles.tp = 0;
end
if isempty(handles.fp)
    handles.fp = 0;
end

evaluation_stats.TP = handles.tp;
evaluation_stats.FP = handles.fp;
evaluation_stats.FN = handles.fn;

uisave('evaluation_stats', [handles.dispChoice{1,1}(1:end-4), '_evaluation_stats']);  

% --- Executes on selection change in seed_location.
function seed_location_Callback(hObject, eventdata, handles)
% hObject    handle to seed_location (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns seed_location contents as cell array
%        contents{get(hObject,'Value')} returns selected item from seed_location
contents = cellstr(get(hObject,'String'));
choice = contents(get(hObject, 'Value'));
if strcmp(choice, 'Regional maxima')
    handles.seed_loc = 'Regional maxima';
elseif strcmp(choice, 'Watershed centroids')
    handles.seed_loc = 'Watershed centroids';
elseif strcmp(choice, 'Manual selection')
    handles.seed_loc = 'Manual selection';
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function seed_location_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seed_location (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function grid_size_Callback(hObject, eventdata, handles)
% hObject    handle to grid_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of grid_size as text
%        str2double(get(hObject,'String')) returns contents of grid_size as a double
grid = get(handles.grid_size,'String');
handles.morph.grid = str2double(grid);
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function grid_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to grid_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in morphometry_button.
function morphometry_button_Callback(hObject, eventdata, handles)
% hObject    handle to morphometry_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.morph.im) && isempty(handles.morph.myelin)
    errordlg('Select the label or IAS or myelin','Label not found');
    return
end

idx = get(handles.morph_IAS_lbl, 'Value')-1;
vox_size = handles.vox_size(idx);

res = vox_size{1,1};
grid_size = handles.morph.grid;
myelin = handles.morph.myelin;
myelin = myelin{1,1};
label = handles.morph.im;
label = label{1,1};
name = handles.morph.name{1,1};
name = ['morphometry_' name(5:end)];

% intraaxon same label as myelin
tic
label(label>0)=1;
label = bwlabeln(label);
h = myelin;

myelin = zeros(size(label,1),size(label,2),size(label,3));

f = waitbar(0,'Mophometry in process...');
NumLbl = max(label(:));

for i = 1:NumLbl
    %1:max(label(:))
    a = false(size(label,1),size(label,2),size(label,3));
    a(label==i) = 1;
    if sum(a(:))<3000
        label(label==i)=0;
    end
    if sum(a(:))>=3000
        a = imdilate(a,true(10,10,10));
        k = immultiply(a,h);
        s = regionprops(k,'Area');
        Areas = cat(1, s.Area);
        [val, idx] = max(Areas);
%         if idx == 42
%             4
%         end
        myelin(h==idx)=i;
%         k(k>0)=true;
%         k = logical(k);
        
    end
    frac = i/(NumLbl/100)/100;
    waitbar(frac,f,sprintf('Myelin+IAS - same label: %d/%d axons done', i, NumLbl));
end
%




for i = 1:NumLbl
    a = false(size(label,1),size(label,2),size(label,3));
    a(label==i) = 1;
    IAS = [];
    Myelin = [];
    if sum(a(:))>500
        b = false(size(myelin,1),size(myelin,2),size(myelin,3));
        b(myelin==i) = 1;
        final_axon = a;
        if sum(b(:))>500
           [IAS, Myelin] = util_axonal_quantification(final_axon,res,b);
        end
    end
%            isosurface(a)  ; hold on ; p = patch(isosurface(b,0.5)); 
%          p.FaceColor = 'red'; p.EdgeColor = 'none'; p.FaceAlpha = 0.3;
%           set(gca,'PlotBoxAspectRatio',[1 1 3.6])
    frac = i/(NumLbl/100)/100;
    waitbar(frac,f,sprintf('Morphometry %d/%d axons done', i, NumLbl));

    axon(i).IAS = IAS;
    axon(i).Myelin = Myelin;
end

close(f)
toc

axon(1).LabelIAS = label;
axon(1).LabelMyelin = myelin;



[filename, filepath] = uiputfile('*.mat', 'Save the project file:',name);
FileName = fullfile(filepath, filename);
save(FileName, 'axon', '-v7.3');


handles.morph.axon = axon;
guidata(hObject, handles);


% --- Executes on selection change in morph_IAS_lbl.
function morph_IAS_lbl_Callback(hObject, eventdata, handles)
% hObject    handle to morph_IAS_lbl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns morph_IAS_lbl contents as cell array
%        contents{get(hObject,'Value')} returns selected item from morph_IAS_lbl
contents = cellstr(get(hObject,'String'));
morphIM_choice = contents(get(hObject, 'Value'));

if strcmp(morphIM_choice, 'Select IAS label...')
    return
else
    indx = strcmp(handles.listFiles.String, morphIM_choice);
    indx = indx(2:end);
    handles.morph.im = handles.image(indx);
    handles.morph.name = morphIM_choice;
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function morph_IAS_lbl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to morph_IAS_lbl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Morph_myelin_lbl.
function Morph_myelin_lbl_Callback(hObject, eventdata, handles)
% hObject    handle to Morph_myelin_lbl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Morph_myelin_lbl contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Morph_myelin_lbl
contents = cellstr(get(hObject,'String'));
morphMyel_choice = contents(get(hObject, 'Value'));

if strcmp(morphMyel_choice, 'Select myelin label...')
    return
else
    indx = strcmp(handles.listFiles.String, morphMyel_choice);
    indx = indx(2:end);
    handles.morph.myelin = handles.image(indx);
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Morph_myelin_lbl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Morph_myelin_lbl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function ruler_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ruler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.dispChoice)
    return
end

if isempty(handles.res)
    prompt = {'Resolution of file in 2D in \mum (Example: 0.0137498)'};
    dlgtitle = 'Resolution';
    dims = [1 50];
    definput = {'0.0137498'};
    opts.Interpreter = 'tex';
    answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
    if isempty(answer)
        return
    end

    res = answer{1,1};
    res = regexprep(res,',','.');
    res = str2double(res);

    while isnan(res)
        dlgtitle = 'Not a number. Input number:';
        dims = [1 70];
        answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
        if isempty(answer)
            %res = 0;
            return
        end

        res = answer{1,1};
        res = regexprep(res,',','.');
        res = str2double(res);
    end
    handles.res = res;
end

res = handles.res;

h = imline;
pos = getPosition(h);
D = round(pdist(pos)*res,2);
f = msgbox(['Distance = ' num2str(D) ' um']);

delete(h)

guidata(hObject, handles);


% --------------------------------------------------------------------
function subplot_OffCallback(hObject, eventdata, handles)
% hObject    handle to subplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.dispChoice)
    return
end

set(handles.SlidesText,'Visible','on')
set(handles.sliderSlices,'Visible','on')
set (gcf, 'WindowScrollWheelFcn', {@mouseScroll, handles}); 

handles.sbplt = 0;

set(handles.Picture, 'Position', [0.165 0.315 0.83 0.705]);
set(handles.hs2,'Visible','off')
set(handles.hs3,'Visible','off')
set(handles.hs4,'Visible','off')
set(handles.subplot4,'Visible','off')
set(handles.slider_sbplt1,'Visible','off')
set(handles.slider_sbplt2,'Visible','off')
set(handles.slider_sbplt3,'Visible','off')
set(get(gca,'children'),'Visible','off')
Disp(hObject, eventdata, handles);

guidata(hObject, handles);

% --------------------------------------------------------------------
function subplot_OnCallback(hObject, eventdata, handles)
% hObject    handle to subplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.dispChoice)
    return
end

handles.sbplt = 1;
set(handles.SlidesText,'Visible','off')
set(handles.sliderSlices,'Visible','off')

indx = handles.indx;
image = handles.image(indx);
data = uint8(image{1,1});
[r,c,s] = size(data);

%subplo1
set(handles.Picture, 'Position', [0.185 0.67 0.38 0.3]);
%slider
set(handles.slider_sbplt1,'Visible','on')
set(handles.slider_sbplt1,'Units', 'normalized');
set(handles.slider_sbplt1, 'Position', [0.165 0.67 0.015 0.3]);

%subplot2
set(handles.subplot2,'Visible','on')
set(handles.subplot2,'Units', 'normalized');
set(handles.subplot2, 'Position', [0.57 0.67 0.38 0.3]);
%slider
set(handles.slider_sbplt2,'Visible','on')
set(handles.slider_sbplt2,'Units', 'normalized');
set(handles.slider_sbplt2, 'Position', [0.97 0.67 0.015 0.3]);

%subplot3
set(handles.subplot3,'Visible','on')
set(handles.subplot3,'Units', 'normalized');
set(handles.subplot3, 'Position', [0.185 0.35 0.38 0.3]);
%slider
set(handles.slider_sbplt3,'Visible','on')
set(handles.slider_sbplt3,'Units', 'normalized');
set(handles.slider_sbplt3, 'Position', [0.165 0.34 0.015 0.3]);

%subplot4
set(handles.subplot4,'Visible','on')
set(handles.subplot4,'Units', 'normalized');
set(handles.subplot4, 'Position', [0.57 0.35 0.4 0.3]);

guidata(hObject, handles);
subplots(hObject, eventdata, handles, handles.subplots.slider1, handles.subplots.slider2, handles.subplots.slider3)

function subplots(hObject, eventdata, handles, s1, s2, s3)
indx = handles.indx;
image = handles.image(indx);
data = uint8(image{1,1});
[r,c,s] = size(data);

%subplot1
axes(handles.Picture)
pic = squeeze(data(:,:,s1));
hs1 = imshow(pic);
handles.hs1 = hs1;

%subplot2
axes(handles.subplot2)
pic = squeeze(data(:,s2,:))';
hs2 = imshow(pic);
handles.hs2 = hs2;

%subplot3
axes(handles.subplot3)
pic = squeeze(data(s3,:,:))';
hs3 = imshow(pic);
handles.hs3 = hs3;

%subplot4
axes(handles.subplot4)
alpha=0.3;

X = [0 0 0 0 0 1; 1 0 1 1 1 1; 1 0 1 1 1 1; 0 0 0 0 0 1];
Y = [0 0 0 0 1 0; 0 1 0 0 1 1; 0 1 1 1 1 1; 0 0 1 1 1 0];
Z = [0 0 1 0 0 0; 0 0 1 0 0 0; 1 1 1 0 1 1; 1 1 1 0 1 1];

C='black';   

X = r*(X) ;
Y = c*(Y) ;
Z = s*(Z) ; 

hs4 = plot3(X,Y,Z,C);
axis equal

extra = 50;
step1 = round(s/4);

hold on
[x y] = meshgrid(-extra:step1:r+extra,-extra:step1:c+extra);
z = s1*ones(size(x));
C = x.*y;
C(:,:,1) = 0.95*ones(size(x)); 
C(:,:,2) = zeros(size(x));
C(:,:,3) = zeros(size(x));
surf(x,y,z,C,'FaceAlpha',alpha,'EdgeColor', 'none');

step2 = round(c/4);
[x y] = meshgrid(-extra:step2:r+extra,-extra:step2:s+extra);
z = s2*ones(size(x));
C = x.*y;
C(:,:,1) = zeros(size(x)); 
C(:,:,2) = zeros(size(x));
C(:,:,3) = zeros(size(x));
surf(x,z,y,C,'FaceAlpha',alpha,'EdgeColor', 'none');

step3 = round(r/4);
[x y] = meshgrid(-extra:step3:s+extra,-extra:step3:c+extra);
z = s3*ones(size(x));
C = x.*y;
C(:,:,1) = zeros(size(x)); 
C(:,:,2) = zeros(size(x));
C(:,:,3) = 0.95*ones(size(x));
surf(z,y,x,C,'FaceAlpha',alpha,'EdgeColor', 'none');


hold off

handles.hs4 = hs4;

guidata(hObject, handles);


function mouseScrollsubplot1(hObject, eventdata,handles)
handles = guidata(hObject);

S  = get(handles.slider_sbplt1,'Value');

indx = handles.indx;
image = handles.image(indx);
data = uint8(image{1,1});
[r,c,s] = size(data);

UPDN = eventdata.VerticalScrollCount;
S = S - UPDN;
if (S < 1); S = 1; elseif (S > s); S = s; end

set(handles.slider_sbplt1,'Value',S);

handles.subplots.slider1 = ceil(S);

guidata(hObject, handles);
subplots(hObject, eventdata, handles, handles.subplots.slider1, handles.subplots.slider2, handles.subplots.slider3)




% --- Executes on slider movement.
function slider_sbplt1_Callback(hObject, eventdata, handles)
% hObject    handle to slider_sbplt1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = handles.indx;
image = handles.image(indx);
data = uint8(image{1,1});
[r,c,s] = size(data);

set(handles.slider_sbplt1, 'Max',s);
set(handles.slider_sbplt1, 'Min', 1);
set(handles.slider_sbplt1, 'SliderStep' , [0.01,0.01] );
S = get(handles.slider_sbplt1,'Value');
if (S < 1); S = 1; elseif (S > s); S = s; end

handles.subplots.slider1 = ceil(S);
set(handles.slider_sbplt1,'Value', S)

set (gcf, 'WindowScrollWheelFcn', {@mouseScrollsubplot1, handles});

guidata(hObject, handles);
subplots(hObject, eventdata, handles, handles.subplots.slider1, handles.subplots.slider2, handles.subplots.slider3)



% --- Executes during object creation, after setting all properties.
function slider_sbplt1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_sbplt1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function mouseScrollsubplot2(hObject, eventdata,handles)
handles = guidata(hObject);

S  = get(handles.slider_sbplt2,'Value');

indx = handles.indx;
image = handles.image(indx);
data = uint8(image{1,1});
[r,c,s] = size(data);

UPDN = eventdata.VerticalScrollCount;
S = S - UPDN;
if (S < 1); S = 1; elseif (S > c); S = c; end

set(handles.slider_sbplt2,'Value',S);

handles.subplots.slider2 = ceil(S);

guidata(hObject, handles);
subplots(hObject, eventdata, handles, handles.subplots.slider1, handles.subplots.slider2, handles.subplots.slider3)

% --- Executes on slider movement.
function slider_sbplt2_Callback(hObject, eventdata, handles)
% hObject    handle to slider_sbplt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = handles.indx;
image = handles.image(indx);
data = uint8(image{1,1});
[r,c,s] = size(data);

set(handles.slider_sbplt2, 'Max',c);
set(handles.slider_sbplt2, 'Min', 1);
set(handles.slider_sbplt2, 'SliderStep' , [0.01,0.01] );
S = get(handles.slider_sbplt2,'Value');
if (S < 1); S = 1; elseif (S > c); S = c; end

handles.subplots.slider2 = ceil(S);
set(handles.slider_sbplt2,'Value', S)

set (gcf, 'WindowScrollWheelFcn', {@mouseScrollsubplot2, handles});

guidata(hObject, handles);
subplots(hObject, eventdata, handles, handles.subplots.slider1, handles.subplots.slider2, handles.subplots.slider3)


% --- Executes during object creation, after setting all properties.
function slider_sbplt2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_sbplt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function mouseScrollsubplot3(hObject, eventdata,handles)
handles = guidata(hObject);

S  = get(handles.slider_sbplt3,'Value');

indx = handles.indx;
image = handles.image(indx);
data = uint8(image{1,1});
[r,c,s] = size(data);

UPDN = eventdata.VerticalScrollCount;
S = S - UPDN;
if (S < 1); S = 1; elseif (S > r); S = r; end

set(handles.slider_sbplt3,'Value',S);

handles.subplots.slider3 = ceil(S);

guidata(hObject, handles);
subplots(hObject, eventdata, handles, handles.subplots.slider1, handles.subplots.slider2, handles.subplots.slider3)

% --- Executes on slider movement.
function slider_sbplt3_Callback(hObject, eventdata, handles)
% hObject    handle to slider_sbplt3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
indx = handles.indx;
image = handles.image(indx);
data = uint8(image{1,1});
[r,c,s] = size(data);

set(handles.slider_sbplt3, 'Max',r);
set(handles.slider_sbplt3, 'Min', 1);
set(handles.slider_sbplt3, 'SliderStep' , [0.01,0.01] );
S = get(handles.slider_sbplt3,'Value');
if (S < 1); S = 1; elseif (S > r); S = r; end

handles.subplots.slider3 = ceil(S);
set(handles.slider_sbplt3,'Value', S)

set (gcf, 'WindowScrollWheelFcn', {@mouseScrollsubplot3, handles});

guidata(hObject, handles);
subplots(hObject, eventdata, handles, handles.subplots.slider1, handles.subplots.slider2, handles.subplots.slider3)


% --- Executes during object creation, after setting all properties.
function slider_sbplt3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_sbplt3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function pos = draw(hObject, eventdata, handles,image)
    h = drawfreehand('LineWidth',5,'Closed',0,'HandleVisibility','off', 'FaceAlpha', 0,'Deletable',false, 'InteractionsAllowed', 'none');
    h.Waypoints(h.Waypoints==1)=0;
    mask = createMask(h,image);
    mask = imdilate(mask,true(4));
    pos2 = find(mask==1);
    [pos(:,2), pos(:,1)] = ind2sub(size(image),pos2);
    pos(:,3) = handles.S*ones(length(pos2),1);
    delete(h)
    
    %linPos = sub2ind(size(image3D), pos(:,1), pos(:,2), pos(:,3));
    
    
    %intensities = image(pos2);

    guidata(hObject, handles);


% --- Executes on button press in myelin_ML.
function myelin_ML_Callback(hObject, eventdata, handles)
% hObject    handle to myelin_ML (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
indx = strcmp(handles.listFiles.String, handles.opt.read_im);
indx = indx(2:end);
handles.indx = indx;
image = handles.image(indx);
image = image{1,1}(:,:,handles.S);

pos = draw(hObject, eventdata,handles,image);
f = features(hObject, eventdata,handles,image,pos);
fs = [f ones(size(f,1),1)];
handles.ML.train = [handles.ML.train; fs];

handles.ML.myelin.pos = [handles.ML.myelin.pos; pos];

Disp(hObject, eventdata, handles);

Im_rslt = Random_forest(hObject, eventdata,handles,image,handles.S);

fig = figure(2);
imshow(uint8(image))
%set(fig, 'Position', [1000 600 700 700]);
hold on
% Lrgb = label2rgb(Im_rslt,'jet', 'k', 'shuffle');
Lrgb = label2rgb(Im_rslt,'lines','k');
himage = imshow(Lrgb); himage.AlphaData = 0.3;


handles = guidata(hObject);
guidata(hObject, handles);


% --- Executes on button press in background_ML.
function background_ML_Callback(hObject, eventdata, handles)
% hObject    handle to background_ML (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
indx = strcmp(handles.listFiles.String, handles.opt.read_im);
indx = indx(2:end);
handles.indx = indx;
image = handles.image(indx);
image = image{1,1}(:,:,handles.S);

pos = draw(hObject, eventdata,handles,image);
f = features(hObject, eventdata,handles,image,pos);
fs = [f 2*ones(size(f,1),1)];
handles.ML.train = [handles.ML.train; fs];

handles.ML.bg.pos = [handles.ML.bg.pos; pos];

Disp(hObject, eventdata, handles);

Im_rslt = Random_forest(hObject, eventdata,handles,image,handles.S);

figure(2)
imshow(uint8(image))
hold on
% Lrgb = label2rgb(Im_rslt,'jet', 'k', 'shuffle');
Lrgb = label2rgb(Im_rslt,'lines', 'k');
himage = imshow(Lrgb); himage.AlphaData = 0.3;

handles = guidata(hObject);
guidata(hObject, handles);


% --- Executes on button press in ML_done.
function ML_done_Callback(hObject, eventdata, handles)
% hObject    handle to ML_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.mach_learn,'Visible','off')
set(handles.myelin_ML,'Visible','off')
set(handles.background_ML,'Visible','off')
set(handles.ML_done,'Visible','off')

indx = strcmp(handles.listFiles.String, handles.opt.read_im);
indx = indx(2:end);
handles.indx = indx;
image = handles.image(indx);
image = image{1,1};

f = waitbar(0,'Myelin segmentation of the whole volume in progress...');

for i = 1:size(image,3)
    waitbar(i/size(image,3),f,'Segmentation of the whole volume in progress...');
    lbl(:,:,i) = Random_forest(hObject, eventdata,handles,image(:,:,i),i);
end

close(f)

handles.ML.label = lbl;

handles.ML.myelin.pos = [];
handles.ML.bg.pos = [];
handles.ML.train = [];
guidata(hObject, handles);
Disp(hObject, eventdata, handles);




function [result] = Random_forest(hObject, eventdata,handles,im,slice)

%code
if isempty(handles.ML.myelin.pos) || isempty(handles.ML.bg.pos)
    result = ones(size(im));
    return
end
    
%train
train = handles.ML.train;

X = single(train(:,1:end-1));
Y = train(:,end);

TreeObject=TreeBagger(100,X,Y,'Method','classification','NVarToSample','all');

%test
pos1 = handles.ML.myelin.pos;
idxPos = pos1(:,3)==slice;
pos1 = pos1(idxPos,1:2);
linPos1 = sub2ind(size(im), pos1(:,2), pos1(:,1));
Y1 = Y(idxPos);

pos2 = handles.ML.bg.pos;
idxPos = pos2(:,3)==slice;
pos2 = pos2(idxPos,1:2);
linPos2 = sub2ind(size(im), pos2(:,2), pos2(:,1));
Y2 = Y(idxPos);

mask = ones(size(im));
mask(linPos1) = 0;
mask(linPos2) = 0;
test_pos = find(mask);
[tPos(:,2),tPos(:,1)] = ind2sub(size(im),test_pos);

test = features(hObject, eventdata,handles,im,tPos);
test = single(test);

tic
L = size(test,1);
step = round(L/50);
seq = 1:step:L;
YFIT = cell(1,10);

parfor i = 1:length(seq)-1
    YFIT{i} = predict(TreeObject,test(seq(i):(seq(i+1)-1),:));
end

YFIT{length(seq)} = predict(TreeObject,test(seq(end):L,:));
toc

pred = [];
for i = 1:length(YFIT)
    pred = [pred;YFIT{1,i}];
end

result = zeros(size(im));
result(test_pos) = str2num(cell2mat(pred));
result(linPos1) = 1;
result(linPos2) = 2;
%result = ~logical(result);


% --- Executes on selection change in listFiles.
function listFiles_Callback(hObject, eventdata, handles)
% hObject    handle to listFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listFiles contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listFiles

index_selected = get(handles.listFiles, 'Value');
file_list = get(handles.listFiles, 'String');
if ~iscell(file_list)
    handles.sel_filename = file_list;    
else    
    filename = file_list{index_selected};
    handles.sel_filename = filename;
    set(handles.save_from_list_file, 'Enable', 'on')
end
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function Picture_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Picture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate Picture

% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Picture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate Picture


% --- Executes on key press with focus on listFiles and none of its controls.
function listFiles_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to listFiles (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
keyPressed = eventdata.Key;
if strcmpi(keyPressed,'delete')
    file = handles.sel_filename;
    
    if strcmp(file, 'File list:')
        return
    end
    
    indx = strcmp(handles.listFiles.String, file);
    
    answer = questdlg(['Do you want to remove ', file, '?'] ,'Delete','Yes', 'No','Yes');
    % Handle response
    switch answer
        case 'No'
            return
        case 'Yes'
            %data delete
            handles.image(indx(2:end)) = [];
            handles.size(indx(2:end)) = [];
            handles.vox_size(indx(2:end)) = [];
            
            %file list delete
            txt = get(handles.listFiles,'String');
            txt(indx) = [];
            set(handles.listFiles,'String',txt);
            set(handles.listFiles,'Value',size(txt,1));
            
            
            % popup menus delete
            poptxt = get(handles.DisplayPopUp,'String');
            poptxt(indx) = [];
            set(handles.DisplayPopUp,'String',poptxt);
            
            poptxt2 = get(handles.OverlayPopUp,'String');
            poptxt2(indx) = [];
            set(handles.OverlayPopUp,'String',poptxt2);

            %read file pop-up menu
            rdfile = get(handles.readImage,'String');
            rdfile(indx) = [];
            set(handles.readImage,'String',rdfile);

            %morphometry
            morphfile1 = get(handles.morph_IAS_lbl,'String');
            morphfile1(indx) = [];
            set(handles.morph_IAS_lbl,'String',morphfile1);

            morphfile2 = get(handles.Morph_myelin_lbl,'String');
            morphfile2(indx) = [];
            set(handles.Morph_myelin_lbl,'String',morphfile2);
            
            if strcmp(file,handles.dispChoice)
                handles.dispChoice = 'No Display';
                set(handles.DisplayPopUp, 'Value', 1)
                
            elseif strcmp(file,handles.overChoice)
                handles.overChoice = 'No Overlay';
                set(handles.OverlayPopUp, 'Value', 1)
            end
            guidata(hObject, handles);
            Disp(hObject, eventdata, handles);
                
    end
  
end




% --------------------------------------------------------------------
function save_from_list_file_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to save_from_list_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename = handles.sel_filename;
if strcmp(filename, 'File list:')
    return
end

indx = strcmp(handles.listFiles.String, filename);
file = handles.image(indx(2:end));
file = file{1,1};
[file_name, filepath] = uiputfile('*.mat', 'Save file:',filename);
if logical(file_name)
    FileName = fullfile(filepath, file_name);
    save(FileName, 'file', '-v7.3');
end



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Erase_button.
function Erase_button_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Erase_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on Erase_button and none of its controls.
function Erase_button_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Erase_button (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
