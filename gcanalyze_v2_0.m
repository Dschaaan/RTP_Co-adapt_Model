%--------------------------------------------------------------------------
% MATLAB GUI source code of a tool developed to count axon fibers/bundles 
% that grow (in parallel) from retinal explant strips for the 
% quantification of in vitro gap-assays. It is associated with the paper 
% Fiedeling, F., et al., „Ephrin-A/EphA specific co-adaptation as a novel 
% mechanism in topographic axon guidance“.

% Copyright (C) 2017, Franco Weth (franco.weth@kit.edu)

% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU general Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version. This program is distributed in the hope that 
% it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
% GNU General Public License for more details.
%--------------------------------------------------------------------------

%%
function varargout = gcanalyze_v2_0(varargin)
% GCANALYZE_V2_0 MATLAB code for gcanalyze_v2_0.fig
%      GCANALYZE_V2_0, by itself, creates a new GCANALYZE_V2_0 or raises the existing
%      singleton*.
%
%      H = GCANALYZE_V2_0 returns the handle to a new GCANALYZE_V2_0 or the handle to
%      the existing singleton*.
%
%      GCANALYZE_V2_0('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GCANALYZE_V2_0.M with the given input arguments.
%
%      GCANALYZE_V2_0('Property','Value',...) creates a new GCANALYZE_V2_0 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gcanalyze_v2_0_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gcanalyze_v2_0_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gcanalyze_v2_0

% Last Modified by GUIDE v2.5 03-Feb-2017 13:51:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gcanalyze_v2_0_OpeningFcn, ...
                   'gui_OutputFcn',  @gcanalyze_v2_0_OutputFcn, ...
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


% --- Executes just before gcanalyze_v2_0 is made visible.
function gcanalyze_v2_0_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gcanalyze_v2_0 (see VARARGIN)

% Choose default command line output for gcanalyze_v2_0
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gcanalyze_v2_0 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gcanalyze_v2_0_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%______________BROWSE 1____________________________________________________

% --- Executes on button press in browse_1.
function browse_1_Callback(hObject, ~, handles) %#ok<*DEFNU>
% hObject    handle to browse_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output=hObject;
[filename,path]=uigetfile('*.tif');    
filegap=fullfile(path,filename);
savepath ('padthdef.m')

%endung='*.tif';
%path=[pathdef endung];
imshow(filegap,'Parent',handles.preview_whole);
set(handles.checkbox1,'Value', 1); 

global unterkante_positionen;
unterkante_positionen=ginput(2); 
unterkante=polyfit(unterkante_positionen(:,1),unterkante_positionen(:,2),1);

x=0:0.1:100;
y=polyval(unterkante,x);
hold on
plot(x,y,'k')
plot(unterkante_positionen(:,1),unterkante_positionen(:,2),'r') 
set(handles.checkbox2,'Value', 1); 

global oberkante_positionen 
oberkante_positionen=ginput(2);
oberkante=polyfit(oberkante_positionen(:,1),oberkante_positionen(:,2),1);

w=0:0.1:100;
z=polyval(oberkante,w);
hold on
plot(w,z,'k')
plot(oberkante_positionen(:,1),oberkante_positionen(:,2),'r') 
set(handles.checkbox3,'Value', 1); 

global dist
global dist_txt
global pixscale

if get(handles.button_5Ob,'Value')==1
    pixscale=1.29;
elseif get(handles.button_10Ob,'Value')==1
    pixscale=0.65;
else
    pixscale=0;
end

dist1=abs((unterkante_positionen(3))-(oberkante_positionen(3)));
dist1_txt=num2str(dist);
dist2=abs((unterkante_positionen(4))-(oberkante_positionen(4)));
dist2_txt=num2str(dist);
dist=((dist1+dist2)/2)*pixscale;
dist_txt=num2str(dist);

set(handles.gap_width,'String',dist_txt);



function gap_width_Callback(hObject, eventdata, handles)
% hObject    handle to gap_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gap_width as text
%        str2double(get(hObject,'String')) returns contents of gap_width as a double
global dist_txt
set(hObject,'String',dist_txt);


% --- Executes during object creation, after setting all properties.
function gap_width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gap_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in button_5Ob.
function button_5Ob_Callback(hObject, eventdata, handles)
% hObject    handle to button_5Ob (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global pixscale

get(hObject,'Value') 
if get(handles.button_5Ob,'Value')==1
    set(handles.button_10Ob,'Value',0);
    pixscale=1.29;
    set(handles.checkobjective,'Value',1);
end


% --- Executes on button press in button_10Ob.
function button_10Ob_Callback(hObject, eventdata, handles)
% hObject    handle to button_10Ob (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global pixscale

get(hObject,'Value') 
if get(handles.button_10Ob,'Value')==1
    set(handles.button_5Ob,'Value',0);
    pixscale=0.65;
    set(handles.checkobjective,'Value',1);
end


%____________BROWSE 2______________________________________________________

% --- Executes during object creation, after setting all properties.
function browse_2_Callback(~, ~, handles)
% hObject    handle to browse_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% load pathdef.m;
[filename,path]=uigetfile('*.tif');

file=fullfile(path,filename);
imshow(file,'Parent',handles.preview_whole);

global unterkante_positionen;
global oberkante_positionen;

hold on
plot(unterkante_positionen(:,1),unterkante_positionen(:,2),'r') 
plot(oberkante_positionen(:,1),oberkante_positionen(:,2),'r') 
set(handles.checkbox4,'Value', 1);

global ingap;
ingap=imcrop(gcf);
set(handles.checkbox5,'Value', 1); 

global aftergap;
aftergap=imcrop(gcf);
set(handles.checkbox6,'Value', 1);
ingapbw=im2bw(ingap,0.2);
aftergapbw=im2bw(aftergap,0.2);

imshow(ingapbw,'Parent',handles.preview_select1);
imshow(aftergapbw,'Parent',handles.preview_select2);

threshvalue=0.2;
set(handles.thresh1,'Value',0.2);
set(handles.thresh1_text,'String',num2str(threshvalue));
set(handles.thresh2,'Value',0.2);
set(handles.thresh2_text,'String',num2str(threshvalue));

%___________THRESH 1_______________________________________________________

% --- Executes on slider movement.
function thresh1_Callback(~, ~, handles)
% hObject    handle to thresh1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global ingap;
global ingapbw;
threshvalue = get(handles.thresh1,'Value');
set(handles.thresh1_text,'String',num2str(threshvalue));

ingapbw=im2bw(ingap,threshvalue);
imshow(ingapbw,'Parent',handles.preview_select1);

% --- Executes during object creation, after setting all properties.
function thresh1_CreateFcn(hObject, ~, ~)
% hObject    handle to thresh1_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function thresh1_text_Callback(~, ~, ~)
% hObject    handle to thresh1_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresh1_text as text
%        str2double(get(hObject,'String')) returns contents of thresh1_text as a double


% --- Executes during object creation, after setting all properties.
function thresh1_text_CreateFcn(hObject, ~, ~)
% hObject    handle to thresh1_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%___________THRESH 2_______________________________________________________
% --- Executes on slider movement.
function thresh2_Callback(~, ~, handles)
% hObject    handle to thresh2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global aftergap;
global aftergapbw;
threshvalue = get(handles.thresh2,'Value');
set(handles.thresh2_text,'String',num2str(threshvalue));
aftergapbw=im2bw(aftergap,threshvalue);
imshow(aftergapbw,'Parent',handles.preview_select2);

% --- Executes during object creation, after setting all properties.
function thresh2_CreateFcn(hObject, ~, ~)
% hObject    handle to thresh1_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function thresh2_text_Callback(~, ~, ~)
% hObject    handle to thresh2_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresh2_text as text
%        str2double(get(hObject,'String')) returns contents of thresh2_text as a double


% --- Executes during object creation, after setting all properties.
function thresh2_text_CreateFcn(hObject, ~, ~)
% hObject    handle to thresh2_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function preview_whole_CreateFcn(hObject, ~, ~)
% hObject    handle to thresh1_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function figure1_CreateFcn(hObject, ~, ~)
% hObject    handle to thresh1_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function uipanel1_CreateFcn(hObject, ~, ~)
% hObject    handle to thresh1_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%_______________START______________________________________________________
% --- Executes on button press in start.
function start_Callback(~, ~, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global ingapbw;
global aftergapbw;

if get(handles.button_10Ob,'Value')==1
    gcsize=20;
elseif get(handles.button_5Ob,'Value')==1
    gcsize=12;
end

axsize=2;



sizeingapbw=size(ingapbw);
colingapbw=sizeingapbw(:,1);
lineingapbw=length(ingapbw);

sizeaftergapbw=size(aftergapbw);
colaftergapbw=sizeaftergapbw(:,1);
lineaftergapbw=length(aftergapbw);

ingapfiber=zeros(colingapbw,1);
aftergapfiber=zeros(colaftergapbw,1);



for ii=1:colingapbw
    signlin(ii)=0;
    
    for jj=2:lineingapbw
        if ingapbw(ii,jj)==1 && ingapbw(ii,jj-1)==0
            if signlin(ii)==0;
                sigcountin=1;
                signlin(ii)=1;
                
            else
                sigcountin=sigcountin+1;
                signlin(ii)=signlin(ii)+1;
                
            end
            stimein(signlin)=jj;
    
            intiin=1;
            
        elseif ingapbw(ii,jj)==0 && ingapbw(ii,jj-1)==1
            entiin=1;
            if exist('intiin','var')==0
                intiin=0;
            end
            if intiin==1 && entiin==1
                
                signallengthin(ii,sigcountin)=jj-stimein(sigcountin);
                
                
                entiin=0;
                intiin=0;
                
                if signallengthin(ii,sigcountin) >= axsize 
                    if signallengthin(ii,sigcountin) <= gcsize
                        ingapfiber(ii)=ingapfiber(ii)+(ceil((signallengthin(ii,sigcountin))/axsize));
                    else
                        ingapfiber(ii)=ingapfiber(ii)+1;
                    end
                else
                    ingapfiber(ii)=ingapfiber(ii)+1;
                            
                end
            end
        end
            
    end
end


for ii=1:colaftergapbw
    signlaf(ii)=0;
    
    for jj=2:lineaftergapbw
        if aftergapbw(ii,jj)==1 && aftergapbw(ii,jj-1)==0
            if signlaf(ii)==0;
                sigcountaf=1;
                signlaf(ii)=1;
                
            else
                sigcountaf=sigcountaf+1;
                signlaf(ii)=signlaf(ii)+1;
                
            end

            stimeaf(signlaf)=jj;
    
            intiaf=1;
            
        elseif aftergapbw(ii,jj)==0 && aftergapbw(ii,jj-1)==1
            entiaf=1;
            if intiaf==1 && entiaf==1
                
                signallengthaf(ii,sigcountaf)=jj-stimeaf(sigcountaf);
                
                
                entiaf=0;
                intiaf=0;
                
                if signallengthaf(ii,sigcountaf) >= axsize 
                    if signallengthaf(ii,sigcountaf) <= gcsize
                        aftergapfiber(ii)=aftergapfiber(ii)+(ceil((signallengthaf(ii,sigcountaf))/axsize));
                    else
                        aftergapfiber(ii)=aftergapfiber(ii)+1;
                    end
                else
                    aftergapfiber(ii)=aftergapfiber(ii)+1;
                            
                end
            end
        end
            
    end
end


fibers_in_gap=mean(ingapfiber);
fibers_after_gap=mean(aftergapfiber);
fibers_in_gap_med=median(ingapfiber);
fibers_after_gap_med=median(aftergapfiber);

fibers_crossing_gap_percent=fibers_after_gap/(fibers_in_gap/100);
fibers_crossing_gap_percent_med=fibers_after_gap_med/(fibers_in_gap_med/100);



res1=num2str(fibers_crossing_gap_percent);
set(handles.result,'String',num2str(res1));

res2=num2str(fibers_in_gap);
set(handles.ingap,'String',num2str(res2));

res3=num2str(fibers_after_gap);
set(handles.aftergap,'String',num2str(res3));

res4=num2str(fibers_crossing_gap_percent_med);
set(handles.median1,'String',num2str(res4));

res5=num2str(fibers_in_gap_med);
set(handles.median2,'String',num2str(res5));

res6=num2str(fibers_after_gap_med);
set(handles.median3,'String',num2str(res6));

bar(ingapfiber,'Parent',handles.result1);
xlabel(handles.result1,'# analyzed row');
ylabel(handles.result1,'fibers counted');

hist(ingapfiber,'Parent',handles.result2);
xlabel(handles.result2,'fibers counted');
ylabel(handles.result2,'relative frequency');

bar(aftergapfiber,'Parent',handles.result3);
xlabel(handles.result3,'# analyzed row');
ylabel(handles.result3,'fibers counted');

hist(aftergapfiber,'Parent',handles.result4);
xlabel(handles.result4,'fibers counted');
ylabel(handles.result4,'relative frequency');

% save('gap_analyzer');

%___________% OF AXONS OVERGROWING_________________________________________

function result_Callback(~, ~, ~)
% hObject    handle to result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of result as text
%        str2double(get(hObject,'String')) returns contents of result as a double


% --- Executes during object creation, after setting all properties.
function result_CreateFcn(hObject, ~, ~)
% hObject    handle to result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%_____________CHECKBOXES___________________________________________________

% --- Executes on button press in checkbox1.
function checkbox1_Callback(~, ~, ~)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1

% --- Executes on button press in checkbox2.
function checkbox2_Callback(~, ~, ~)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2

% --- Executes on button press in checkbox3.
function checkbox3_Callback(~, ~, ~)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3

% --- Executes on button press in checkbox4.
function checkbox4_Callback(~, ~, ~)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4

% --- Executes on button press in checkbox5.
function checkbox5_Callback(~, ~, ~)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5

% --- Executes on button press in checkbox6.
function checkbox6_Callback(~, ~, ~)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6

% --- Executes on button press in checkobjective.
function checkobjective_Callback(hObject, eventdata, handles)
% hObject    handle to checkobjective (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkobjective


%__________AXONS IN GAP____________________________________________________
function ingap_Callback(~, ~, ~)
% hObject    handle to ingap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ingap as text
%        str2double(get(hObject,'String')) returns contents of ingap as a double


% --- Executes during object creation, after setting all properties.
function ingap_CreateFcn(hObject, ~, ~)
% hObject    handle to ingap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%________AXON AFTER GAP____________________________________________________
function aftergap_Callback(~, ~, ~)
% hObject    handle to aftergap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of aftergap as text
%        str2double(get(hObject,'String')) returns contents of aftergap as a double


% --- Executes during object creation, after setting all properties.
function aftergap_CreateFcn(hObject, ~, ~)
% hObject    handle to aftergap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%________RESULTS 1 MEDIAN__________________________________________________
function median1_Callback(~, ~, ~)
% hObject    handle to median1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of median1 as text
%        str2double(get(hObject,'String')) returns contents of median1 as a double


% --- Executes during object creation, after setting all properties.
function median1_CreateFcn(hObject, ~, ~)
% hObject    handle to median1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%_______RESULTS 2 MEDIAN___________________________________________________
function median2_Callback(~, ~, ~)
% hObject    handle to median2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of median2 as text
%        str2double(get(hObject,'String')) returns contents of median2 as a double


% --- Executes during object creation, after setting all properties.
function median2_CreateFcn(hObject, ~, ~)
% hObject    handle to median2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%________RESULTS 3 MEDIAN__________________________________________________
function median3_Callback(~, ~, ~)
% hObject    handle to median3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of median3 as text
%        str2double(get(hObject,'String')) returns contents of median3 as a double


% --- Executes during object creation, after setting all properties.
function median3_CreateFcn(hObject, ~, ~)
% hObject    handle to median3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%_______________RESET______________________________________________________
% --- Executes on button press in reset.
function reset_Callback(~, ~, ~)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


evalin('base','clear all'); 
close(gcbf)
gcanalyze_v2_0
