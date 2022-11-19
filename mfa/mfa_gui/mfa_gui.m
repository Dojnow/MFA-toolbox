 %'202007222033.22
function varargout = mfa_gui(varargin)

 % MultiFractal Analysis with Graphical User interface.
 %
 % Firstly, push the button "Update", then the button "Browse", then select a file with time series, then the buttons "Next" and "Finish". Use the mouse/touchpad for cursor positioning and choose appropriate time series. From the pop-up menu choose multifractal results and parameters - Fq(s), H(q), Tau(q), F(a). A detrending interval can be selected with the mouse/touchpad for cursor positioning until a singular spectrum without defects is obtained.

 % Copyright (C) 2020 Peter Dojnow

 % This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

 % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.


% MFA_GUI M-file for mfa_gui.fig
%      MFA_GUI, by itself, creates a new MFA_GUI or raises the existing
%      singleton*.
%
%      H = MFA_GUI returns the handle to a new MFA_GUI or the handle to
%      the existing singleton*.
%
%      MFA_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MFA_GUI.M with the given input arguments.
%
%      MFA_GUI('Property','Value',...) creates a new MFA_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mfa_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mfa_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mfa_gui

% Last Modified by GUIDE v2.5 22-Jul-2010 20:23:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mfa_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @mfa_gui_OutputFcn, ...
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

% --- Executes just before mfa_gui is made visible.
function mfa_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mfa_gui (see VARARGIN)

% Choose default command line output for mfa_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using mfa_gui.
if strcmp(get(hObject,'Visible'),'off')
    %plot(rand(5));
end

% UIWAIT makes mfa_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mfa_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;
global D Fsq Hq Tauq Fa I int intf q
popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1
        D = uiimport; D=D.(cell2mat(fieldnames(D)));
        plot(D);
        [x,y] = ginput(1);
        for i=1:size(D,2), z(i)=abs(D(round(x),i)-y); end
        I=find(min(z)==z);
    case 2
        [Fsq,int,q]=mdfast(D(:,I),-12:0.5:12,round(logspace(log10(10),log10(length(D)/4),32)),length(D),'single');
        loglog(Fsq)
        button = 'No';%(5)
        while strcmp(button,'No')%(5)
            loglog(int,Fsq)%,title([filename,num2str(j)])%(3) loglog(int,Dat.(Var{1})(:,:,j))%(2) loglog(int,squeeze(D(:,j,:)))%(1)
            [fitint,y] = ginput(1);% fitint(fitint<int(1)) = int(1);%(2) [Fitint(:,j,i),y]=ginput(2)%(1)
            [fitint(2),y] = ginput(1); fitint(fitint>int(end)) = int(end); fitint(fitint<int(1)) = int(1); %(2) %(6)
            fitint = sort(floor(fitint)); intf = int>=fitint(1) & int<=fitint(end);
            [a,f] = spect(tau(fit(Fsq(intf,:),int(intf),'log',1),q),q); plot(a,f)%(5.2)
            button = questdlg('Next?'); assert(~strcmp(button,'Cancel'),'Cancel')%(5)
        end%(5)
    case 3
        button = 'No';%(5)
        while strcmp(button,'No')%(5)
            loglog(int,Fsq)%,title([filename,num2str(j)])%(3) loglog(int,Dat.(Var{1})(:,:,j))%(2) loglog(int,squeeze(D(:,j,:)))%(1)
            [fitint,y] = ginput(1);% fitint(fitint<int(1)) = int(1);%(2) [Fitint(:,j,i),y]=ginput(2)%(1)
            [fitint(2),y] = ginput(1); fitint(fitint>int(end)) = int(end); fitint(fitint<int(1)) = int(1); %(2) %(6)
            fitint = sort(floor(fitint)); intf = int>=fitint(1) & int<=fitint(end);
            Hq = fit(Fsq(intf,:),int(intf),'log',1); plot(q,Hq,'-o')%(5.2)
            button = questdlg('Next?'); assert(~strcmp(button,'Cancel'),'Cancel')%(5)
        end%(5)
    case 4
        button = 'No';%(5)
        while strcmp(button,'No')%(5)
            loglog(int,Fsq)%,title([filename,num2str(j)])%(3) loglog(int,Dat.(Var{1})(:,:,j))%(2) loglog(int,squeeze(D(:,j,:)))%(1)
            [fitint,y] = ginput(1);% fitint(fitint<int(1)) = int(1);%(2) [Fitint(:,j,i),y]=ginput(2)%(1)
            [fitint(2),y] = ginput(1); fitint(fitint>int(end)) = int(end); fitint(fitint<int(1)) = int(1); %(2) %(6)
            fitint = sort(floor(fitint)); intf = int>=fitint(1) & int<=fitint(end);
            Tauq = tau(fit(Fsq(intf,:),int(intf),'log',1),q); plot(q,Tauq,'-o')%(5.2)
            button = questdlg('Next?'); assert(~strcmp(button,'Cancel'),'Cancel')%(5)
        end%(5)
    case 5
        button = 'No';%(5)
        while strcmp(button,'No')%(5)
            loglog(int,Fsq)%,title([filename,num2str(j)])%(3) loglog(int,Dat.(Var{1})(:,:,j))%(2) loglog(int,squeeze(D(:,j,:)))%(1)
            [fitint,y] = ginput(1);% fitint(fitint<int(1)) = int(1);%(2) [Fitint(:,j,i),y]=ginput(2)%(1)
            [fitint(2),y] = ginput(1); fitint(fitint>int(end)) = int(end); fitint(fitint<int(1)) = int(1); %(2) %(6)
            fitint = sort(floor(fitint)); intf = int>=fitint(1) & int<=fitint(end);
            [a,f] = spect(tau(fit(Fsq(intf,:),int(intf),'log',1),q),q); plot(a,f,'-o')%(5.2)
            button = questdlg('Next?'); assert(~strcmp(button,'Cancel'),'Cancel')%(5)
        end%(5)
end


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


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'Time Series', 'Fq(s)', 'H(q)', 'Tau(q)', 'F(a)'});
