function varargout = firstpanel(varargin)
% FIRSTPANEL MATLAB code for firstpanel.fig
%      FIRSTPANEL, by itself, creates a new FIRSTPANEL or raises the existing
%      singleton*.
%
%      H = FIRSTPANEL returns the handle to a new FIRSTPANEL or the handle to
%      the existing singleton*.
%
%      FIRSTPANEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIRSTPANEL.M with the given input arguments.
%
%      FIRSTPANEL('Property','Value',...) creates a new FIRSTPANEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before firstpanel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to firstpanel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help firstpanel

% Last Modified by GUIDE v2.5 11-Dec-2017 14:38:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @firstpanel_OpeningFcn, ...
                   'gui_OutputFcn',  @firstpanel_OutputFcn, ...
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


% --- Executes just before firstpanel is made visible.
function firstpanel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to firstpanel (see VARARGIN)

% Choose default command line output for firstpanel
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes firstpanel wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = firstpanel_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn_wavefield.
function btn_wavefield_Callback(hObject, eventdata, handles)
% hObject    handle to btn_wavefield (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%	template_beamform.m
%	Ultrasound beamforming project template
%
%	This script requires the variables
%	Data [ntime, nelem]	RF signals from each array element
%	ed_f0			Transducer center frequency (MHz)
%	ed_fs			sampling frequency (MHz)
%	ed_c			speed of sound (mm/usec)
%	dx			transducer element spacing (mm)
%
% Get needed variables
%
load data06;
f0 = str2double(get(handles.ed_f0, 'String'));  %4 MHz
fs = str2double(get(handles.ed_fs, 'String'));	%16 MHz
c =  str2double(get(handles.ed_c,  'String'));  %1.54 mm/us
db = str2double(get(handles.ed_db, 'String'));
% 
% Output Array Parameters
%
% Determine the array spacing dx in mm
global lambda;
lambda = c/f0;
dx = lambda/2;
deltat=1/fs;

[ntime, nelem] = size(Data);		% # time samples, # array elements
%fprintf('f0=%g MHz, deltat=%g usec, dx=%g mm\n', f0, deltat, dx)
%fprintf('# of Time Samples=%g,  # of Array Elements=%g\n',ntime,nelem)

%
% --> QUESTION A <--
% Make a btn_wavefield plot of the raw Data
% Comment this out while you debug other parts of your program
%
% showimage(Data, -4)
% disp 'hit key', pause;
showimage(Data,-4,db)
disp 'hit key';  
return;



function ed_db_Callback(hObject, eventdata, handles)
% hObject    handle to ed_db (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_db as text
%        str2double(get(hObject,'String')) returns contents of ed_db as a double


% --- Executes during object creation, after setting all properties.
function ed_db_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_db (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ed_f0_Callback(hObject, eventdata, handles)
% hObject    handle to ed_f0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_f0 as text
%        str2double(get(hObject,'String')) returns contents of ed_f0 as a double


% --- Executes during object creation, after setting all properties.
function ed_f0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_f0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_fs_Callback(hObject, eventdata, handles)
% hObject    handle to ed_fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_fs as text
%        str2double(get(hObject,'String')) returns contents of ed_fs as a double


% --- Executes during object creation, after setting all properties.
function ed_fs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_c_Callback(hObject, eventdata, handles)
% hObject    handle to ed_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_c as text
%        str2double(get(hObject,'String')) returns contents of ed_c as a double


% --- Executes during object creation, after setting all properties.
function ed_c_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_calc_beam.
function btn_calc_beam_Callback(hObject, eventdata, handles)
% hObject    handle to btn_calc_beam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --> QUESTION B <--
%   Compute the number of total beams and btn_calc_beam spacing
%   I used variables called
%	txt_nbeam = number of beams
%	sin_theta = vector of btn_calc_beam positions

nelem = str2double(get(handles.elem, 'String'));
delta_sin_theta = sin(1/nelem);
nbeam0 = (((sqrt(2)/2)-(-(sqrt(2)/2))))/(delta_sin_theta);
nbeam = nbeam0/2;
fprintf('nbeam = %g\n', nbeam); %(((sqrt(2)/2)-(-(sqrt(2)/2))))/(1/nspace);
sin_theta = linspace(sin(-pi/4),sin(pi/4),nbeam);
set (handles.txt_nbeam, 'String', nbeam);
return;

% modulus = mod(nbeam,2);
% nbeam = nbeam + (modulus == 0);
% fprintf('nbeam = %g\n', nbeam)
% linspace(sin(-pi/4),sin(pi/4),nbeam);

function elem_Callback(hObject, eventdata, handles)
% hObject    handle to elem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of elem as text
%        str2double(get(hObject,'String')) returns contents of elem as a double


% --- Executes during object creation, after setting all properties.
function elem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_nbeam_Callback(hObject, eventdata, handles)
% hObject    handle to txt_nbeam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_nbeam as text
%        str2double(get(hObject,'String')) returns contents of txt_nbeam as a double


% --- Executes during object creation, after setting all properties.
function txt_nbeam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_nbeam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_clean_beam.
function btn_clean_beam_Callback(hObject, eventdata, handles)
% hObject    handle to btn_clean_beam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set (handles.elem, 'String', '');
set (handles.txt_nbeam, 'String', '');


% --- Executes on button press in btn_clean_wavefield.
function btn_clean_wavefield_Callback(hObject, eventdata, handles)
% hObject    handle to btn_clean_wavefield (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set (handles.ed_f0, 'String', '');
set (handles.ed_fs, 'String', '');
set (handles.ed_c,  'String', '');
set (handles.ed_db, 'String', '');



% --- Executes on button press in btn_calc_delay.
function btn_calc_delay_Callback(hObject, eventdata, handles)
% hObject    handle to btn_calc_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% h = waitbar(0,'Initializing waitbar...');
% waitbar(0.5,h,'Halfway there...');
% perc = 75;
% waitbar(perc/100,h,sprintf('%d%% along...',perc));
% close(h);

f0 = str2double(get(handles.ed_f0, 'String'));  %4 MHz
fs = str2double(get(handles.ed_fs, 'String'));	%16 MHz
c =  str2double(get(handles.ed_c,  'String'));  %1.54 mm/us
% db = str2double(get(handles.ed_db, 'String'));
nelem = 65;
delta_sin_theta = sin(1/nelem);
nbeam0 = (((sqrt(2)/2)-(-(sqrt(2)/2))))/(delta_sin_theta);
nbeam = ceil(nbeam0/2);
sin_theta = linspace(sin(-pi/4),sin(pi/4),nbeam);
lambda = c/f0;
dx = lambda/2;
dr = lambda/2;% r = ct/2, ct = lambda
r=linspace(dr,60,60/dr);% r = 60 mm (beam length)
xn = linspace(-32*dx,32*dx,65);% element position
nrsamp = round(length(r)); %number of sample space along r
ir = 0;
ib = 0;
ie = 0;

tdelay = zeros(nbeam, nrsamp, nelem);
for ib = 1:nbeam
    for ir = 1:nrsamp
        for ie = 1:nelem % ie = 1 :65;
            tdelay(ib,ir,ie) = floor (fs*(-xn(ie)*sin_theta(ib)/c+(xn(ie)^2)*(1-sin_theta(ib)^2)/(2*c*r(ir))));
        end
    end
end
set (handles.txt_show_delay, 'String', 'CALCULATED!');
set (handles.txt_ib, 'String', ib);
set (handles.txt_ir, 'String', ir);
set (handles.txt_ie, 'String', ie);
% set (handles.tab_data, 'Data', tdelay);
% set (handles.txt_show_delay, 'String', tdelay );

%========================= SOLITARY REFLECTORS ==========================

%%%% r=8mm
load data06;
b1 = Data(:,1)';% size(b1) = 1 65
delay = ceil (tdelay (ceil (nbeam/2), ceil(8/dr),1)); %nbeam/2 means theta0,
    if delay >= 0
        b1 = [b1(delay+1 : end) zeros(1, delay)];
    else
        b1 = [zeros(1, -delay) b1(1:end+delay)];
    end
    axes(handles.ax1)
    plot(b1);
    %xlabel('Frequency [Hz]''color','w');
    %ylabel('Normalized Amplitude''color','w');
    title('r = 8mm, element 01','color','k');

    %figure, plot(b1);


b33 = Data(:,33)';
delay = ceil(tdelay(ceil(nbeam/2), ceil(8/dr),33));
    if delay >= 0
        b33=[b33(delay+1:end) zeros(1,delay)];
    else
        b33 = [zeros(1,-delay) b33(1:end+delay)];
    end
    axes(handles.ax2)
    plot(b33,'r');
    title('r = 15mm, element 01','color','k');
    %figure, plot(b33);


b65 = Data(:,65)';
delay = ceil(tdelay(ceil(nbeam/2),ceil(8/dr),65));
    if delay >=0
        b65 = [b65(delay+1:end) zeros(1,delay)];
    else
        b65 = [zeros(1,-delay) b65(1:end+delay)];
        end
    axes(handles.ax3)
    plot(b65);
    title('r = 8mm, element 33','color','k');
    %figure,plot(b65);


%%%% r=15mm
b1 = Data(:,1)';
delay = ceil (tdelay (ceil (nbeam/2), ceil(15/dr),1));%nbeam/2 means theta0,
    if delay >= 0
        b1 = [b1(delay+1 : end) zeros(1, delay)];
    else
        b1 = [zeros(1, -delay) b1(1:end+delay)];
    end
    axes(handles.ax4)
    plot(b1,'r');
    title('r = 15mm, element 33','color','k');
    %figure, plot(b1);


b33 = Data(:,33)';
delay = ceil(tdelay(ceil(nbeam/2), ceil(15/dr),33));
    if delay >= 0
        b33=[b33(delay+1:end) zeros(1,delay)];
    else
        b33 = [zeros(1,-delay) b33(1:end+delay)];
    end 
    axes(handles.ax5)
    plot(b33);
    title('r = 8mm, element 65','color','k');
    %figure, plot(b33);


b65 = Data(:,65)';
delay = ceil(tdelay(ceil(nbeam/2),ceil(15/dr),65));
    if delay >=0
        b65 = [b65(delay+1:end) zeros(1,delay)];
    else
        b65 =[ zeros(1,-delay) b65(1:end+delay)];
    end 
    axes(handles.ax6)
    plot(b65,'r');
    title('r = 15mm, element 65','color','k');
    %figure,plot(b65);
    %isBusy = CmdWinTool('isBusy');

% beam theta 0 formed by summing across all elements


[ntime, nelem] = size(Data);	
beam24 = zeros(1,nrsamp);
for ir = 1:nrsamp
    tempdata = zeros(ntime,1);
    for ie = 1:nelem
    delay = tdelay(24,ir,ie);
        if delay >= 0
    tempdata = tempdata + [Data(delay+1:end, ie);zeros(delay,1)];
        else
    tempdata = tempdata + [zeros(-delay,1);Data(1:end+delay,ie)];
    end
end
    tempdata = tempdata/nelem;
    beam24(ir) = tempdata(ceil(ir/nrsamp*1500));
end


axes(handles.ax20)
plot(beam24,'k');
title('Beam 24 (theta 0)','color','k');

msgbox('STEP 01 | Please, wait for a while. We are calculating...');
        

% create data for the r-sintheta
% buffer using all receiving elements for all beams
% r-sin(theta) data buffer after coherent sum across all elements

rsdata = zeros(nrsamp, nbeam);
for ib = 1:1:nbeam
    fprintf('Beam %d of %d\n', ib, nbeam);
        for i = 1:1:nrsamp
            tempdata = zeros(ntime,1);
        for ie = 1:1:nelem
            delay = tdelay(ib,i,ie);
        if delay >= 0
            tempdata = tempdata + [Data(delay+1:end,ie);zeros(delay,1)];
        else
            tempdata = tempdata + [zeros(-delay,1);Data(1:end+delay,ie)];
    end
end
tempdata = tempdata / nelem;
rsdata(i,ib) = tempdata(ceil(i/nrsamp*1500));
end
end

showimage(abs(rsdata),-4)
showimage(abs(rsdata),-4,40)

msgbox('STEP 02 | Please, wait for a while. We are calculating...');
                                    % END %
%============================ SOLITARY REFLECTORS ========================

databb = baseband(Data, f0/fs);
% Plot the signals for one element (Element 33) before and after baseband
% conversion
% Polot both I channel and Q channel after baseband conversion
axes(handles.axes16)
plot(Data(:,33),'k');

axes(handles.axes17)
plot(real(databb(:,33)),'k');

axes(handles.axes18)
plot(-1i*imag(databb(:,33)),'k');

axes(handles.axes19)
plot(abs(databb(:,33)),'k');

showimage(abs(databb), -4,40)

msgbox('STEP 03 | Please, wait for a while. We are calculating...');

%----------------------------- QUESTION G -------------------------

% ---------- BEGINNING OF WAIT BAR -----------

% h = waitbar(0,'Please wait...',...
%                'Name','Calculating time Delay',...
%                'CreateCancelBtn',...
%                'setappdata(gcbf,''cancel_callback'',1)');
% setappdata(h,'cancel_callback',0);
% steps = 10000;
% for step = 1:steps 
%     
%     %CALCULATIONS HERE
% 
% if getappdata(h,'cancel_callback')
%    break;
% end
%     waitbar(step / steps, h, sprintf('Calculating...%.2f%%',step/steps*100));
% end
% delete(h)
%---------------WAITAR----------------


rsdata2 = zeros(nrsamp, nbeam);
tdelaybb = zeros(nbeam, length(r), nelem);
for ib = 1:nbeam
    for ir = 1:nrsamp 
        for ie = 1:nelem % ie = 1 :65;
            tdelaybb(ib,ir,ie) = floor ...,
                (20*(-xn(ie)*sin_theta(ib)/c+(xn(ie)^2)*...,
                (1-sin_theta(ib)^2)/(2*c*r(ir))));
        end
    end
end


for ib = 1:1:nbeam
    fprintf('Beam %d of %d\n', ib, nbeam);
    for i = 1:1:nrsamp
        tempdata = zeros(ntime,1);
        for ie = 1:1:nelem
            delay = tdelaybb(ib,i,ie);
            if delay >= 0 
                tempdata = tempdata + [databb(delay+1:end,ie);zeros(delay,1)];
            else
                tempdata = tempdata + [zeros(-delay,1);databb(1:end+delay,ie)];
            end
        end
        tempdata = tempdata / nelem;
        rsdata2(i,ib) = tempdata(ceil(i/nrsamp*1500));
    end
end
showimage(abs(rsdata2),-4,40);
showimage(abs(rsdata2),-4,20);

rsdata3 = zeros(nrsamp, nbeam);
for ib = 1:1:nbeam %for odd channel
    fprintf('Beam %d of %d\n', ib, nbeam);
    for i = 1:1:nrsamp
        tempdata = zeros(ntime,1);
        for ie = 1:2:nelem
            delay = tdelaybb(ib,i,ie);
            if delay >= 0 
                tempdata = tempdata + [databb(delay+1:end,ie);zeros(delay,1)];
            else
                tempdata = tempdata + [zeros(-delay,1);databb(1:end+delay,ie)];
            end
        end
        tempdata = tempdata / nelem;
        rsdata3(i,ib) = tempdata(ceil(i/nrsamp*1500));
    end
end


%------------------------------END QUESTION G -------------------------

% ------------------------ BEGGIN SCAN CONVERTER -----------------------------

%	template_scancon.m
%	template script for converting from r-sin(theta) data to x-y image
%	must load r-sin(theta) data: rsdata

% --> QUESTION I <--
% Scan convert the r-sin(theta) buffer to produce a sector scan image.
% Use bilinear interpolation to compute the image values on the
% sector scan image grid.  Matlab's "interp2" function will help you
% do bilinear interpolation.

% compute values needed for interpolation
[R,S] = ndgrid(r,sin_theta);
xx = linspace(-30,30,512);
zz = linspace(0,60,512);
[Z,X] = ndgrid(zz,xx);
RI = (X.*X+Z.*Z).^0.5;
SI = X./RI;

% Create image w/ bilinear interpolation

image = interp2(S, R, abs(rsdata), SI, RI, 'bilinear');
t = find(isnan(image));
image(t) = zeros(size(t));

image2 = interp2(S, R, abs(rsdata2), SI, RI, 'bilinear');
t = find(isnan(image2));
image2(t) = zeros(size(t));

image3 = interp2(S, R, abs(rsdata3), SI, RI, 'bilinear');
t = find(isnan(image3));
image3(t) = zeros(size(t));

% --> QUESTION I <--
% Use two images on a logarithmic scale to answer this question:
% one on a 40dB scale, the other on a 20dB scale
showimage(image, 0, 20);		% Display 20 dB scale image
showimage(image, 0, 40);		% Display 40 dB scale image
showimage(image2, 0, 20);		% Display 20 dB scale image
showimage(image2, 0, 40);		% Display 40 dB scale image
showimage(image3, 0, 20);		% Display 20 dB scale image
showimage(image3, 0, 40);		% Display 40 dB scale image
%
% ALso, to answer this question, go back and look at the
% image PRIOR to scan conversion, especially in contrasting
% Artifacts for the Full Aperture and Decimated Aperture

msgbox('CALCULATION FINISHED!');
return;


% ------------------------ SCAN CONVERTER ----------------------------


function txt_show_delay_Callback(hObject, eventdata, handles)
% hObject    handle to txt_show_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_show_delay as text
%        str2double(get(hObject,'String')) returns contents of txt_show_delay as a double


% --- Executes during object creation, after setting all properties.
function txt_show_delay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_show_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_clean_delay.
function btn_clean_delay_Callback(hObject, eventdata, handles)
% hObject    handle to btn_clean_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function txt_ib_Callback(hObject, eventdata, handles)
% hObject    handle to txt_ib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_ib as text
%        str2double(get(hObject,'String')) returns contents of txt_ib as a double


% --- Executes during object creation, after setting all properties.
function txt_ib_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_ib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_ie_Callback(hObject, eventdata, handles)
% hObject    handle to txt_ie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_ie as text
%        str2double(get(hObject,'String')) returns contents of txt_ie as a double


% --- Executes during object creation, after setting all properties.
function txt_ie_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_ie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_ir_Callback(hObject, eventdata, handles)
% hObject    handle to txt_ir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_ir as text
%        str2double(get(hObject,'String')) returns contents of txt_ir as a double


% --- Executes during object creation, after setting all properties.
function txt_ir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_ir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function tab_data_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tab_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object deletion, before destroying properties.
function tab_data_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to tab_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btn_graphs.
function btn_graphs_Callback(hObject, eventdata, handles)
% hObject    handle to btn_graphs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function scr_Callback(hObject, eventdata, handles)
% hObject    handle to scr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function about_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('This Software was designed by Eduardo Martins as a final Ultrasound Class Project in the Biomedical Engineering Master Course at Yonsei University. Special thanks to PIAO JINDAN that helped a lot with the code. Contact: eduardo.aeronautics@outlook.com', 'About the Software');

% --------------------------------------------------------------------
function beam_theta_Callback(hObject, eventdata, handles)
% hObject    handle to ax20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
beam_theta


% --------------------------------------------------------------------
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button = questdlg('Are you sure do you want to quit Ultrasound Simulation?','Exit Dialog','Yes','No','No');
switch button
    case 'Yes'
        exit;
    case 'No'
        quit cancel;
end

% --- Executes during object creation, after setting all properties.
function ax1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ax3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object deletion, before destroying properties.
function ax1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to ax1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function ax3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ax3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate ax3


% --- Executes when uipanel1 is resized.
function uipanel1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_instructions_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_instructions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('In order to get the results of Calculated Beam Number and generate the Graphs and the Images, you need to fill the fields with text label in GREY. Filds with text label in YELLOW, will be calculated by the software. All the images will appear in your screen during the process. There are 3 STEPS that you need to wait the calculation. When the calculation is finished, a message will inform you.','INSTRUCTIONS');
