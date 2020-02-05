function varargout = beam_theta(varargin)
% BEAM_THETA MATLAB code for beam_theta.fig
%      BEAM_THETA, by itself, creates a new BEAM_THETA or raises the existing
%      singleton*.
%
%      H = BEAM_THETA returns the handle to a new BEAM_THETA or the handle to
%      the existing singleton*.
%
%      BEAM_THETA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BEAM_THETA.M with the given input arguments.
%
%      BEAM_THETA('Property','Value',...) creates a new BEAM_THETA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before beam_theta_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to beam_theta_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help beam_theta

% Last Modified by GUIDE v2.5 03-Dec-2017 13:03:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @beam_theta_OpeningFcn, ...
                   'gui_OutputFcn',  @beam_theta_OutputFcn, ...
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


% --- Executes just before beam_theta is made visible.
function beam_theta_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to beam_theta (see VARARGIN)

% Choose default command line output for beam_theta
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes beam_theta wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = beam_theta_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --- Executes on button press in btn_beam_theta.
function btn_beam_theta_Callback(hObject, eventdata, handles)
% hObject    handle to btn_beam_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% beam theta 0 formed by summing across all elements

%variables declaration
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
figure, plot(beam24);
title('Beam 24 (theta 0)');
        

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
msgbox('DONE!');

                                    % END %
%============================ SOLITARY REFLECTORS ========================
       
