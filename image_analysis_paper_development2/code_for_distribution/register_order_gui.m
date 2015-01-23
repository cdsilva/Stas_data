function varargout = register_order_gui(varargin)
% REGISTER_ORDER_GUI MATLAB code for register_order_gui.fig
%      REGISTER_ORDER_GUI, by itself, creates a new REGISTER_ORDER_GUI or raises the existing
%      singleton*.
%
%      H = REGISTER_ORDER_GUI returns the handle to a new REGISTER_ORDER_GUI or the handle to
%      the existing singleton*.
%
%      REGISTER_ORDER_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REGISTER_ORDER_GUI.M with the given input arguments.
%
%      REGISTER_ORDER_GUI('Property','Value',...) creates a new REGISTER_ORDER_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
% %      applied to the GUI before register_order_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to register_order_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help register_order_gui

% Last Modified by GUIDE v2.5 22-Jan-2015 15:43:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @register_order_gui_OpeningFcn, ...
    'gui_OutputFcn',  @register_order_gui_OutputFcn, ...
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


% --- Executes just before register_order_gui is made visible.
function register_order_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to register_order_gui (see VARARGIN)

% Choose default command line output for register_order_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes register_order_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = register_order_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function output_image_dir_Callback(hObject, eventdata, handles)
% hObject    handle to output_image_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output_image_dir as text
%        str2double(get(hObject,'String')) returns contents of output_image_dir as a double

handles.image_dir = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function output_image_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_image_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_images.
function save_images_Callback(hObject, eventdata, handles)
% hObject    handle to save_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function number_images_Callback(hObject, eventdata, handles)
% hObject    handle to number_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_images as text
%        str2double(get(hObject,'String')) returns contents of number_images as a double

nimages = str2double(get(hObject,'String'));
handles.nimages = nimages;

handles.subplot_dim1 = round(sqrt(nimages));
handles.subplot_dim2 = ceil(nimages/handles.subplot_dim1);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function number_images_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function image_extension_Callback(hObject, eventdata, handles)
% hObject    handle to image_extension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of image_extension as text
%        str2double(get(hObject,'String')) returns contents of image_extension as a double

handles.image_ext = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function image_extension_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_extension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in read_images.
function read_images_Callback(hObject, eventdata, handles)
% hObject    handle to read_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.images_raw, handles.nchannels] = read_images(handles.image_dir, handles.image_name, handles.image_ext, handles.stack_name, handles.nimages, handles.nstack, handles.npixels, handles.dim);

if handles.nchannels == 1
    
    for i=2:3
        
        set(handles.weight_slider(i), 'enable','off')
        set(handles.blur_slider(i), 'enable','off')
        set(handles.weight_numbers(i), 'enable','off')
        set(handles.blur_numbers(i), 'enable','off')
        set(handles.normalize_checkbox(i), 'enable','off')
        set(handles.mean_center_checkbox(i), 'enable','off')
        set(handles.normalize_checkbox(i), 'value',0)
        set(handles.mean_center_checkbox(i), 'value',0)
        set(handles.color_label(i), 'enable','off')
    end
end


guidata(hObject,handles);


% --- Executes on selection change in image_dim.
function image_dim_Callback(hObject, eventdata, handles)
% hObject    handle to image_dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns image_dim contents as cell array
%        contents{get(hObject,'Value')} returns selected item from image_dim

str = get(hObject, 'String');
val = get(hObject,'Value');

% Set current data to the selected data set.
switch str{val};
    case '2D'
        handles.dim = 2;
        set(handles.stack_prefix, 'enable','off')
        set(handles.number_stack, 'enable','off')
        set(handles.stack_prefix, 'value','')
        set(handles.number_stack, 'value','')
        set(handles.register_button, 'enable', 'on')
        set(handles.order_button, 'enable', 'on')
        set(handles.register_order_button, 'enable', 'on')
        set(handles.register_order_zstack_button, 'enable', 'off')
    case '3D'
        handles.dim = 3;
        set(handles.stack_prefix, 'enable','on')
        set(handles.number_stack, 'enable','on')
        set(handles.register_button, 'enable', 'off')
        set(handles.order_button, 'enable', 'on')
        set(handles.register_order_button, 'enable', 'off')
        set(handles.register_order_zstack_button, 'enable', 'on')
end
% Save the handles structure.
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function image_dim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.dim = 2;
guidata(hObject,handles)


function image_dir_Callback(hObject, eventdata, handles)
% hObject    handle to image_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of image_dir as text
%        str2double(get(hObject,'String')) returns contents of image_dir as a double

handles.image_dir = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function image_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function image_prefix_Callback(hObject, eventdata, handles)
% hObject    handle to image_prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of image_prefix as text
%        str2double(get(hObject,'String')) returns contents of image_prefix as a double

handles.image_name = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function image_prefix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function number_stack_Callback(hObject, eventdata, handles)
% hObject    handle to number_stack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_stack as text
%        str2double(get(hObject,'String')) returns contents of number_stack as a double

handles.nstack = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function number_stack_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_stack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'enable','off');
handles.number_stack = hObject;
handles.nstack = 0;
guidata(hObject,handles);



function stack_prefix_Callback(hObject, eventdata, handles)
% hObject    handle to stack_prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stack_prefix as text
%        str2double(get(hObject,'String')) returns contents of stack_prefix as a double

handles.stack_name = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function stack_prefix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stack_prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'enable','off');
handles.stack_prefix = hObject;
handles.stack_name = '';
guidata(hObject,handles);



function number_pixels_Callback(hObject, eventdata, handles)
% hObject    handle to number_pixels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_pixels as text
%        str2double(get(hObject,'String')) returns contents of number_pixels as a double

handles.npixels = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function number_pixels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_pixels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in show_images1.
function show_images1_Callback(hObject, eventdata, handles)
% hObject    handle to show_images1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_images(handles.images_raw, handles.dim, handles.subplot_dim1, handles.subplot_dim2)

% --- Executes on button press in apply_image_functions.
function apply_image_functions_Callback(hObject, eventdata, handles)
% hObject    handle to apply_image_functions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.images = apply_image_functions(handles.images_raw, handles.dim, handles.channel_weight, handles.channel_blur, handles.channel_normalize, handles.channel_mean_center);

guidata(hObject, handles);


% --- Executes on button press in show_images2.
function show_images2_Callback(hObject, eventdata, handles)
% hObject    handle to show_images2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_images(handles.images, handles.dim, handles.subplot_dim1, handles.subplot_dim2)

function kernel_scale_Callback(hObject, eventdata, handles)
% hObject    handle to kernel_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kernel_scale as text
%        str2double(get(hObject,'String')) returns contents of kernel_scale as a double

handles.eps_scale = str2double(get(hObject,'String')) ;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function kernel_scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kernel_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', '1');
handles.eps_scale = 1;
guidata(hObject, handles);


% --- Executes on button press in register_order_vdm.
function register_order_vdm_Callback(hObject, eventdata, handles)
% hObject    handle to register_order_vdm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.recalc_pairwise_rotations
    show_waitbar = true;
    [handles.R, handles.W] = compute_pairwise_alignments(handles.images, handles.nrot, show_waitbar);
    handles.recalc_pairwise_rotations = false;
end

% compute optimal rotations + embedding coordinates using vector diffusion maps
ncomps = 1;
[R_opt, embed_coord, D2] = vdm(handles.R, handles.W, handles.eps_scale, ncomps);

% register images using optimal rotations
images_registered = register_all_images(handles.images, R_opt);

[~, I] = sort(embed_coord);

if ndims(images_registered) == 3
    handles.images_analyzed = images_registered(:,:,I);
elseif ndims(images_registered) == 4
    handles.images_analyzed = images_registered(:,:,:,I);
elseif ndims(images_registered) == 5
    handles.images_analyzed = images_registered(:,:,:,I);
end

guidata(hObject, handles)



function number_rotations_Callback(hObject, eventdata, handles)
% hObject    handle to number_rotations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_rotations as text
%        str2double(get(hObject,'String')) returns contents of number_rotations as a double
handles.nrot = str2double(get(hObject,'String'));
handles.recalc_pairwise_rotations = true;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function number_rotations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_rotations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String','36');
handles.nrot = 36;
handles.recalc_pairwise_rotations = true;
guidata(hObject, handles)


% --- Executes on button press in order.
function order_Callback(hObject, eventdata, handles)
% hObject    handle to order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


W = squareform(pdist(reshape(double(images),[], handles.nimages)')).^2;

ncomps = 1;
[embed_coord, D2] = dmaps(W, median(W(:))*handles.eps_scale, ncomps);

[~, I] = sort(embed_coord);

if ndims(images) == 3
    handles.images_analyzed = handles.images(:,:,I);
elseif ndims(images) == 4
    handles.images_analyzed = handles.images(:,:,:,I);
elseif ndims(images) == 5
    handles.images_analyzed = handles.images(:,:,:,I);
end

guidata(hObject, handles)

% --- Executes on button press in register.
function register_Callback(hObject, eventdata, handles)
% hObject    handle to register (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.recalc_pairwise_rotations
    show_waitbar = true;
    [handles.R, handles.W] = compute_pairwise_alignments(handles.images, handles.nrot, show_waitbar);
    handles.recalc_pairwise_rotations = false;
end

% compute optimal rotations + embedding coordinates using vector diffusion maps
ncomps = 1;
[R_opt, embed_coord, D2] = vdm(handles.R, handles.W, handles.eps_scale, ncomps);

% register images using optimal rotations
handles.images_analyzed = register_all_images(handles.images, R_opt);

guidata(hObject, handles)


% --- Executes on button press in register_order_zstacks.
function register_order_zstacks_Callback(hObject, eventdata, handles)
% hObject    handle to register_order_zstacks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if handles.nchannels > 1
    max_proj_images = squeeze(max(handles.images, [], 4));
else
    max_proj_images = squeeze(max(handles.images, [], 3));
end

if handles.recalc_pairwise_rotations
    show_waitbar = true;
    [handles.R, handles.W] = compute_pairwise_alignments(max_proj_images, handles.nrot, show_waitbar);
    handles.recalc_pairwise_rotations = false;
end

% compute optimal rotations + embedding coordinates using vector diffusion maps
ncomps = 1;
[R_opt, embed_coord, D2] = vdm(handles.R, handles.W, handles.eps_scale, ncomps);
images_registered = register_all_images(handles.images, R_opt);

W = squareform(pdist(reshape(double(images_registered), [], handles.nimages)')).^2;
eps = median(W(:))*handles.eps_scale;
[V, D] = dmaps(W, eps, 10);

[~, I] = sort(V(:,2));

if ndims(handles.images) == 4
    handles.images_analyzed = images_registered(:,:,:,I);
elseif ndims(handles.images) == 5
    handles.images_analyzed = images_registered(:,:,:,:,I);
end

guidata(hObject, handles)


% --- Executes on slider movement.
function red_weight_Callback(hObject, eventdata, handles)
% hObject    handle to red_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.channel_weight(1) = get(hObject,'Value');
set(handles.weight_numbers(1), 'String', sprintf('%2.2f',get(hObject,'Value')));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function red_weight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject, 'Value', 1)
handles.channel_weight(1) = 1;
handles.weight_slider(1) = hObject;
guidata(hObject, handles);


% --- Executes on slider movement.
function green_weight_Callback(hObject, eventdata, handles)
% hObject    handle to green_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.channel_weight(2) = get(hObject,'Value');
set(handles.weight_numbers(2), 'String', sprintf('%2.2f',get(hObject,'Value')));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function green_weight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject, 'Value', 1)
handles.channel_weight(2) = 1;
handles.weight_slider(2) = hObject;
guidata(hObject, handles);

% --- Executes on slider movement.
function blue_weight_Callback(hObject, eventdata, handles)
% hObject    handle to blue_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.channel_weight(3) = get(hObject,'Value');
set(handles.weight_numbers(3), 'String', sprintf('%2.2f',get(hObject,'Value')));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function blue_weight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject, 'Value', 1)
handles.channel_weight(3) = 1;

handles.weight_slider(3) = hObject;
guidata(hObject, handles);

% --- Executes on slider movement.
function red_blur_Callback(hObject, eventdata, handles)
% hObject    handle to red_blur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.channel_blur(1) = get(hObject,'Value');
set(handles.blur_numbers(1), 'String', sprintf('%2.0f%%',100*get(hObject,'Value')));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function red_blur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_blur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject, 'Value', 0)
handles.channel_blur(1) = 0;

handles.blur_slider(1) = hObject;
guidata(hObject, handles);

% --- Executes on slider movement.
function green_blur_Callback(hObject, eventdata, handles)
% hObject    handle to green_blur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.channel_blur(2) = get(hObject,'Value');
set(handles.blur_numbers(2), 'String', sprintf('%2.0f%%',100*get(hObject,'Value')));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function green_blur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green_blur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject, 'Value', 0)
handles.channel_blur(2) = 0;

handles.blur_slider(2) = hObject;
guidata(hObject, handles);

% --- Executes on slider movement.
function blue_blur_Callback(hObject, eventdata, handles)
% hObject    handle to blue_blur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.channel_blur(3) = get(hObject,'Value');
set(handles.blur_numbers(3), 'String', sprintf('%2.0f%%',100*get(hObject,'Value')));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function blue_blur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue_blur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject, 'Value', 0)
handles.channel_blur(3) = 0;

handles.blur_slider(3) = hObject;
guidata(hObject, handles);

% --- Executes on button press in red_normalize.
function red_normalize_Callback(hObject, eventdata, handles)
% hObject    handle to red_normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of red_normalize

handles.channel_normalize(1) = get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in green_normalize.
function green_normalize_Callback(hObject, eventdata, handles)
% hObject    handle to green_normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of green_normalize

handles.channel_normalize(2) = get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in blue_normalize.
function blue_normalize_Callback(hObject, eventdata, handles)
% hObject    handle to blue_normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of blue_normalize

handles.channel_normalize(3) = get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in red_mean_center.
function red_mean_center_Callback(hObject, eventdata, handles)
% hObject    handle to red_mean_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of red_mean_center

handles.channel_mean_center(1) = get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in green_mean_center.
function green_mean_center_Callback(hObject, eventdata, handles)
% hObject    handle to green_mean_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of green_mean_center

handles.channel_mean_center(2) = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in blue_mean_center.
function blue_mean_center_Callback(hObject, eventdata, handles)
% hObject    handle to blue_mean_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of blue_mean_center

handles.channel_mean_center(3) = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function red_weight_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_weight_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.weight_numbers(1) = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function green_weight_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green_weight_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.weight_numbers(2) = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function blue_weight_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue_weight_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.weight_numbers(3) = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function red_blur_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_blur_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.blur_numbers(1) = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function green_blur_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green_blur_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.blur_numbers(2) = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function blue_blur_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue_blur_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.blur_numbers(3) = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function red_normalize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.normalize_checkbox(1) = hObject;
set(hObject, 'Value',0)
handles.channel_normalize(1) = false;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function green_normalize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green_normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.normalize_checkbox(2) = hObject;
set(hObject, 'Value',0)
handles.channel_normalize(2) = false;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function blue_normalize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue_normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.normalize_checkbox(3) = hObject;
set(hObject, 'Value',0)
handles.channel_normalize(3) = false;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function red_mean_center_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_mean_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.mean_center_checkbox(1) = hObject;
set(hObject, 'Value',1)
handles.channel_mean_center(1) = true;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function green_mean_center_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green_mean_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.mean_center_checkbox(2) = hObject;
set(hObject, 'Value',1)
handles.channel_mean_center(2) = true;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function blue_mean_center_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue_mean_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.mean_center_checkbox(3) = hObject;
set(hObject, 'Value',1)
handles.channel_mean_center(3) = true;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function green_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.color_label(2) = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function blue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


handles.color_label(3) = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function register_CreateFcn(hObject, eventdata, handles)
% hObject    handle to register (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.register_button = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - hand

handles.order_button = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function register_order_vdm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to register_order_vdm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.register_order_button = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function register_order_zstacks_CreateFcn(hObject, eventdata, handles)
% hObject    handle to register_order_zstacks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.register_order_zstack_button = hObject;
set(hObject,'enable','off');
guidata(hObject, handles);


% --- Executes on button press in show_images3.
function show_images3_Callback(hObject, eventdata, handles)
% hObject    handle to show_images3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_images(handles.images_analyzed, handles.dim, handles.subplot_dim1, handles.subplot_dim2)
