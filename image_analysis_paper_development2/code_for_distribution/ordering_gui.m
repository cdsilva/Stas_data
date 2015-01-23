function ordering_gui
% ordering_gui shows options for registering and ordering images



%  Create and then hide the GUI as it is being constructed.
f = figure('Visible','off','Position',[360,500,450,285]);

% Construct the components.
crop    = uicontrol('Style','checkbox',...
             'String','Crop?','Position',[10,220,70,25],...
             'Callback',@surfbutton_Callback);
         
filename_label = uicontrol('Style','text',...
    'String','File Name','position',[10,180,70,25]);
         
filename = uicontrol('Style','edit',...
    'String','','position',[10,160,70,25]);

dir_label = uicontrol('Style','text',...
    'String','Directory','position',[10,130,70,25]);
         
dir = uicontrol('Style','edit',...
    'String','','position',[10,110,70,25]);

weight_labels = uicontrol('Style','text',...
    'String','Channel weights','position',[10,90,70,25]);
         
weight_red = uicontrol('Style','slider',...
    'String','Red','position',[100,90,70,25]);

weight_green = uicontrol('Style','slider',...
    'String','Green','position',[100,70,70,25]);

weight_blue = uicontrol('Style','slider',...
    'String','Blue','position',[100,50,70,25]);


% align([hsurf,hmesh,hcontour,htext,hpopup],'Center','None');
% 
% % Initialize the GUI.
% % Change units to normalized so components resize automatically.
% f.Units = 'normalized';
% hsurf.Units = 'normalized';
% hmesh.Units = 'normalized';
% hcontour.Units = 'normalized';
% htext.Units = 'normalized';
% hpopup.Units = 'normalized';
% 
% % Generate the data to plot.
% peaks_data = peaks(35);
% membrane_data = membrane;
% [x,y] = meshgrid(-8:.5:8);
% r = sqrt(x.^2+y.^2) + eps;
% sinc_data = sin(r)./r;
% 
% % Create a plot in the axes.
% current_data = peaks_data;
% surf(current_data);
% 
% % Assign the GUI a name to appear in the window title.
% f.Name = 'Simple GUI';
% 
% % Move the GUI to the center of the screen.
% % movegui(f,'center')

% Make the GUI visible.
set(f, 'Visible', 'on');

%  Pop-up menu callback. Read the pop-up menu Value property to
%  determine which item is currently displayed and make it the
%  current data. This callback automatically has access to 
%  current_data because this function is nested at a lower level.
   function popup_menu_Callback(source,eventdata) 
      % Determine the selected data set.
      str = get(source, 'String');
      val = get(source,'Value');
      % Set current data to the selected data set.
      switch str{val};
      case 'Peaks' % User selects Peaks.
         current_data = peaks_data;
      case 'Membrane' % User selects Membrane.
         current_data = membrane_data;
      case 'Sinc' % User selects Sinc.
         current_data = sinc_data;
      end
   end

  % Push button callbacks. Each callback plots current_data in the
  % specified plot type.

  function surfbutton_Callback(source,eventdata) 
  % Display surf plot of the currently selected data.
       surf(current_data);
  end

  function meshbutton_Callback(source,eventdata) 
  % Display mesh plot of the currently selected data.
       mesh(current_data);
  end

  function contourbutton_Callback(source,eventdata) 
  % Display contour plot of the currently selected data.
       contour(current_data);
  end
end