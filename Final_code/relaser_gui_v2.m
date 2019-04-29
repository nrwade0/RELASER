%%  RELASER GUI
% Version 2.0 - March 15, 2019 NW
%             - Second try. First used tabs and became unnecessarily 
%             complicated.
%             - Built basic functionality options I/O file push buttons,
%             text, and run button.
%             - Some experimentation with specific datadump file selection
%             but no luck so far (just outputs all for now).
%             - Put in some kind of plotting system?
%
% Version 2.1 - March 20-27, 2019 NW
%             - Added variable list under inputs for convenience, should be
%             updated when program is fully functional.
%             - Some experimentation with a drop-down menu and plots, not
%             sure what exactly the best plots are. Looking at fig 8,9 from
%             "Spectroscopy and modeling of solid state lanthanide lasers:
%             Application..." by B. Walsh, April 2004. Possibly others from
%             Nash's version to make sure it works correctly.
%
% Version 2.2 - March 28-April 3, 2019 NW
%             - Added axes with different kinds of plots from drop-down
%             menu. First one is most important and should be replicated
%             from above paper for success.
%             - Added some manifold populations and phi over time as well
%             for the hell of it - seems good to know.
%             - How do these look compared to Fortran program?
%
% Version 2.3 - April 4-15 , 2019 NW
%             - Looking into data reduction techniques for some population
%             over time plots.
%             - Fixing current program.
%
% FUTURE IDEAS:
%             - Make variables and arguments from both (and how they
%             connect) more comprehensible. Might be worth jotting a
%             relaser gui variable list down, as well.
%
%
% Supress 'not used' error messages in callback functions
%#ok<*DEFNU,*INUSD>
%#ok<*INUSL> 
%% ------------------------------- MAIN GUI CONTROLS ------------------------%
function varargout = relaser_gui_v2(varargin)
    %  RELASER_GUI_V2 MATLAB code for relaser_gui_v2.fig
    %      RELASER_GUI_V2, by itself, creates a new RELASER_GUI_V2 or raises the existing
    %      singleton*.
    %
    %      H = RELASER_GUI_V2 returns the handle to a new RELASER_GUI_V2 or the handle to
    %      the existing singleton*.
    %
    %      RELASER_GUI_V2('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in RELASER_GUI_V2.M with the given input arguments.
    %
    %      RELASER_GUI_V2('Property','Value',...) creates a new RELASER_GUI_V2 or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before relaser_gui_v2_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to relaser_gui_v2_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help relaser_gui_v2

    % Last Modified by GUIDE v2.5 08-Apr-2019 08:55:41

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @relaser_gui_v2_OpeningFcn, ...
                       'gui_OutputFcn',  @relaser_gui_v2_OutputFcn, ...
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
end

% --- Executes just before relaser_gui_v2 is made visible.
function relaser_gui_v2_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to relaser_gui_v2 (see VARARGIN)

    % Choose default command line output for relaser_gui_v2
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes relaser_gui_v2 wait for user response (see UIRESUME)
    % uiwait(handles.figure1);
end


% --- Outputs from this function are returned to the command line.
function varargout = relaser_gui_v2_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;
end



%% ------------------------ OUTPUT FILE COMMANDS ---------------------------%
% Button to open popnum.txt file
function popnum_Callback(hObject, eventdata, handles)
    % hObject    handle to popnum (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    winopen('popnum.txt');

end


% Button to open perfnum.txt file
function perfnum_Callback(hObject, eventdata, handles)
    % hObject    handle to perfnum (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    winopen('perfnum.txt');

end


% Button to open a data dump file
function datadump_Callback(hObject, eventdata, handles)
    % hObject    handle to datadump (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Load in relaser data
    S = load('final_relaser.mat');
    
    % Load user choice from popup menu
    popup_sel_index = get(handles.datadump_popupmenu, 'Value');
    
    % Find choice the user made and open that file
    for i = 1:S.IMAX
        if popup_sel_index == i
            file_name = sprintf('datadump_%d.txt', i);
            winopen(file_name);
        end
    end
    
end


% Popup menu selection callback to select data dump file
% Accessed when a choice is made in the menu
function datadump_popupmenu_Callback(hObject, eventdata, handles)
    % hObject    handle to datadump_popupmenu (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns datadump_popupmenu contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from datadump_popupmenu
end


% Popup menu for data dump file
% Function is used when creating the file, i.e. numbers in the list
function datadump_popupmenu_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to datadump_popupmenu (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: popupmenu controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
    % Load in relaser data
    S = load('final_relaser.mat');
    
    dump_files = cell(1, S.IMAX);
    for i = 1:S.IMAX
        dump_files(1, i) = num2cell(i);
    end
    
    set(hObject, 'String', dump_files);
end



%% -------------------------- INPUT FILE COMMANDS --------------------------%
% Button to open the LASPAR.xlsx file
function laspar_file_Callback(hObject, eventdata, handles)
    % hObject    handle to laspar_file (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    winopen('laspar.xlsx');
    
end


% Button to open the SPECPAR.xlsx file
function specpar_file_Callback(hObject, eventdata, handles)
    % hObject    handle to specpar_file (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    winopen('specpar.xlsx');

end


% Button to open the Variable list text file
function variables_Callback(hObject, eventdata, handles)
    % hObject    handle to variables (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    winopen('variable_list.txt');
end



%% --------------------- EXECUTION AND VERSION CONTROL ----------------------%
% Runs the relaser program and returns to GUI (updates relaser data for
% plots and output files.
function execute_Callback(hObject, eventdata, handles)
    % hObject    handle to execute (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    relasermain
end

% Popup menu to open the previous version history from the Fortran file.
function Fversionhistory_Callback(hObject, eventdata, handles)
    % open version history text file
    winopen('272F_vhist.F');
end



%% ----------------------------- PLOTTING ----------------------------------%
% Popup menu to change the current plot.
% Accessed when a choice is made in the menu
function plot_popup_Callback(hObject, eventdata, handles)
    % hObject    handle to plot_popup (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns plot_popup contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from plot_popup

end


% Popup menu for plot choice
% Function is used when creating the file, i.e. plot names
function plot_popup_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to plot_popup (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: popupmenu controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
    set(hObject, 'String', {'Pump Energy vs. Laser Energy', 't vs. n5', ...
        't vs. n6', 't vs. n7', 't vs. n8', 't vs. phi', ...
        't vs FracUp*n7, FracLo*n8', 't vs FracUp*n7 - FracLo*n8'});
    
end


% Creates axes for the plotting area
% If a default plot is wanted, you can change that here.
% Plots themselves aren't made here, see "update_buton_Callback(...).
function main_plot_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to main_plot (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: place code in OpeningFcn to populate main_plot
    % load in relasermain.m data
    S = load('final_relaser.mat');
% 
%     % options = 1:S.IMAX;
    hold on
    xlim([0 ceil(S.penergy(S.IMAX))]);
    ylim([0 ceil(S.EABS_plot(S.IMAX))]);
    xlabel('Pump Energy');
    ylabel('Laser Energy');
    plot(S.penergy, S.EABS, 'b-');
    hold off

end


% Button to refresh the plot based on current popup menu choice
% Change the plot data here.
function update_button_Callback(hObject, eventdata, handles)
    % hObject    handle to pushbutton1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    axes(handles.main_plot);
    cla;
    
    % Load in relaser data
    S = load('final_relaser.mat');
    
    popup_sel_index = get(handles.plot_popup, 'Value');
    switch popup_sel_index
        
        case 1 % Pump energy vs Laser energy
            hold on;
            xlim([0 ceil(S.penergy(S.IMAX))]);
            ylim([0 ceil(S.EABS_plot(S.IMAX))]);
            xlabel('Pump Energy');
            ylabel('Laser Energy');
            plot(S.penergy, S.EABS_plot, 'r-');
            legend('off')
            hold off;
            
        case 2 % t vs n5
            hold on;
            xlim([-inf inf]);
            ylim([-inf inf]);
            xlabel('Time (microsec)');
            ylabel('Manifold n5 population density');
            plot(S.T_print, S.printing_pops(1,:), 'b-');
            legend('off')
            hold off;
            
        case 3 % t vs n6
            hold on;
            xlim([-inf inf]);
            ylim([-inf inf]);
            xlabel('Time (microsec)');
            ylabel('n6');
            plot(S.T_print, S.printing_pops(2,:), 'b-');
            legend('off')
            hold off;
            
        case 4 % t vs n7
            hold on;
            xlim([-inf inf]);
            ylim([-inf inf]);
            xlabel('Time (microsec)');
            ylabel('n7');
            plot(S.T_print, S.printing_pops(3,:), 'b-');
            legend('off')
            hold off;
            
        case 5 % t vs n8
            hold on;
            xlim([-inf inf]);
            ylim([-inf inf]);
            xlabel('Time (microsec)');
            ylabel('n8');
            plot(S.T_print, S.printing_pops(4,:), 'b-');
            legend('off')
            hold off;
            
        case 6 % t vs phi
            hold on;
            xlim([-inf inf]);
            ylim([-inf inf]);
            xlabel('Time (microsec)');
            ylabel('FLUX');
            plot(S.T_print, S.printing_pops(5,:), 'b-');
            legend('off')
            hold off;
            
        case 7 % t vs Fracup*nt, Fraclo*n8
            hold on;
            xlim([-inf inf]);
            ylim([-inf inf]);
            xlabel('Time');
            ylabel('Population');
            plot(S.T_print, S.FracUp*S.printing_pops(3,:), 'b-');
            plot(S.T_print, S.FracLo*S.printing_pops(4,:), 'r-');
            legend({'FracUp*n7', 'FracLo*n8'}, 'Location', 'Southeast');
            hold off;
            
        case 8 % t vs Fracup*nt - Fraclo*n8
            hold on;
            xlim([-inf inf]);
            ylim([-inf inf]);
            xlabel('Time');
            ylabel('Population');
            plot(S.T_print, S.FracUp*S.printing_pops(3,:)-S.FracLo*S.printing_pops(4,:), 'r-');
            legend('off')
            hold off;
    end
end
