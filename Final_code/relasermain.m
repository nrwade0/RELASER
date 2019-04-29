%%  RELASER MAIN PROGRAM
%
% This program calculates the performance of a rare Earth laser from a
% series of rate equations. Approximately 50 input parameters are required
% for execution. These variables are given from two excel spreadsheets,
% laspar.xlsx (laser parameters) and specpar.xlsx (spectroscopic
% parameters). Details for each variables used in this program can be found
% in the variable list. The results are written to two files: popnum.txt
% containing manifold populations versus time and perfnum.txt containing
% laser energy versus pump energy. 
% 
% Version 1.0 (3.0 Fortran) - February 19, 2019 NW
%             - Transcribed line for line from Fortran program. Did not
%             worry about performance or ugliness, will clean later.
%             - Took all day and I mostly skipped using I/O files, just
%             read what Fortran was reading and declared the variables as
%             they were. Looking like close results so far.
%
% Version 1.1 - February 20-March 1, 2019 NW
%             - Removed a few redundant blocks and started to consolidate
%             into functions - this may just have to be avoided. Calling a 
%             .mat data file for ten functions isn't very efficient.
%             Perhaps another way or simply one long script?
%             - Code came pretty close in preliminary runs. Still off by a 
%             smidge, will need to be ironed out asap.
%             - Still no I/O files involved - was thinking something less
%             pedantic than counting columns in a text file would suffice.
%             Would also make generalizing easier down the road.
%
% Version 1.2 - March 2-6, 2019 NW
%             - Major effort in input being completed. Decided to use excel
%             files for pulling data since it is easy to read/change and
%             significantly easier to generalize. Working well so far
%             although it is a bit slower. I think the trade-off is
%             ultimately worth it.
%             - Started placing comments where I can - getting more
%             familiar with the variable terminology, although some still
%             elude me (THRESH?).
%
% Version 1.3 - March 7-13, 2019 NW
%             - Another major effort in files, but this time in output.
%             popnum and perfnum are outputted to text files. popnum is
%             probably best kept as a text file with N lines.
%             However, perfnum could be given as a spreadsheet in excel 
%             (The power of modern programming!). Or possibly .csv if you
%             need to plot it using MATLAB or python.
%             - Still off about 10 microseconds in pulse width. I'm not sure
%             what I'm missing for it and it's eating me up. As soon as
%             that is fixed, I'll have more time for a GUI skeleton.
%
% Version 1.4 - March 14-25, 2019 NW
%             - Started making moves towards implementing a GUI. See GUI
%             version history for more details.
%             - I still can't find the pwidth error. It may be some
%             difference in precision between Fortran and MATLAB. I'd like
%             to go through line-by-line and comment anyways. Two birds,
%             one stone?
%             - Added the dumpin and dumpin2 functions to output a datadump
%             text file. One file per energy. Same formatting as the 
%             original program.
%
% Version 1.5 - March 26-April 5, 2019 NW
%             - Stopped myself short of changing all the variable names -
%             I'm starting to get used to them. Updated variable list is
%             coming along, but I am having trouble describing them in
%             detail. The science is lacking.
%             - Input excel files now separated into sheets for easier use.
%             This may slow the program, but I would consider it a worthy
%             trade-off.
%             - Still working to generalize input/output files, but haven't
%             ventured into a new testing model until the current one is
%             accurate.
%             - Some updates on GUI.
%             - Considering the END GOAL is to replicate "Spectroscopy and
%             modeling of solid state lanthanide lasers: Application..." by
%             B. Walsh, April 2004, specifically his Fig. 8,9.
%             - I stopped naming the file with version numbers since it is
%             connected to GUI and I don't want to bother tracking down the
%             necessary changes.
%             
% Version 1.6 - April 8-15, 2019 NW  ----------- CURRENT
%             - GUI has, essentially, been completed.
%             - Commenting and relaser variable list was completed.
%             - Basic debugging for Pwidth error continued. Found an error
%             in resetting the initial value of y(9) (photon denisty) to 0.
%             It sets it to 5...? Not sure where it that happens. Fixed 
%             this error at the end of import_specpar. See below for more 
%             information on the Pwidth error.
%
% Future updates/issues
%    - First and foremost, RELASER doesn't give the exact same value as the FORTRAN
%     version. Will need to walkthrough and confirm no small detail is
%     omitted.
%    - Final task is to replicate '04 paper called "Spectroscopy and modeling
%     of solid state lanthanide lasers: Application to trivalent Tm3+ and
%     Ho3+ in YLiF4 and LuLiF4." Will need to scour for parameters from this
%     paper. Generalize this implementation with any number of parameters.
%    - Solve using ode15s instead of piecewise integration. May get rid of
%     Pwidth issue, but requires a solid understanding of the program.
%    - Add extensive commenting (import from FORTRAN, as well)
%
% Pwidth error
%    What I've tried so far: Debugging line-by-line in the CRUNCH() and 
%    RATEEQ() functions of the Fortran program against relasermain.
%    Unfortunately, I have yet to find the error. I know it is something
%    affecting the ydot array of the RATEEQ function but I can't find where
%    exactly. It seems to be always a ten thousandth of a decimal place
%    off, which makes me think it is a precision error but MATLAB default =
%    REAL(8) in Fortran (or so the internet tells me). The dumpit files do
%    not go into enough detail for REACT() and ydot() which make it
%    impossible to track the changes.
% 
%% --------------------------- PROGRAM PREAMBLE ----------------------------% 
clc
clear variables
close all
format long

%% --------------------------- TALK TO THE USER ----------------------------%
% Print program name and retrieve the run title.
fprintf("Program RELASER \n")

% Gives default title as RELASER
title = "RELASER";
disp("Default title is 'RELASER'");

% Asks for a title
%title = input('Enter the title of this run: ', 's');

%% ------------------------- EXPLICITLY DEFINE CONSTANTS ---------------------%
% Number of times program has run (use for output?)
icase = 0;

% Constants: c = speed of light and hc = Planck's const*c
c = 2.99792e4;    % [cm/us]
%c = 3e4;% Interestingly, this value of c gives a closer pwidth than currently
hc = 1.98648e-19; % [J um]
% PI supplied in MATLAB as 'pi'

%% --------------------- IMPORT DATA FROM LASPAR.XLSX -----------------------%
% Import laspar.xlsx as a matrix (instead of table). NOTE: laspar matrix 
% does not include row 1 of laspar.xlsx, since it only contains NaNs.
laspar = xlsread('laspar.xlsx', 'data');
%disp(laspar);

% LASER PARAMETERS
% Take cells from laspar matrix and store into proper variables
cavl     = laspar(1,3);   % [cm] Cavity length
cmrad    = laspar(2,3);   % [cm] Cavty mode radius
rodl     = laspar(3,3);   % [cm] Rod length
rodrad   = laspar(4,3);   % [cm] Rod radius
tq       = laspar(5,3);   % Q-switch delay
tqi      = laspar(6,3);   % Q-switch open interval
restype  = laspar(7,3);   % Ring (1) or Linear (2) resonator
pumptype = laspar(8,3);   % End (1) or Side (2) pumping

% INTEGRATION PARAMETERS
deltat   = laspar(9,3);   % [us] Stepsize
Npoints  = laspar(10,3);  % Number of points
tout     = laspar(11,3);  % [us] Time interval between points

% Reflectance arrays are filled from laspar matrix.
% See notes sheet to see how to add more values.
% Main routine cycles thru each J1MAX, JMAX, & IMAX.
J1MAX    = laspar(13,3);  % Number of REAR reflectors
backrf   = zeros(J1MAX,1);% Rear reflector reflectances array
for i = 1:J1MAX
    backrf(i) = laspar(13,3+i);
end

JMAX     = laspar(15,3);  % Number of FRONT reflectors
frntrf   = zeros(JMAX,1); % Front reflector reflectances array
for i = 1:JMAX
    frntrf(i) = laspar(15,3+i);
end

IMAX      = laspar(17,3);  % Number of pump energies
penergy   = zeros(IMAX,1); % Pump energy array
for i = 1:IMAX
    penergy(i) = laspar(17+i,3);
end

% Calculate cavity area and rod area
cmarea  = pi*cmrad*cmrad;
rodarea = pi*rodrad*rodrad;

% Set pump absorption length depending on pump type
switch pumptype
    case 1 % End pumping
        pabsl = rodl;
        
    case 2 % Side pumping
        pabsl = 2*rodrad;
end

% IQSWITCH IS USED FOR WHAT?????
% Working on implementing Q-switch delay (some commented out in FORTRAN)
IQSWITCH = 1; % Q-switch flag
if(tq <= 0) % negative and zero tq is turned to 0?
    IQSWITCH = 0;
    tq = 0;
end


%% --------------------------- CREATE OUTPUT FILES -------------------------%
% Create popnum.txt and perfnum.txt files to record output
%     popnum.txt contains manifold populations versus time
%     perfnum.txt contains laser energy versus pump energy
% 'wt' signifies writing to and replaces previous file of same name
popnum = fopen('popnum.txt','wt');
perfnum = fopen('perfnum.txt','wt');

% Write title on first line of each
fprintf(popnum,'RELASER Run Title: %s \n',title);
fprintf(perfnum,'RELASER Run Title: %s \n',title);

% Write date on second line of each
time = datestr(clock,'mm/dd/YYYY HH:MM');
fprintf(popnum, 'Date: %2s \n', time);
fprintf(perfnum, 'Date: %2s \n', time);


%% ------------------------- IMPORT DATA FROM SPECPAR -----------------------%
% Import each sheet from specpar.xlsx and store in proper variables.
% Output is filinpt_dump.mat data file holding all variables.
% This data file is loaded on each run of pump energy (third for loop).
% Previously FILINPT subroutine
import_specpar();

% used for printing popnum
printing_pops = zeros(IMAX, Npoints);

%
%% ---------------------------- BEGIN CALCULATION --------------------------%
% Displays onto command line
disp('           Begin Timer          ');tic; % Begin timer
disp(' --------- CALCULATING ---------');


for J1 = 1:J1MAX % for each REAR REFLECTOR
    
    for J = 1:JMAX % for each FRONT REFLECTOR
        
        % Print Rear/Front reflector number above the given values
        fprintf('Front Reflector %d |  Rear Reflector %d \n', J, J1)
        fprintf('OC = %0.2f         |  HR = %0.2f  \n', frntrf(J), backrf(J1));
        
        % Print table header to command line
        fprintf("%8s %8s %8s %8s %8s %8s %10s %10s %10s %10s \n",...
                'TEMP','PUMPE','EABS','TOUT','LOSS','DELAY','WIDTH','EOUTLAM1','EOUTLAM2','EOUTLAM3');
        
        % Arrays used for plots
        EABS_plot = zeros(IMAX,1);
        
        for e = 1:IMAX % for each PUMP ENERGY
            %% -------------------- LOAD DATA and SETUP ---------------------%
            % Store No. of runs in icase
            % Is this necessary for output?
            icase = icase + 1;
            
            % Load specpar parameters.
            load('specpar_dump.mat');
            
            % Set back and front reflectances to current values
            reflb = backrf(J1);
            refl = frntrf(J);
            
            % Fill column arrays with the current reflectances.
            [Reflout, Reflbck] = deal(zeros(Nlaslam,1));
            for K = 1:Nlaslam
                Reflout(K) = refl;
                Reflbck(K) = reflb;
            end
            
            % Current pump energy stored
            pumpe = penergy(e);
            
            % Calculate wavelength dependent terms & pump scale terms
            [optlng, spo, roundtt] = deal(zeros(Nlaslam,1));
            for i = 1:Nlaslam
                optlng(i) = rodl*roddex(i) + cavl - rodl; % Optical cav len
                switch pumptype
                    case 1
                        B = cmarea/(4*pi*optlng(i)^2);
                    case 2
                        B = rodarea/(4*pi*optlng(i)^2);
                end
                spo(i) = B*rodl/optlng(i); % Spontaneous emission factor
                roundtt(i) = restype*optlng(i)/c; % Round trip time
            end

            % No. of rate equations = no. of manifolds + no. of wavelengths
            mancount = sum(Nman);
            NEQ = mancount + Nlaslam;
            
            % Define pump scale parameters depending on pump type
            if(ipump == 0)
                fwhm = Npoints*tout;
            end
            
            % Define pump scale parameter, W1
            switch pumptype
                case 1
                    W1 = pumpeta*pumpe*(pumplam/hc)/(cmarea*rodl);
                case 2
                    W1 = pumpeta*pumpe*(pumplam/hc)/(rodarea*rodl);
            end
            
            
            %% ----------------- INITIALIZE OUTPUT FILES ---------------------%
            % Write column headers on third line of popnum and perfnum txts
            % Prints headers only on first go thru (e == 1 only)
            if(e == 1)
                % perfnum file given same headers as cmd line disp();
                % '%#s' signifies # spaces allotted for each column
                fprintf(perfnum, "%8s %8s %8s %8s %8s %8s %10s %10s %10s %10s \n",...
                'TEMP','PUMPE','EABS','TOUT','LOSS','DELAY','WIDTH','EOUTLAM1','EOUTLAM2','EOUTLAM3');
                
                % popnum file gives 'TIME' column first then cycles thru 
                % and prints the names of upper manifolds in array 'MAN'.
                % Finally, it prints 'FLUX' for last column.
                fprintf(popnum, "%13s", 'TIME');
                for i = 1:Ncol-1
                    fprintf(popnum, "%14s", man(i));
                end
                fprintf(popnum,"%14s \n", 'FLUX');
            end
            
            %% ---------------------- MAIN CALCULATION ----------------------%
            % THE MAIN CALCULATION - CRUNCH()
            % Save workspace for crunch function and load data after it is
            % completed.
            save('crunch_data.mat')
            crunch();
            load('crunchdatadump.mat');
            EABS_plot(e) = EABS;
            
        end
    end
end

% Save final matlab data, goes to GUI
save('final_relaser.mat');

% vars = load('final_relaser.mat');

% Close files
fclose(popnum);
fclose(perfnum);

% Print elapsed time
toc;

% open GUI
relaser_gui_v2


% END



%% ------------------------------- import_specpar() -------------------------------%
% - No formal input parameters
% - Output is saved in a .mat data file and includes imports from specpar.
% - 6 sheets in specpar.xlsx give information regarding the processes
% occuring in the laser simulation. Data is taken in the most general way
% possible to allow for changes to be made on the fly. See notes sheet for
% more detail on how to add new processes. Sheets are named for their
% respective common blocks from original FORTRAN code. Indexing for
% importing is explained in notes of specpar.xlsx. 
function import_specpar()
    %% ------------------------------ /CONCK/ --------------------------------%
    %             GET NUMBER OF IONS, IONS, AND CONCENTRATIONS
    
    % Import conck data into three matrices containing numbers and strings.
    % Raw data matrix is not used and omitted with ~
    [num_conck, txt_conck, ~] = xlsread('specpar.xlsx','conck');
    Nions = num_conck(2,2);
    %disp(raw_conk);
    
    % for each ion given...
    % K1 signifies a dummy ion (i.e. 'TM' when K1 = 1)
    [ion, man] = deal(strings(Nions, 1));
    Nman = zeros(Nions, 1);
    for K1 = 1:Nions
        
        % Save ion name in ION array.
        % Store no. of manifolds array NMAN.
        % Current No. of manifolds stored in MAXMAN.
        ion(K1)  = txt_conck(3+(K1-1)*6,2);
        Nman(K1) = num_conck(4+(K1-1)*6,2);
        maxman   = Nman(K1);
        mancount = sum(Nman);

        % preallocate cell array MAN, array of concentrations Y, and
        % array of total concentrations for each ion TOTLCON
        Y = zeros();
        totlcon = zeros(Nions,1);
        
        % Gather information on each manifold
        for i = 1:maxman
            % K stores current manifold number (B5:B8).
            % MAN stores the string for that manifold number (C5:C8).
            % Y stores initial concentrations in the order of manifold
            %     numbers (E5:E8).
            % TOTLCON sums up those concentrations into an array that holds
            %     total concentrations for each ion (SUM(E5:E8)).
            K           = num_conck(4+i+(K1-1)*6, 2);
            man(K)      = txt_conck(4+i+(K1-1)*6, 3);
            Y(K)        = num_conck(4+i+(K1-1)*6, 5);
            totlcon(K1) = totlcon(K1) + Y(K);
        end
    end

    %% ------------------------------ /DECAY1/ -------------------------------%
    %                      GET SINGLE ION DECAY TIMES
    
    % Import DECAY1 data from specpar xlsx. Only number data is used.
    [num_decay1, ~, ~] = xlsread('specpar.xlsx','decay1');
    
    % Number of decay times (B2)
    Ndktime = num_decay1(2,2);
    
    % Collect data for each decay time
    [manup,manlo,dkrate] = deal(zeros(Ndktime,1));
    for i = 1:Ndktime
        % MANUP is upper manifold array
        % MANLO is lower manifold array
        % DKTIME is decay time for jump from upper to lower manifold
        %    DKTIME is converted to DKRATE
        K1        = num_decay1(2+i, 2);
        manup(K1) = num_decay1(2+i, 4);
        manlo(K1) = num_decay1(2+i, 6);
        dktime    = num_decay1(2+i, 8);
        % 
        if(dktime > 0)
            dkrate(K1) = 1/dktime;
        else
            dkrate(K1) = 0;
        end
    end

    %% ------------------------------ /DECAY2/ -------------------------------%
    %                      GET BILIEAR TRANSFER PARAMETERS
    
    % Import data from specpar. Only number data is used
    [num_decay2, ~, ~] = xlsread('specpar.xlsx','decay2');
    
    % Get the number of bilinear transfers (B2)
    Nbilin = num_decay2(2,2);
    
    % for each bilinear transfer process... collect data
    % K1 represents the index, could be replaced by iterator i
    [ManLos1,ManLos2,ManGan1,ManGan2, ...
        Naceptr,bilinrt] = deal(zeros(Nbilin,1));
    for i = 1:Nbilin
        %
        K1          = num_decay2(2+i, 2);
        ManLos1(K1) = num_decay2(2+i, 4);
        ManLos2(K1) = num_decay2(2+i, 6);
        ManGan1(K1) = num_decay2(2+i, 8);
        ManGan2(K1) = num_decay2(2+i, 10);
        bilintm     = num_decay2(2+i, 12);
        Naceptr(K1) = num_decay2(2+i, 14);
        %
        %K = Naceptr(K1);
        % bilinrt cannot be less than 0
        if(bilintm <= 0)
            bilinrt(K1) = 0;
        else
            bilinrt(K1) = bilintm;
        end
    end
    
    %% ----------------------------- /WAVLNG/ --------------------------------%
    %                  GET LASER WAVELENGTHS & TEMPERATURE
    
    % Read wavlng data. Only number matrix is used
    [num_wavlng, ~, ~] = xlsread('specpar.xlsx','wavlng');
    
    % Get number of laser wavelengths and the specified temperature
    % TEMPER is only printed as of now.
    % Will need to change around to allow for differing temperatures
    Nlaslam = num_wavlng(2,4);
    temper = num_wavlng(2,10);

    % Fill arrays for each laser wavelength
    [ManUpLz,ManLoLz,emislam,FracUp,FracLo,stim,spon,cavalfa,...
        Reflout,Reflbck,surfacet,roddex] = deal(zeros(Nlaslam, 1));
    for i = 1:Nlaslam
        %
        % K1 is the current wavelength number
        % Can be replaced with i
        K1 = num_wavlng(3+(i-1)*10, 2);
        %
        % See RELASER list for details on each variable
        ManUpLz(K1)  = num_wavlng(3+(i-1)*10, 4);
        ManLoLz(K1)  = num_wavlng(3+(i-1)*10, 6);
        emislam(K1)  = num_wavlng(3+(i-1)*10, 10);
        FracUp(K1)   = num_wavlng(4+(i-1)*10, 10);
        FracLo(K1)   = num_wavlng(5+(i-1)*10, 10);
        stim(K1)     = num_wavlng(6+(i-1)*10, 10);
        sponlife     = num_wavlng(7+(i-1)*10, 10);
        if(sponlife > 0)
            spon(K1) = 1/sponlife;
        else
            spon(K1) = 0;
        end
        cavalfa(K1)  = num_wavlng(8+(i-1)*10, 10);
        Reflout(K1)  = num_wavlng(9+(i-1)*10, 10);
        Reflbck(K1)  = num_wavlng(10+(i-1)*10, 10);
        surfacet(K1) = num_wavlng(11+(i-1)*10, 10);
        roddex(K1)   = num_wavlng(12+(i-1)*10, 10);
    end

    %% ----------------------------- /PUMPK/ ---------------------------------%
    %                       GET PUMPING PARAMETERS
    
    % Number data is only used, num_pumpk is matrix holding data
    [num_pumpk, ~, ~] = xlsread('specpar.xlsx','pumpk');
    
    % Get pump wavelength and no. of light absorption manifolds
    pumplam = num_pumpk(2,2);
    Npmplv = num_pumpk(3,2);

    % Get values for each light absorbing manifold
    manpmp = zeros(Npmplv, 2);
    pumpsig = zeros(Npmplv, 1);
    for i = 1:Npmplv
        % indexing: MANPMP cell starts in row 4 but every time
        % there is another wavelength (ie. NPMPLV+1) row 4 and 5
        % will presumably be push down. Therefore, shoving it down
        % by 2 for every iteration of this for loop. This is
        % assuming the bottom rows ("Has Energy... Pump alignment
        % efficiency" are scooted down as well.
        manpmp(i,1) = num_pumpk(4+(i-1)*2, 4);
        manpmp(i,2) = num_pumpk(4+(i-1)*2, 6);
        pumpsig(i)  = num_pumpk(5+(i-1)*2, 7);
    end

    % indexing depends on how many wavelengths are given... 
    % if more than one is given, these values will be shoved down
    % by 2 for each beyond the typical one. Hopefully that suffices
    % Why pumpe not used?
    pumpe   = num_pumpk(6+(Npmplv-1)*2, 7);
    fwhm    = num_pumpk(7+(Npmplv-1)*2, 7);
    ipump   = num_pumpk(8+(Npmplv-1)*2, 7);
    pumpeta = num_pumpk(9+(Npmplv-1)*2, 7);
    pmpalgn = num_pumpk(10+(Npmplv-1)*2, 7);
    
    if pumpe <= 0
        fwhm = 0;
    end

    %% ------------------------------ /COLS/ ---------------------------------%
    %                 GATHER CONCENTRATION COLUMN INFORMATION
    
    % Matrix num_cols has number data from specpar
    [num_cols, ~, ~] = xlsread('specpar.xlsx','cols');
    
    % Number of concentration columns
    Ncol = num_cols(2,2);

    % Gather column output information
    [icol,jcol,isum] = deal(zeros(Ncol,1));
    for i = 1:Ncol
        %
        icol(i) = num_cols(2+i,2);
        jcol(i) = num_cols(2+i,4);
        isum(i) = num_cols(2+i,6);
    end

    % Add photon initial values to Y array.
    % Found below column data in specpar.
    for i = (mancount+1):(mancount+Nlaslam+1)
        % Indexing refers to adding after the last Y value (in Y(mancount))
        % and adding one for each laser wavelength. (i.e. B8 & B9)
        Y(i) = num_cols(2+Ncol+i);
    end
    
    Y(9) = 0;
    
    %% --------------------------- SAVE IMPORTED DATA -----------------------%
    % Save data from import_specpar function.
    save('specpar_dump.mat');
    
end



%% ------------------------------- CRUNCH() --------------------------------%
% - No formal input parameters
% - Output is workspace variables saved in .mat data file
% - The three final functions are called in crunch with seperate .mat data
% files for input/output themselves.
% - The rate equations are solved simultaneously using a piecewise
% integration technique?
function crunch()
    
    % Load workspace data
    % Is there a more efficient way? Structure?
    load('crunch_data.mat');
    
    %% ----------------- INITIALIZE VARIABLES AND DUMP PARAMETERS -----------%
    % What does iflag do?
    % plsmin/plsmax redeclared
    iflag = 0;
    plsmin = 0;
    plsmax = 0;
    
    % Nreact is an array for the number of reactions occuring.
    % filled with zeros. T = T change.
    Nreact = Ndktime + Nbilin + Npmplv + 8*Nlaslam;
    react = zeros(Nreact,1);
    T = deltat;
    
    % Threshold density?
    % Unsure as to where these equations originated?
    thresh = zeros(Nlaslam);   
    thrsmn = 1e30;
    for i = 1:Nlaslam
        T1 = 2 - Reflout(i) - Reflbck(i) - 2*surfacet(i) + 2*rodl*cavalfa(i);
        GAMMA = 1 + (FracLo(i)/FracUp(i));
        T2 = restype*stim(i)*rodl*(GAMMA - 1);
        T3 = restype*stim(i)*rodl*GAMMA;
%         T2 = 2.0*stim(i)*rodl*FracLo(i);
%         T3 = 2.0*stim(i)*rodl*(FracLo(i) + FracUp(i));
        thresh(i) = totlcon(1);
%         if(T3 > 0)
%             thresh(i) = (T1 + T2*totlcon(1))/T3;
%         end
        thresh(i) = (T1 + T2*totlcon(1))/T3;
%         if(thresh(i) < thrsmn)
            thrsmn = thresh(i);
            T1S = T1;
            T2S = T2;
            T3S = T3;
%         end
        
    end
    
    % Passing thru a struct to dumpin()instead of .mat saves on runtime.
    S = struct('e',e, 'time',time, 'deltat',deltat, 'pumpe',pumpe, ...
        'tout',tout, 'Npoints',Npoints, 'cavl',cavl, 'cmarea',cmarea, ...
        'rodl',rodl, 'rodrad',rodrad, 'pabsl',pabsl, 'tq',tq, 'tqi',tqi, ...
        'ipump',ipump, 'fwhm',fwhm, 'Npmplv',Npmplv, 'pumpsig',pumpsig, ...
        'pumplam',pumplam, 'W1', W1, 'pumpeta',pumpeta, 'pmpalgn',pmpalgn, ...
        'manpmp',manpmp, 'Nions',Nions, 'ion',ion, 'Nman',Nman, ...
        'mancount',mancount, 'man',man, 'Y',Y , 'totlcon',totlcon, ...
        'Ndktime',Ndktime, 'manup',manup, 'manlo',manlo, 'dkrate',dkrate, ...
        'Nbilin',Nbilin, 'ManLos1',ManLos1, 'ManGan1',ManGan1, ...
        'ManLos2',ManLos2, 'ManGan2',ManGan2, 'bilinrt',bilinrt, ...
        'Naceptr',Naceptr, 'temper',temper, 'Nlaslam',Nlaslam, ...
        'ManUpLz',ManUpLz, 'ManLoLz',ManLoLz, 'emislam',emislam, ...
        'FracUp',FracUp, 'FracLo',FracLo, 'stim',stim, 'spon',spon, ...
        'cavalfa',cavalfa, 'Reflout',Reflout, 'Reflbck',Reflbck, ...
        'surfacet',surfacet, 'roddex',roddex, 'optlng',optlng, 'spo',spo, ...
        'roundtt',roundtt, 'thresh',thresh);
    dumpin(S);
    
    % Converts reflectance, transmittace, and bulk absorption to
    % equivalent loss factors. Keeping them seperate rather than a
    % single cavity lifetime factor allows individual bookkeeping
    % of their impact on laser performance.
    outsav = zeros(Nlaslam);
    for i = 1:Nlaslam
        outsav(i) = -log(Reflout(i));
        cavalfa(i) = restype*rodl*cavalfa(i)/roundtt(i);
        Reflout(i) = -log(Reflout(i))/roundtt(i);
        Reflbck(i) = -log(Reflbck(i))/roundtt(i);
        surfacet(i) = 2*(1 - surfacet(i))/roundtt(i);
    end
    
    % Number of integrations each slice
    tspan = floor(tout/deltat+0.001);
    
    T_print = zeros(1, Npoints);
        
    %% ---------------- FOR EACH INTEGRATION SLICE (1 TO Npoints) ---------------%
    % N goes for each Npoints (10,000)
    for N = 1:Npoints
        
        % Set up manifold densities
        % For each rate equation, set Ymax, why?
        Ymax = zeros(NEQ,1);
        for K = 1:NEQ 
            Ymax(K) = -1e20;
        end
        
        %% ---------------- FOR EACH INTEGRATION SLICE (1 TO Npoints) ---------------%
        % t represents a small step in (tout/deltat) inside each jump of N.
        % i.e. do this tspan times each integration slice.
        for t = 1:tspan
            
            % Current reaction term
            nterm = 0;
               
            % preallocate change in densities array
            ydot = zeros(NEQ,1);
            
            %  Begin integration of rate equations
            %% ------------------- LINEAR DECAY PROCESSES -------------------%
            for i = 1:Ndktime
                % Upper/lower manifold densities
                K1 = manup(i);
                K2 = manlo(i);
                
                % R1 is a decay term from upper manifolds
                R1 = Y(K1)*dkrate(i);
                % Reaction @ nterm = [1-4] is +R1
                nterm = nterm + 1;
                react(nterm) = react(nterm) + R1;
                % Upper/lower manifolds affected by decay
                % Upper decreases, lower increases
                ydot(K1) = ydot(K1) - R1;
                ydot(K2) = ydot(K2) + R1;
            end
            
            
            %% ------------ BILINEAR ENERGY TRANSFER PROCESSES ------------%
            for i = 1:Nbilin
                % Upper/lower manifolds loss/gain
                K1 = ManLos1(i);
                K2 = ManLos2(i);
                K3 = ManGan1(i);
                K4 = ManGan2(i);
                % Each transfer given variable R1
                R1 = Y(K1)*Y(K2)*bilinrt(i);
                % Reaction @ nterm = 5,6 is +R1
                nterm = nterm + 1;
                react(nterm) = react(nterm) + R1;
                % Energy transfers are applied to each manifold.
                % First two are losses, second two are gains.
                ydot(K1) = ydot(K1) - R1;
                ydot(K2) = ydot(K2) - R1;
                ydot(K3) = ydot(K3) + R1;
                ydot(K4) = ydot(K4) + R1;
            end
            
            
            %% ------------------- THE PUMPING TERMS -----------------------%
            % Temp variables ssum, ssum2, arg, and RT are constituents of
            % the larger equation for the pumping terms RR, R1.
            % Perhaps this can be condensed to avoid wasting space on those
            % variables?
            ssum = 0;
            ssum2 = 0;
            for i = 1:Npmplv
                ssum = ssum - pumpsig(i)*pabsl*Y(manpmp(i,1));
                ssum2 = ssum2 + 1 - exp(-pumpsig(i)*pabsl*Y(manpmp(i,1)));
            end
            
            RT = W1*pumpval(ipump,fwhm,T)*(1 - exp(ssum));
            
            % Compute more temp variables for each pump.
            for i = 1:Npmplv
                arg = 1 - exp(-pumpsig(i)*pabsl*Y(manpmp(i,1)));
                R1 = RT*arg/ssum2;
                RR = R1;
                % Reaction @ nterm = 7 is +R1
                nterm = nterm +1;
                react(nterm) = react(nterm) + R1;
                % Pumping terms are applied to the manifold densities.
                ydot(manpmp(i,1)) = ydot(manpmp(i,1)) - RR;
                ydot(manpmp(i,2)) = ydot(manpmp(i,2)) + RR;
            end
            
            
            %% --- CAVITY LOSS, SPONTANEOUS AND STIMULATED EMISSION TERMS ---%
            % Compute for each wavelength into laser (from WAVLNG)
            for i = 1:Nlaslam
                % Laser manifolds and photons dummy variables
                %   K1 = lower laser manifold (7), K2 = upper laser
                %   manifold (8), K3 = photon density (9).
                K1 = ManLoLz(i);
                K2 = ManUpLz(i);
                K3 = mancount + i;
                
                % Effects by stimulated emission on upper/lower manifolds
                R1UP = stim(i)*c*Y(K3)*FracUp(i)*Y(K2);
                R1LO = stim(i)*c*Y(K3)*FracLo(i)*Y(K1);
                
                % Other parameter effects
                R1 = R1UP - R1LO;
                R2 = R1*rodl/optlng(i);
                R3 = spo(i)*FracUp(i)*Y(K2)*spon(i);
                % Parameter effects on photon density
                R4 = Y(K3)*cavalfa(i);
                R5 = Y(K3)*surfacet(i);
                R6 = Y(K3)*Reflbck(i);
                R7 = Y(K3)*Reflout(i);
                
                % Cavity gains by stimulated emission
                % Cavity losses by parameters
                gain = stim(i)*c*(FracUp(i)*Y(K2) - FracLo(i)*Y(K1));
                xloss = cavalfa(i) + surfacet(i) + Reflbck(i) + Reflout(i);
                
                %% --------------------- Q-switch routine -------------------%
                % Determine Q-switch term (R8?)
                R8prime = 0;
                %R8 = 0;
                
                % tq = Q-switch delay, T = Change in time, tqi = Q-switch 
                % open inteval
                
                % "If there is a Q-switch delay and it is less than the
                % current change in time..."
                if(tq > 0 && T < tq)
                    R8prime = 1e8;
                end
                % "If the Q-switch open interval exists and the sum of the
                % interval and the delay is greater than the current change
                % in time..."
                if(tqi > 0 && T > tq+tqi)
                    R8prime = 1e8;
                end
                
                % Q-switch term
                R8 = Y(K3)*R8prime;
                
                % outsav is -ln(reflectance) so what does this mean?
                % Is this some kind of extreme boost/loss from spontaneous
                % emission? Does not occur in current run: -ln(0.4) = 0.9
                if(outsav(i) < 0.01)
                    xloss = xloss + 1e8;
                    Ydum = Y(K3) + R3*deltat;
                    R4 = Ydum*cavalfa(i);
                    R5 = Ydum*surfacet(i);
                    R6 = Ydum*Reflbck(i);
                    R7 = Ydum*Reflout(i);
                    R8 = Ydum*1e8;
                end
                
                % Cavity gain/loss is recalculated
                gain = gain*(Y(K3) + R3*deltat)*rodl/optlng(i);
                xloss = (xloss + R8prime)*(Y(K3) + R3*deltat);
                
                % Total of cavity gain/loss terms
                RT = R4 + R5 + R6 + R7 + R8;
                
                % What is this checking for?
                % If the photon density is decreasing over deltat, do...
                if((Y(K3) + (gain - xloss)*deltat) < 0)
                    % Stimulated emission is added to current photon
                    % density into 'hm'
                    hm = Y(K3) + (R2 + R3)*deltat;
                    % denom is given total cavity losses over deltat
                    denom = RT*deltat;
                    recip = 0;
                    if(denom ~= 0)
                        recip = 1/denom;
                    end
                    
                    % Loss terms are minimized (?) if photon density is
                    % decreasing.
                    R4 = R4*hm*recip;
                    R5 = R5*hm*recip;
                    R6 = R6*hm*recip;
                    R7 = R7*hm*recip;
                    R8 = R8*hm*recip;
                    RT = ydot(K3) + R3 + R2 + Y(K3)/deltat;
                end

                % R terms are applied to each manifold population
                RR = RT;
                ydot(K1) = ydot(K1) + R1;
                ydot(K2) = ydot(K2) - R1;
                ydot(K3) = ydot(K3) + R2 + R3 - RR;
                
                % Reaction @ nterm = 7-15 is +R1, +R2, ..., +R8, respect.
                react(nterm + 1) = react(nterm + 1) + R1;
                react(nterm + 2) = react(nterm + 2) + R2;
                react(nterm + 3) = react(nterm + 3) + R3;
                react(nterm + 4) = react(nterm + 4) + R4;
                react(nterm + 5) = react(nterm + 5) + R5;
                react(nterm + 6) = react(nterm + 6) + R6;
                react(nterm + 7) = react(nterm + 7) + R7;
                react(nterm + 8) = react(nterm + 8) + R8;
                
                % Increment nterm after 9 uses
                % Unnecessary
                nterm = nterm + 9;
            end
            
            %% ------------- STEPWISE INTEGRATION OF TERMS -----------------%
            % For each rate equation...
            for i = 1:NEQ
                % Euler's method on steroids?
                Y(i) = Y(i) + ydot(i)*deltat;
                
                % Population cannot dip to a negative
                if(Y(i) < 0)
                    Y(i) = 0;
                end
                
                % Photon count less that 1e-10 then just set equal to zero.
                % Cannot see 1e-10 photons anyway?
                if(Y(i) < 1e-10 && i > mancount)
                    Y(i) = 0;
                end
                
                % Set new max population value.
                if(Y(i) > Ymax(i))
                    Ymax(i) = Y(i);
                end
                
                % Set new lowest possible max population value to 1e-32.
                if(Ymax(i) <= 1e-32)
                    Ymax(i) = 1e-32;
                end
            end
            
            %% ----------------- LASER WAVELENGTH FIX? ---------------------%
            % For each possible laser wavelength...
            for i = 1:Nlaslam
                
                % Upper laser manifold dummy value K1, photon index K2
                K1 = ManUpLz(i);
                K2 = mancount + i;
                
                % "If upper laser manifold pop. is less than [some kind of
                % threshold] OR photon density is less than 1e10
                if(Y(K1) < .95*(T1S + T2S*(totlcon(1) - Y(10)))/T3S...
                        || Y(K2) < 1e10)
                    continue;
                end
                
                if(Y(K1) < .95*(T1S + T2S*(totlcon(1)))/T3S)
                    continue;
                end
                
                plsmax = T;
                if(iflag == 1)
                    continue;
                end
                iflag = 1;
                plsmin = T;
            end
            T = T + deltat;
        end
        
        %% ----------------- NEW YMAX AND PRINT TO POPNUM ------------------%
        % Set new max values
        % Unnecessary
        for i = 1:mancount
            Ymax(i) = Y(i);
            
            % Saves Y densities for plotting
            %Y_plot(i,N) = Y(1,i);
        end
        
        % for writing purposes
        if(e == IMAX)
            Yprint = zeros(Ncol,1);
            fprintf(popnum, '%13.5e ', T); 
            T_print(N) = T;
            for i = 1:Ncol % for each required column
                ijcol = jcol(i);
                Yprint(i) = 0;
                for j = 1:ijcol % for each sum of how many manifolds (1 usually)
                    Yprint(i) = Yprint(i) + Ymax(isum(i,j));
                end
                fprintf(popnum, '%13.5e ', Yprint(i));
                printing_pops(i, N) = Yprint(i);
            end
            fprintf(popnum,"\n");
        end
    end

    % Passing thru a struct instead of .mat saves a second on runtime.
    S = struct('e',e, 'deltat',deltat, 'pumpe',pumpe, ...
        'Nreact',Nreact, 'Ndktime',Ndktime, 'Nbilin',Nbilin, 'Npmplv',Npmplv, 'Nlaslam',Nlaslam, ...
        'rodl',rodl, 'cmarea',cmarea, 'mancount',mancount, 'pumptype',pumptype, 'pmpalgn',pmpalgn, ...
        'rodarea',rodarea, 'ManLos1',ManLos1, 'ManLos2',ManLos2, 'ManGan1',ManGan1, 'ManGan2',ManGan2, ...
        'manpmp',manpmp, 'hc', hc, 'pumplam',pumplam, 'optlng',optlng, ...
        'emislam',emislam, 'plsmin',plsmin, 'Reflout',Reflout, 'Reflbck',Reflbck, ...
        'cavalfa',cavalfa, 'surfacet',surfacet, 'Y',Y , 'roundtt',roundtt, ...
        'temper',temper, 'react',react, 'manup',manup, 'manlo',manlo, 'NEQ',NEQ, 'man',man, 'plsmax',plsmax,...
        'perfnum',perfnum);
    
    dumpin2(S);
	load('dataoffdump2.mat');
    save('crunchdatadump.mat');
    
end



%% --------------------------- DUMPIN(struct S) -----------------------------%
% - S is a structure of required variables to dump
% - Output is new data dump txt file, no .mat data file here
% - Data is just dumped into file line by line in a format similar to that
% of the original Fortran file
function dumpin(S)
    %% ----------------------------- LOAD DATA ------------------------------%
    %
    %
    %% ---------------------- CREATE FILE AND HEADER ------------------------%
    % open new file called datadump_e where e is the current energy
    % iterator and give this the name 'f'.
    file_name = sprintf('datadump_%d.txt', S.e);
    f = fopen(file_name,'wt');
    
    % Write header: file number and date
    fprintf(f, 'dump.txt output file %d. \n', S.e);
    fprintf(f, 'date: %2s \n \n', S.time);
    %
    
    %% ---------------------- INTEGRATION VALUES ------------------------%
    %
    % Print Integration values
    fprintf(f, 'Smallest integration step, Delta T    = %2.3f \n', S.deltat);
    fprintf(f, 'Time interval between plotted points  = %2.3f \n', S.tout);
    fprintf(f, 'Number of points written to plot file = %d \n \n', S.Npoints);
    %
    %% ---------------------- LASPAR PARAMETERS ------------------------%
    %
    % Print physical parameters - not using smallest step (iqswitch)
    fprintf(f, 'Cavity length  = %3.1f \n', S.cavl);
    fprintf(f, 'Cavity area    = %2.4f \n', S.cmarea);
    fprintf(f, 'Rod length     = %3.2f \n', S.rodl);
    fprintf(f, 'Rod radius     = %3.2f \n', S.rodrad);
    fprintf(f, 'Pump abs. path = %3.2f \n', S.pabsl);
    fprintf(f, 'Q-switch delay = %2.2f \n', S.tq);
    fprintf(f, 'Q-switch flag  = %d \n \n', S.tqi);
    %fprintf(f, 'Smallest integration step, Delta T = %2.3f \n', iqswitch);
    % 
	%% ---------------------- PUMPING PARAMETERS ------------------------%
    %
    % Print pump parameters
    fprintf(f, 'Pump shape function         = %d \n', S.ipump);
    fprintf(f, 'Pump FWHM                   = %5.1f \n', S.fwhm);
    for i = 1:S.Npmplv
        fprintf(f, 'Pump abs. cross section     = %3.3e \n', S.pumpsig(i));
    end
    fprintf(f, 'Pump energy                 = %3.3f \n', S.pumpe);
    fprintf(f, 'Pump wavelength             = %3.3f \n', S.pumplam);
    fprintf(f, 'Pump scale factor, W1       = %2.3e \n', S.W1);
    fprintf(f, 'Net pump efficiency         = %1.2f \n', S.pumpeta);
    fprintf(f, 'Pump alignment efficiency   = %1.2f \n', S.pmpalgn);
    fprintf(f, 'Lower & upper manifolds     = ');
    for i = 1:S.Npmplv
        for j = 1:2
            fprintf(f, ' %d ', S.manpmp(i,j));
        end
    end
    fprintf(f, '\n \n');
    %
	%% ------------------------- ION PARAMETERS ----------------------------%
    %
    % Print ion/manifold parameters
    fprintf(f, 'Number of ions being considered = %d \n', S.Nions);
    for i = 1:S.Nions
        fprintf(f, 'Ion number %d (%s) has %d manifolds \n', i, S.ion(i), S.Nman(i));
    end
    for i = 1:S.mancount
        fprintf(f, 'Manifold number %d (%s) has an initial concentration = %2.3e \n', i, S.man(i), S.Y(i));
    end
    for i = 1:S.Nions
        fprintf(f, 'Ion number %d has a total concentration = %2.3e \n', i, S.totlcon(i));
    end
    %
	%% ------------------------- DECAY PARAMETERS --------------------------%
    %
    % Print decay parameters
    fprintf(f, '\nNumber of decay times = %d \n', S.Ndktime);
    for i = 1:S.Ndktime
        fprintf(f, 'Decay time %d goes from manifold %d to %d with rate = %2.3e \n', i, S.manup(i), S.manlo(i), S.dkrate(i));
    end
    %
	%% ------------------- BILINEAR PROCESS PARAMETERS ----------------------%
    %
    % Print bilinear transfer parameters
    fprintf(f, '\nNumber of bilinear processes = %d \n', S.Nbilin);
    for i = 1:S.Nbilin
        fprintf(f, 'Process %d: loss manifolds = %d %d, gain manifolds = %d %d. Rate = %2.3e. Acceptor = %d \n', ...
            i, S.ManLos1(i), S.ManLos2(i), S.ManGan1(i), S.ManGan2(i), S.bilinrt(i), S.Naceptr(i));
    end
    %
	%% --------------------------- TEMPERATURE -----------------------------%
    %
    % Print temperature
    fprintf(f, '\nTemperature = %2.1f \n \n', S.temper);
    %
	%% ---------------------- WAVELENGTH PARAMETERS ------------------------%
    %
    % Print wavelength parameters
    fprintf(f, 'Number of laser wavelengths = %d \n', S.Nlaslam);
    for i = 1:S.Nlaslam
        fprintf(f, '   Upper laser manifold = %d       Lower laser manifold = %d \n', S.ManUpLz(i), S.ManLoLz(i));
        fprintf(f, '   Emission wavelength                   = %2.3f \n', S.emislam(i));
        fprintf(f, '   Upper Stark fraction = %2.3f   Lower Stark fraction = %2.3f \n', S.FracUp(i), S.FracLo(i));
        fprintf(f, '   Stimulated emission cross section     = %2.3e \n', S.stim(i));
        fprintf(f, '   Spontaneous emission rate             = %2.3e \n', S.spon(i));
        fprintf(f, '   Bulk absorption coefficient (per cm)  = %2.3e \n', S.cavalfa(i));
        fprintf(f, '   Output coupler reflectivity           = %2.2f \n', S.Reflout(i));
        fprintf(f, '   Other mirror net reflectivity         = %2.2f \n', S.Reflbck(i));
        fprintf(f, '   Net intra-cavity surface transmission = %2.2f \n', S.surfacet(i));
        fprintf(f, '   Rod index of refraction               = %2.3f \n', S.roddex(i));
        fprintf(f, '   Cavity optical length (cm)            = %3.2f \n', S.optlng(i));
        fprintf(f, '   Spontaneous emission geometric factor = %2.3e \n', S.spo(i));
        fprintf(f, '   Cavity round trip time (micro sec)    = %2.3f \n', S.roundtt(i));
        fprintf(f, '   Threshold density assuming ion 1      = %2.3e \n \n', S.thresh(i));
    end
    
    %% ---------------------- TOP HALF CONCLUSION ------------------------%
    % The rest of dump it file is completed in DUMPIN2();
    fprintf(f,' \n \n ----------------------- \n \n \n');
    
end



%% ------------------------ PUMPVAL(ipump, fwhm, T) ------------------------%
% - ipump is the pump shape function, FWHM is the pump pulse duration, and
%   T is the current change in time.
% - Output a value for pump that fits these parameters. This value is given
%   by the chart below.
% - Structure of pumpval(...) is a large switch statement that is dependent
%   on the pump shape function. Smaller conditions are given inside each
%   case, which are usually dependent on FWHM, T.
function pump = pumpval(ipump, fwhm, T)
%% ----------------------- NORMALIZED PUMP FUNCTION ------------------------%
%
%                         IF  I=0  THEN CW PUMP
%                            1  THEN RECTANGULAR
%                            2  THEN  SIN(T)*SIN(T)
%                            3  THEN  T * EXP(-T*T)
%                            4  THEN  T * EXP(-T)
%                            5  THEN  1-T
%                            6  THEN   T
%                            7  THEN  SQUARE PULSE TRAIN WITH DUTY CYCLE 0.1 LENGTH FWHM
%                            8  THEN    "      "    "          "     "   0.5   "     "
%                          Read fortran code: 
%                          GO TO (label1, label2, …, labeln), integer-expression
%                          It is equivalent to this block IF statement:
% 
%                          IF (integer-expression .EQ. 1) THEN
%                             GO TO label1
%                          ELSE IF (integer-expression .EQ. 2) THEN
%                             GO TO label2
%                          …
%                          ELSE IF (integer-expression .EQ. n) THEN
%                             GO TO labeln
%                          END IF 
%
%% ---------------------- SWITCH CASE BASED ON IPUMP ------------------------%
    % 9 possiblities for ipump (= 0 - 9). See above for table
    switch ipump

        case 0 % ipump = 0: CW pump
            pump = 1/fwhm;

        case 1 % ipump = 1: Rectangular
            if(T > fwhm)
                pump = 0;
            else
                pump = 1/fwhm;
            end

        case 2 % ipump = 2: sin(T)*sin(T)
            width = 2*fwhm;
            if(T > width)
                pump = 0;
            else
                rarg = pi*T/width;
                pump = (sin(rarg)*sin(rarg))/fwhm;
            end

        case 3 % ipump = 3: T*exp(-T*T)
            if(fwhm <= 0)
                pump = 0;
            else
                alpha = 1.1331512/fwhm;
                arg = -(alpha*alpha*T*T);
                pump = 2*alpha*alpha*T*exp(arg);
            end

        case 4 % ipump = 4: T*exp(-T)
            if(fwhm <= 0)
                pump = 0;
            else
                alpha = 2.446386/fwhm;
                arg = -(alpha*T);
                pump = alpha*alpha*T*exp(arg);
            end

        case 5 % ipump = 5: 1 - T
            width = 2*fwhm;
            if(T > width)
                pump = 0;
            else
                pump = 1/fwhm - T/(2*fwhm*fwhm);
            end    

        case 6 % ipump = 6: T
            width = 2*fwhm;
            if(T > width)
                pump = 0;
            else
                pump = T/(2*fwhm*fwhm);
            end

        case 7 % ipump = 7: Square pulse train w/ duty cycle 0.1 len fwhm
            duty = 0.1;
            id = fix(T*16/fwhm);
            if((T >= (fwhm*id/16)) && (T < ((id + duty)*fwhm/16)))
                pump = 1/(duty*fwhm);
            elseif(T > fwhm)
                pump = 0;
            end

        case 8 % ipump = 8: Square pulse train w/ duty cycle 0.5 len fwhm
            duty = 0.1;
            id = fix(T*16/fwhm);
            if((T >= (fwhm*id/16)) && (T < ((id + duty)*fwhm/16)))
                pump = 1/(duty*fwhm);
            elseif(T > fwhm)
                pump = 0;
            end
    end
end



%% -------------------------- DUMPIN2(struct S) ---------------------------%
% - S is a structure of variables needed to compute and dump
% - Output is the same data dump file, .mat data as well
% - Data is just dumped into file line by line in a specific format similar
% to that of the original Fortran file
function dumpin2(S)
    %% ------------------------- APPEND TO FILE ---------------------------%
    % Append to same current data dump file
    % This section details where excitations occur...
    file_name = sprintf('datadump_%d.txt', S.e);
    f = fopen(file_name, 'a');
    
    %% ------------ INITIALIZE SOME VARIABLES & PRINT HEADER ----------------%
    % preallocate loss/gains and multiply reactions by stepsize
    [Rloss, Rgain] = deal(zeros(S.Nreact, S.mancount));
    for i = 1:S.Nreact
        S.react(i) = S.react(i)*S.deltat;
    end

    % three is an oddly specific number here - why?
    % output in cmd line is Eout(1,2,3)?
    Eout = zeros(3,1); % allocate array EOUT
    
    % Determine rod volume depending on end/side pump
    switch S.pumptype
        case 1 % pumptype == 1 (end resonator)
            rodvol = S.cmarea*S.rodl*S.pmpalgn;
        case 2 % pumptype == 2 (linear resonator)
            rodvol = S.rodarea*S.rodl*S.pmpalgn;
    end
    
    % Second half of dumpit file header
    fprintf(f, 'These are the contributions for each rate term. \n');
    fprintf(f, '   There are %d individual terms. \n \n', S.Nreact);
    
    %% ------------------- SINGLE ION DECAY PROCESSES ----------------------%
    % Print single ion decay process information
    fprintf(f, 'For single ion decay processes: \n');
    for i = 1:S.Ndktime
        % Reactions of decays are affected by rod volume.
        Rnumber = S.react(i)*rodvol;
        % Sum of the losses/gains
        Rloss(i, S.manup(i)) = Rloss(i, S.manup(i)) + Rnumber;
        Rgain(i, S.manlo(i)) = Rgain(i, S.manlo(i)) + Rnumber;
        % Print information
        fprintf(f, '   Process %d contributes %2.3e per cm^3 = %2.3e total ions.\n', i, S.react(i), Rnumber);
    end
    fprintf(f,' \n');
    
    %% ------------------- BILINEAR EXCHANGE PROCESSES ----------------------%
    % Find iteration interval (from Ndktime + 1 to Ndktime + Nbilin).
    NLO = S.Ndktime + 1;
    NHI = S.Ndktime + S.Nbilin;
    
    % Print process details
    fprintf(f, 'For bilinear exchange processes: \n');
    for i = NLO:NHI
        % K1 is the dummy ion number
        K1 = i - NLO + 1;
        % Reactions of processes affected by rod volume.
        Rnumber = S.react(i)*rodvol;
        % Sum of the reactions loss/gain
        Rloss(i, S.ManLos1(K1)) = Rloss(i, S.ManLos1(K1)) + Rnumber;
        Rgain(i, S.ManGan1(K1)) = Rgain(i, S.ManGan1(K1)) + Rnumber;
        Rloss(i, S.ManLos2(K1)) = Rloss(i, S.ManLos2(K1)) + Rnumber;
        Rgain(i, S.ManGan2(K1)) = Rgain(i, S.ManGan2(K1)) + Rnumber;
        % Print information
        fprintf(f, '   Process %d contributes %2.3e per cm^3 = %2.3e total ions.\n', K1, S.react(i), Rnumber);
    end
    fprintf(f,' \n');

    %% ------------------------- PUMPING PROCESSES --------------------------%
    EABS = 0;
    
    for i = 1:S.Npmplv
        % Lowest iterator required for pumping reactions
        NLO = S.Ndktime + S.Nbilin + i;
        % Print header
        fprintf(f, 'For the pumping process. \n');
        % Pumping reaction affected by rod volume and pumping alignment eff
        Rnumber = S.react(NLO)*rodvol/S.pmpalgn;
        % Calculate losses/gains
        Rloss(NLO, S.manpmp(i,1)) = Rloss(NLO, S.manpmp(i,1)) + Rnumber;
        Rgain(NLO, S.manpmp(i,2)) = Rgain(NLO, S.manpmp(i,2)) + Rnumber;
        % Calculate absorbed energy
        AbsEnergy = Rnumber*S.hc/S.pumplam;
        % Print information
        fprintf(f, '   Density absorbed = %2.3e \n', S.react(NLO));
        fprintf(f, '   Photons absorbed = %2.3e \n', Rnumber);
        fprintf(f, '   Energy absorbed  = %2.3e \n', AbsEnergy);
        % Sum to total absorbed energy
        EABS = AbsEnergy + EABS;
    end
    fprintf(f,'\n');

    %% --------------------- LASING-RELATED PROCESSES -----------------------%
    % NLO = Lowest value on interval
    NLO = S.Ndktime + S.Nbilin + S.Npmplv;
    
    % Print header
    fprintf(f, '%20s %28s %15s %14s \n','For the lasing process', 'Density', 'Number', 'Energy');
    for i = 1:S.Nlaslam
        fprintf(f, '   Wavelength = %2.3f \n', S.emislam(i));
        optvol = S.cmarea*S.optlng(i)*S.pmpalgn;
        %
        % ION LOSS FROM STRIMULATED EMISSION
        NLO = NLO + 1;
        Rnumber = S.react(NLO)*rodvol;
        energy = Rnumber*S.hc/S.emislam(i);
        fprintf(f, '   Lower level gain due to stim. emis.    = %2.3e      %2.3e      %2.3e \n', S.react(NLO), Rnumber, energy);
        %
        % PHOTON GAIN FROM STIMULATED EMISSION
        NLO = NLO + 1;
        Rnumber = S.react(NLO)*optvol;
        energy = Rnumber*S.hc/S.emislam(i);
        fprintf(f, '   Photon increase due to stim. emis.     = %2.3e      %2.3e      %2.3e \n', S.react(NLO), Rnumber, energy);
        %
        % LOSS FROM BULK CAVITY ABSORPTION
        NLO = NLO + 1;
        Rnumber = S.react(NLO)*optvol;
        energy = Rnumber*S.hc/S.emislam(i);
        fprintf(f, '   Bulk cavity absorption losses extract  = %2.3e      %2.3e      %2.3e \n', S.react(NLO), Rnumber, energy);
        %
        % GAIN FROM SPONTANEOUS EMISSION
        NLO = NLO + 1;
        Rnumber = S.react(NLO)*optvol;
        energy = Rnumber*S.hc/S.emislam(i);
        fprintf(f, '   Photon increase due to spon. emis.     = %2.3e      %2.3e      %2.3e \n', S.react(NLO), Rnumber, energy);
        %
        % LOSS FROM SURFACE TRANSMISSIONS
        NLO = NLO + 1;
        Rnumber = S.react(NLO)*optvol;
        energy = Rnumber*S.hc/S.emislam(i);
        fprintf(f, '   Surface transmission losses extract    = %2.3e      %2.3e      %2.3e \n', S.react(NLO), Rnumber, energy);
        %
        % LOSS FROM OTHER MIRRORS
        NLO = NLO + 1;
        Rnumber = S.react(NLO)*optvol;
        energy = Rnumber*S.hc/S.emislam(i);
        fprintf(f, '   Other mirror reflection losses extract = %2.3e      %2.3e      %2.3e \n', S.react(NLO), Rnumber, energy);
        %
        % LOSS FROM OUT COUPLING
        NLO = NLO + 1;
        Rnumber = S.react(NLO)*optvol;
        Eout(i) = Rnumber*S.hc/S.emislam(i);
        fprintf(f, '   Output coupling losses extract         = %2.3e      %2.3e      %2.3e \n', S.react(NLO), Rnumber, Eout(i));
        %
        % LOSS FROM Q-SWITCH
        NLO = NLO + 1;
        Rnumber = S.react(NLO)*optvol;
        energy = Rnumber*S.hc/S.emislam(i);
        fprintf(f, '   Q-switch losses extract                = %2.3e      %2.3e      %2.3e \n', S.react(NLO), Rnumber, energy);
    end
    fprintf(f,'\n');
    
    %% ----------------- PROCESS NUMBER CLASSIFICATIONS --------------------%
    % Single ion decay processes
    NLO = 1;
    NHI = S.Ndktime;
    fprintf(f,'Processes %d - %d are single ion decay. \n', NLO, NHI);
    %
    % Bilinear exchange processes
    NLO = NHI + 1;
    NHI = NHI + S.Nbilin;
    fprintf(f,'Processes %d - %d are bilinear exchange. \n', NLO, NHI);
    %
    % Pumping processes
    NLO = NHI + 1;
    NHI = NHI + S.Npmplv;
    fprintf(f,'Processes %d - %d are pumping. \n \n', NLO, NHI);

    %% ------------------- ION CONCENTRATION INCREASES ----------------------%
    % CURRENTLY NOT WORKING
    fprintf(f,'The budget for increasing ion concentrations. \n');
    %fprintf(f,'Process %d manifold %d manifold %d density %d number %d \n',);
    %fprintf(f,'%d %d %d %d %d \n');
    for i = 1:NHI
        for j = 1:S.mancount
            Rdensity = Rloss(i,j)/rodvol;
        end
    end
    fprintf(f,' \n');
    
    %% ------------------- FINAL POPULATION DENSITIES ----------------------%
    fprintf(f,'Final values for the population densities are... \n');
    for i = 1:S.NEQ-1
        fprintf(f,'   %d    %2s         %2.3e \n', i, S.man(i), S.Y(i));
    end
    fprintf(f, '\n');
    
    %% -------------------------- PULSE DATA ------------------------------%
    % Write plsmin
    fprintf(f,'Pulse starts at time = %2.3e \n', S.plsmin);
    % Write plsmax
    fprintf(f,'Pulse stops at time  = %2.3e \n', S.plsmax);
    % Write pwidth
    Pwidth = S.plsmax - S.plsmin;
    fprintf(f,'I.E. pulse width     = %2.3e \n', Pwidth);
    %
    
    %% ------------------------- TOTAL LOSSES ----------------------------%
    for i = 1:S.Nlaslam
        S.Reflout(i)  = S.Reflout(i)*S.roundtt(i);
        S.Reflbck(i)  = S.Reflbck(i)*S.roundtt(i);
        S.cavalfa(i)  = S.cavalfa(i)*S.roundtt(i);
        S.surfacet(i) = S.surfacet(i)*S.roundtt(i);
        xloss(i)      = S.Reflbck(i) + S.cavalfa(i) + S.surfacet(i);
    end
    
    %% ------------------- PRINT TO CMD LINE / PERFNUM ----------------------%
    % build matrix M holding perfnum data
    M = [S.temper; S.pumpe; EABS; S.Reflout(1); xloss(1); S.plsmin; Pwidth; Eout(1); Eout(2); Eout(3)];
    % print to command terminal
    fprintf("%8.1f %8.3f %8.4f %8.3f %8.3f %8.3f %10.2f %10.3f %10.3f %10.3f \n", M);
    %
    % print same data to perfnum file
    fprintf(S.perfnum, "%8.1f %8.3f %8.4f %8.3f %8.3f %8.3f %10.2f %10.3f %10.3f %10.3f \n", M);
    
    %% --------------------------- SAVE FINAL DATA --------------------------%
    save('dataoffdump2.mat');
    
    % Close current data dump file.
    fclose(f);
    
end

