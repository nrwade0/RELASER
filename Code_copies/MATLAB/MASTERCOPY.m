clc, clear all, close all, format shorte;  
tic; 


% FIRST PROMPT AT START OF PROGRAM 
prompt={'Do you want to delete to previous record file "RELASER.txt"?'};
title='RELASER'; 
definput={'No'}; 
answer=newid(prompt,title,[1 40],definput); 
RECORD=answer{1};



%SWITCH TO DETERMINE WHETHER OR NOT A RECORD FILE NEEDS TO BE KEPT 
switch RECORD 
    case 'No'
        disp('Previous record file not deleted'); 
    case 'Yes'
        if (exist('RELASER.txt')==2) 
            delete 'RELASER.txt';  
        end 
        disp('Previous record file deleted');
end 

% RECORD FILE OF COMMAND WINDOW 
diary('RELASER.txt'); 

% OPEN EXCEL FILE TO SEE LASER PARAMETERS
parameters='parameters.xlsx';               
laspar=1;                  
specpar=2; 
[~,~,raw]=xlsread(parameters,laspar)     

% DISPLAY PARAMETERS SO USER CAN VIEW THEM FIRST BEFORE DECIDING TO CHANGE

disp('SUMMARY OF PARAMETERS')  
                             
Q_delay=xlsread(parameters,laspar,'C5:C5');
if (Q_delay<=0)
    disp('Q-SWITCH DELAY: No Q-Switch')
elseif (Q_delay>0)
    z=['Q-SWITCH DELAY: ', num2str(Q_delay)];
    disp(z)
end 

Q_interval=xlsread(parameters,laspar,'C6:C6'); 
if (Q_interval==0)
    disp('Q-INTERVAL: No Q-Switch Interval')
elseif (Q_interval~=0)
    z=['Q-INTERVAL: ', num2str(Q_interval)];
    disp(z)
end 

Resonator=xlsread(parameters,laspar,'C7:C7'); 
if (Resonator==1)
    disp('RESONATOR: RING (1)')
elseif (Resonator==2)
    disp('RESONATOR: LINEAR (2)')
end

Pumping=xlsread(parameters,laspar,'C8:C8');  
if (Pumping==1)
    disp('PUMPING: END (1)')
elseif (Pumping==2)
    disp('PUMPING: SIDE (2)')
end 

laspardata='Yes'; 

% CHUNK OF CODE THAT ASKS THE USER IF THEY WANT TO CHANGE ANY OF THE
% PARAMTERS
% BEGIN WHILE LOOP SO THAT IF THE USER DOES CHANGE A PARAMETER THEY ARE
% GIVEN THE OPTION TO THEN CHANGE ANOTHER, AND IF THEY DO NOT WANT TO THE
% PROGRAM WITH MOVE ON TO THE REST OF THE CODE 

while (laspardata=='Yes') 

prompt={'Do you want to change any parameters?'}; 
title='laspar'; 
definput={'No'};
answer=newid(prompt,title,[1 40],definput);
laspardata=answer{1}; 
switch laspardata
    case 'Yes'
        prompt={'Which entry do you want to change?'};
        title='laspar';
        dims=[1 40];
        definput={'Input the number corresponding to the entry'};

        answer=newid(prompt,title,dims,definput);
        user_val=str2num(answer{1}); 
       

            if (user_val==1)
            prompt={'Input new value for CAVITY LENGTH'};
            title='CL';
            dims=[1 30];
            definput={num2str(xlsread(parameters,laspar,'C1:C1'))};
            answer=newid(prompt,title,dims,definput);
            new_val=str2double(answer{1});
            xlswrite(parameters,new_val,laspar,'C1:C1');
            [~,~,raw]=xlsread(parameters,laspar)
            end

            if (user_val==2)
            prompt={'Input new value for CAVITY RADIUS'};
            title='CR';
            dim=[1 30];
            definput={num2str(xlsread(parameters,laspar,'C2:C2'))};
            anwser=newid(prompt,title,dims,definput);
            new_val=str2double(answer{1});

            xlswrite(parameters,new_val,laspar,'C2:C2');
            [~,~,raw]=xlsread(parameters,laspar)
            end

            if (user_val==3)
            prompt={'Input new value for ROD LENGTH(CM)'};
            title='RL';
            dim=[1 30];
            definput={num2str(xlsread(parameters,laspar,'C3:C3'))};
            answer=newid(prompt,title,dims,definput);
            new_val=str2double(answer{1});
            xlswrite(parameters,new_val,laspar,'C3:C3');
            [~,~,raw]=xlsread(parameters,laspar)
            end

            if (user_val==4)
            prompt={'Input new value for ROD RADIUS (CM)'};
            title='RR';
            dim=[1 30];
            definput={num2str(xlsread(parameters,laspar,'C4:C4'))};
            answer=newid(prompt,title,dims,definput);
            new_val=str2double(answer{1});
            xlswrite(parameters,new_val,laspar,'C4:C4');
            [~,~,raw]=xlsread(parameters,laspar)
            end

            if (user_val==5)
            prompt={'Input new value for Q-SWITCH DELAY'};
            title='QD';
            dim=[1 30];
            definput={num2str(xlsread(parameters,laspar,'C5:C5'))};
            answer=newid(prompt,title,dims,definput);
            new_val=str2double(answer{1});
            xlswrite(parameters,new_val,laspar,'C5:C5');
            [~,~,raw]=xlsread(parameters,laspar)
            end

            if (user_val==6)
            prompt={'Input new value for Q-SWITCH INTERVAL'};
            title='QI';
            dim=[1 30];
            definput={num2str(xlsread(parameters,laspar,'C6:C6'))};
            answer=newid(prompt,title,dims,definput);
            new_val=str2double(answer{1});
            xlswrite(parameters,new_val,laspar,'C6:C6');
            [~,~,raw]=xlsread(parameters,laspar)
            end

            if (user_val==7)
            prompt={'Input new value for RESONATOR'};
            title='R';
            dim=[1 30];
            definput={num2str(xlsread(parameters,laspar,'C7:C7'))};
            answer=newid(prompt,title,dims,definput);
            new_val=str2double(answer{1});
            xlswrite(parameters,new_val,laspar,'C7:C7');
            [~,~,raw]=xlsread(parameters,laspar)
            end

            if (user_val==8)
            prompt={'Input new value for PUMPING'};
            title='P';
            dim=[1 30];
            definput={num2str(xlsread(parameters,laspar,'C8:C8'))};
            answer=newid(prompt,title,dims,definput);
            new_val=str2double(answer{1});
            xlswrite(parameters,new_val,laspar,'C8:C8');
            [~,~,raw]=xlsread(parameters,laspar)
            end

            if (user_val==9)
            prompt={'Input new value for STEP SIZE'};
            title='SS';
            dim=[1 30];
            definput={num2str(xlsread(parameters,laspar,'C9:C9'))};
            answer=newid(prompt,title,dims,deinput);
            new_val=str2double(answer{1});
            xlswrite(parameters,new_val,laspar,'C9:C9');
            [~,~,raw]=xlsread(parameters,laspar)
            end

            if (user_val==10)
            prompt={'Input new value for NUMNBER OF POINTS'};
            title='NP';
            dim=[1 30];
            definput={num2str(xlsread(parameters,laspar,'C10:C10'))};
            answer=newid(prompt,title,dims,definput);
            new_val=str2double(answer{1});
            xlswrite(parameters,new_val,laspar,'C10:C10');
            [~,~,raw]=xlsread(parameters,laspar)
            end

            if (user_val==11)
            prompt={'Input new value for TIME INTERVAL BETWEEN POINTS'};
            title='TI';
            dim=[1 30];
            definput={num2str(xlsread(parameters,laspar,'C11:C11'))};
            answer=newid(prompt,title,dims,definput);
            new_val=str2double(answer{1});
            xlswrite(parameters,new_val,laspar,'C11:C11');
            [~,~,raw]=xlsread(parameters,laspar)
    end
    case 'No'
        break 
end 
end 


% HERE WE BEGIN GETTING ALL OF OUR DUCKS IN A ROW AND DEFINING ALL
% VARIABLES BEFORE MOVING ON THE EQUATIONS

% BEGIN DEFINING AND EXTRACTING ALL VARIABLES FROM EXCEL FILE AFTER THE
% USER HAS MADE ALL OF THEIR CHANGES TO THE LASER PARAMETER FILE. ALL
% VARIABLES COMING FROM EXCEL FILE WILL BE IN CAPTIAL LETTERS TO STAY
% SIMILAR TO FORTRAN PROGRAM THAT THIS PROGRAM IS BASED OFF OF. 
% ALL UNITS ARE IN METERS AND SECONDS 

%----------------------------LASPAR PARAMETERS----------------------------%
CAVL=xlsread(parameters,laspar,'C1:C1');    % [cm] CAVITY LENGTH 
CMRAD=xlsread(parameters,laspar,'C2:C2');   % [cm] CAVITY MODE RADIUS 
RODL=xlsread(parameters,laspar,'C3:C3');    % [cm] ROD LENGHT   
RODRAD=xlsread(parameters,laspar,'C4:C4');  % [cm] ROD RADIUS 
TQ=xlsread(parameters,laspar,'C5:C5');      % [us] Q-SWITCH DELAY  
IQ=xlsread(parameters,laspar,'C6:C6');      % [dimensionless] Q-SWITCH INTERVAL 
RESTYPE=xlsread(parameters,laspar,'C7:C7'); % [dimensionless] LINEAR(1) OR (2) RING RESONATOR  
PUMPTYPE=xlsread(parameters,laspar,'C8:C8');% [dimensionless] END (1) OR (2) SIDE PUMPING  
DELTAT=xlsread(parameters,laspar,'C9:C9');  % [us] STEP SIZE 
NPOINTS=xlsread(parameters,laspar,'C10:C10'); % [dimensionless] TOTAL NUMBER OF POINTS  
TOUT=xlsread(parameters,laspar,'C11:C11');  % [us] TIME INTERVAL BETWEEN POINTS  


%---------------------------SPECPAR PARAMETERS----------------------------%
NIONS=xlsread(parameters,specpar,'B2:B2');  % [dimensionless] NUMBER OF IONS   
CONS=xlsread(parameters,specpar,'E5:E8');   % [vector] MANIFOLD CONCENTRATIONS
NMAN=xlsread(parameters,specpar,'B3:B3');   % [dimensionless] NUMBER OF MANIFOLD BEING CONSIDERED OF ION TYPE
DKTIME=xlsread(parameters,specpar,'B11:B14'); % [vector] DECAY TIMES 
NDKTIME=xlsread(parameters,specpar,'B10:B10'); % [dimensionless] THE NUMBER OF DECAY TIMES 
NBILIN=xlsread(parameters,specpar,'B15:B15');  % Number of bilinear processes
UPSTARK=xlsread(parameters,specpar,'B20:B20'); % [dimensionless] UPPER STARK POPULATION 
LOSTARK=xlsread(parameters,specpar,'B21:B21'); % [dimensionless] LOWER STARK POPULATION
SPONLIFE=xlsread(parameters,specpar,'B23:B23'); % [us] SPONTANEOUS EMISSION LIFETIME
CAVALFA=xlsread(parameters,specpar,'B24:B24'); % [cm] BULK CAVITY LOSSES (PER CM)
REFLOUT=xlsread(parameters,specpar,'B25:B25'); % [dimensionless] OUTPUT COUPLER REFLECTIVITY 
REFLBCK=xlsread(parameters,specpar,'B26:B26'); % [dimensionless] OTHER MIRROR NET REFLECTIVITY 
SURFACET=xlsread(parameters,specpar,'B27:B27'); % [dimensionless] INTRA-CAVITY SURFACE NET TRANSMISSION 
RODDEX=xlsread(parameters,specpar,'B28:B28'); % [dimensionless] ROD INDEX OF REFRACTION 
PUMPLAM=xlsread(parameters,specpar,'B30:B30'); % [cm] PUMP WAVELENGTH 
MANAPL=xlsread(parameters,specpar,'B31:B31'); % [dimensionless] NUMBER OF MANIFOLDS CAPABLE OF ABSORBING LIGHT 
PUMPE=xlsread(parameters,specpar,'B33:B33'); % [mJ] AMOUNT OF ENERGY THE PUMP HAS 
FWHM=xlsread(parameters,specpar,'B34:B34');  % [us] FULL WIDTH HALF MAX (ALSO KNOWN AS PUMP PULSEWIDTH IN NASH CODE)
IPUMP=xlsread(parameters,specpar,'B35:B35'); % [dimensionless] PUMP SHAPE FUNCTION 
PUMPETA=xlsread(parameters,specpar,'B36:B36'); % [dimensionless] NET PUMP COUPLING EFFICIENCY 
NPMPLV=xlsread(parameters,specpar,'B37:B37'); % [dimensionless] PUMP ALIGNMENT EFFICIENCY


% OTHER VARIABLES THAT WERE NOT FOUND IN THE LASPAR AND SPECPAR SHEETS 
% PI=pi FOR CONVIENCE IN CAPS LOCK TYPING 
PI=pi; 
CMAREA=PI*CMRAD^2;                  % [cm^3] CAVITY MODE AREA 
RODAREA=PI*RODRAD^2;                % [cm^2] ROD AREA
RODVOL=RODAREA*RODL;                % [cm^3] ROD VOLUME 

OPTLNG=(RODDEX-1)*RODL+CAVL;        % [cm] CAVITY OPTICAL PATHLENGTH
% OPTLNG=RODDEX*RODL+CAVL-RODL (*from fortran program, ask Dr. Walsh)

FRACUP=0.0874;                      % [dimensionless] UPPER LASER LEVEL POPULATION (boltzmann f7 in paper)
FRACLO=0.0238;                      % [dimensionless] LOWER LASER LEVEL POPULATION (boltzmann f8 in paper)
STIM=1.6E-20;                       % [dimensionless] STIMULATED EMISSION (this is sigma_e in the paper)
SIGMA_SE=STIM/FRACUP;               % [cm^2] STIMULATED EMISSION CROSS SECTION
PUMPSIG=8.36E-21;                   % [cm^2] PUMP ABSORPTION CROSS SECTION 
ETA_P=1.3;                          % [dimensionless] PUMP DISTRIBUTION EFFICIENCY (FUDGE FACTOR) 
LASLAM=2.063E-4;                    % [cm] LASER WAVELENGTH 
R_PASSIVE=0.98;                     % [dimensionless] PASSIVE CAVITY LOSSES 
TS=1;                               % [dimensionless] LASER ROD SURFACE TRANSMISSION
CHo=5E-3;                           % [dimensionless] CONCENTRATION OF HOLMIUM IS 0.5%
rho_LuLF = 6.170;                   % [g/cm^3], density of LuLF crystal 
mLu = 174.967;                      % [amu], atomic mass of Lutetium 
mLi = 6.941;                        % [amu], atomic mass of Lithium
mF = 18.9984032;                    % [amu], atomic mass of Fluorine
mp = 1.6726e-24;                    % [g], proton mass
Lunum = 1;                          % [dimensionless], number of Lutetium ions in LuLiF unit
Linum = 1;                          % [dimensionless], number of Lithium ions in LuLiF unit
Fnum = 4;                           % [dimensionless], number of Fluorine ions in LuLiF unit
NS=1.441E22;                        % [cm^-3] SITE DENSITY 
n8_INIT=NS*CHo;                     % [cm^-3] INITIAL CONCENTRATION IN GROUNDSTATE MANIFOLD n8


% TRANSFER RATES [cm^3/us]
p77=3.2E-24;
p58=4.9E-24; 

% LIFETIMES [us]
t7=16000;
t6=2200;
t5=20;

%BRANCHING RATIOS/PROBABILITIES
B56=1;
B67=1;

%PHYSICAL CONSTANTS 
H = 6.626e-25;                      % [mJ*us] planck's constant 
C = 2.99792e4;                      % [cm/us] speed of light 

CAVDIA=4;
CAVAREA=PI*CAVDIA^2/4; 
% IF STATEMENT FOR RESONATOR TYPE THAT DETERMINES WHICH B (SPONTANEOUS
% EMISSION) WILL BE USED 
if RESTYPE==1 
    B=CMAREA/(4*PI*OPTLNG^2);      % [dimensionless] SPONTANEOUS EMISSION 
elseif RESTYPE==2 
    B=RODAREA/(4*PI*OPTLNG^2);      % [dimensionless] SPONTANEOUS EMISSION 
end 

SPO=B*RODL/OPTLNG;

% HERE WE BEGIN BY BUILDING THE EQUATIONS WITH THE DEFINED VARIABLES ABOVE
Rp=(ETA_P*PUMPLAM*PUMPE)/(H*C*PI*CMRAD^2*RODL*FWHM); % [cm^-3*us^-1] pump rate 

CAVLIFE=1/((C/(RESTYPE*OPTLNG))*(log(REFLOUT*R_PASSIVE)+2*(RESTYPE^2)*(1-TS))); %[us] CAVITY LIFETIME 

% OTHER VARIABLES AND CAVITY LIFETIMES FROM NASH'S PROGRAM, DOES NOT CHANGE
% RESULTS IN TESTING, SO FAR
% PASSIVE_LOSS=0;
% CAVLIFE=1/((C/(RESTYPE*OPTLNG))*(-log(REFLOUT)-log(R_PASSIVE)-log((1-PASSIVE_LOSS)^2)));

ROUNDTT=(RESTYPE*OPTLNG)/C;         % CAVITY ROUNDTRIP TIME [us]

V_RES=PI*CMRAD^2*OPTLNG;            % [cm^3] RESONATOR MODE VOLUME 

PHI_INIT=1/V_RES;                   % [cm^-3] INITIAL PHOTON DENSITY 

T_MAX=1E3;                         % [us] TOTAL CAVITY EVOLUTION TIME 
TT=(0:DELTAT:T_MAX)';               % TOTAL TIME DIVIDED INTO NPOINTS FROM LASPAR 
% THIS WAS LEFT OUT DUE TO THE REDUNDANCE OF IT
% ITS PURPOSE IS FOR pump_rate IN THE RATE EQUATIONS PROGRAM, IF YOU USE
% THAT, THEN USE THIS LINE OF CODE AND CHANGE R_p TO Rp_LEVEL
% Rp=(TT<=1E4)*Rp_LEVEL 

%---------------FORMATION OF THE DIFFERENTIAL EQUATIONS-------------------%


% FOR A START WE BEGIN BY DEFINING THE INTIAL CONDITIONS OF THE RATE
% EQUATIONS 

y_init=[0,0,0,n8_INIT, PHI_INIT];

% SETS THE TOLERANCES OF THE ODE SOLVER
options=odeset('RelTol',1e-3,'AbsTol',1e-6,'OutputFcn',@odeplot,'OutputSel',[1 2 3 4 5],'NonNegative',[1 2 3 4 5],'InitialStep',DELTAT,'MaxStep',NPOINTS); 

% VECTOR THAT CONTAINS VARIABLES THAT WILL BE USED IN DIFF EQS. 
p=[Rp,p77,p58,t5,t6,t7,B56,B67,FRACUP,FRACLO,OPTLNG,PUMPSIG,SIGMA_SE,CAVLIFE,B,REFLOUT,ROUNDTT,RODL,C,FWHM,RODDEX,n8_INIT];

clear title; % PLEASE DO NOT DELETE 

% ACTUAL SOLVING OF THE DIFFERENTIAL EQUATIONS 
% FIGURE COMMAND MAKES GRAPH OF SOLVING OF DIFF EQS.

figure; 
[t,y,pump_rate]=ode15s('rate_equations_end_pump',[0 T_MAX],y_init,options,p,TT,Rp); 
savefig('RateEQ.fig');

n7=y(:,1); 
n6=y(:,2); 
n5=y(:,3); 
in8=y(:,4); 
phi=y(:,5); 

n8=n8_INIT-n7-n6-n5;

Eout_paper = ((H*C)/LASLAM)*PI*(CMRAD^2)*OPTLNG*phi(end);  %  [J], OUTPUT ENERGY OF THE PAPER

yTotal = y(:,1:4); 
test = sum(yTotal');
yHo = y(:,1:4);
testHo = sum(yHo');

disp('                                                                ');
disp('**************************RESULTS*******************************');
disp(['Laser Output Pulse Energy = ',num2str(Eout_paper*1000),'mJ']);
disp('                                                                ');

figure; subplot(4,1,1); plot(TT,Rp); legend('pump'); title(['v3 EndPump, CAVALFA=',num2str(CAVALFA)]);
subplot(4,1,2); plot(t,FRACUP*n7,t,FRACLO*n8); legend('f7n7','f8n8');
subplot(4,1,3); plot(t,FRACUP*n7-FRACLO*n8); legend('f7n7-f8n8');
subplot(4,1,4); plot(t,phi); legend('phi');
savefig('subplots.fig'); 

figure; plot(t,n5); title('n5'); 
figure; plot(t,n6); title('n6'); 
figure; plot(t,n7); title('n7'); 
figure; plot(t,in8); title('n8'); 
figure; plot(t,phi); title('n9'); 

toc

disp('                                                                ');


% TURNS OFF THE RECORDING OF THE COMMAND WINDOW
diary off 

% FINAL PROMPT TO SEE IF THE USER WANTS TO RUN THE PROGRAM AGAIN 
%(MAYBE FOR A TRIAL WITH DIFFERENT PARAMETERS)

prompt={'Do you want to run the program again?'};
title='Run?';
definput={'No'};
answer=newid(prompt,title,[1 40],definput);
RUNAGAIN=answer{1}; 

% SWITCH COMMAND THAT IS BASED OFF OF THE ANSWER TO THE PROMPT ABOVE

switch RUNAGAIN
    case 'Yes'
        run('MASTERCOPY.m')
    case 'No'
end 
    
% FINAL CODE THAT TELLS HOW LONG ENTIRE PROCESS TOOK (INCLUDING TIME TAKEN
% TO ANSWER THE PROMPTS)
toc;




