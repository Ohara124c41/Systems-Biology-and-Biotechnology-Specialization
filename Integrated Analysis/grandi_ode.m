
close all
clear all

colors = repmat('krgbmc',1,500) ;


global R Frdy Temp FoRT Cmem Qpow 


global DcaJuncSL DcaSLcyto DnaJuncSL DnaSLcyto
global Vcell Vmyo Vsr Vsl Vjunc SAjunc SAsl 
global J_ca_juncsl J_ca_slmyo J_na_juncsl J_na_slmyo 

global Fjunc Fsl Fjunc_CaL Fsl_CaL 

global Cli Clo Ko Nao Cao Mgi 

global ecl 



global GNa GNaB IbarNaK KmNaip KmKo Q10NaK Q10KmNai 

global pNaK gkp 

global GClCa GClB KdClCa 

global pNa pCa pK Q10CaL 

global IbarNCX KmCai KmCao KmNai KmNao ksat nu Kdact 
global Q10NCX IbarSLCaP KmPCa GCaB Q10SLCaP 

global Q10SRCaP Vmax_SRCaP Kmf Kmr hillSRCaP ks 
global koCa kom kiCa kim ec50SR 

global Bmax_Naj Bmax_Nasl koff_na kon_na Bmax_TnClow koff_tncl kon_tncl 
global Bmax_TnChigh koff_tnchca kon_tnchca koff_tnchmg kon_tnchmg Bmax_CaM 
global koff_cam kon_cam Bmax_myosin koff_myoca kon_myoca koff_myomg kon_myomg 
global Bmax_SR koff_sr kon_sr Bmax_SLlowsl Bmax_SLlowj koff_sll kon_sll 
global Bmax_SLhighsl Bmax_SLhighj koff_slh kon_slh Bmax_Csqn koff_csqn kon_csqn 


global epi AF ISO RA
%% Model Parameters
%% EPI or ENDO?
epi=1;
%% AF
AF=0;
%% ISO
ISO=0;
%% Right ATRIUM
RA=0;



% Constants
R = 8314;       % [J/kmol*K]  
Frdy = 96485;   % [C/mol]  
Temp = 310;     % [K]
FoRT = Frdy/R/Temp;
Cmem = 1.1e-10;   % [F] membrane capacitance 1.3810e-10;%
Qpow = (Temp-310)/10;

% Cell geometry
cellLength = 100;     % cell length [um]113;%100
cellRadius = 10.25;   % cell radius [um]12;%10.25
junctionLength = 160e-3;  % junc length [um]
junctionRadius = 15e-3;   % junc radius [um]

DcaJuncSL = 1.64e-6;  % Dca junc to SL [cm^2/sec]
DcaSLcyto = 1.22e-6; % Dca SL to cyto [cm^2/sec]
DnaJuncSL = 1.09e-5;  % Dna junc to SL [cm^2/sec]
DnaSLcyto = 1.79e-5;  % Dna SL to cyto [cm^2/sec] 


Vcell = pi*cellRadius^2*cellLength*1e-15;    % [L]
Vmyo = 0.65*Vcell; 
Vsr = 0.035*Vcell; 
Vsl = 0.02*Vcell; 
Vjunc = 1*0.0539*.01*Vcell; 
SAjunc = 20150*pi*2*junctionLength*junctionRadius;  % [um^2]
SAsl = pi*2*cellRadius*cellLength;          % [um^2]


%J_ca_juncsl = DcaJuncSL*SAjunc/distSLcyto*1e-10;% [L/msec] = 1.1074e-13
%J_ca_slmyo = DcaSLcyto*SAsl/distJuncSL*1e-10;  % [L/msec] = 1.5714e-12
%J_na_juncsl = DnaJuncSL*SAjunc/distSLcyto*1e-10;% [L/msec] = 7.36e-13
%J_na_slmyo = DnaSLcyto*SAsl/distJuncSL*1e-10;  % [L/msec] = 2.3056e-11
%J_ca_juncsl = DcaJuncSL*SAjunc/distJuncSL*1e-10;% [L/msec] = 9.9664e-014
%J_ca_slmyo = DcaSLcyto*SAsl/distSLcyto*1e-10;  % [L/msec] = 1.7460e-012
%J_na_juncsl = DnaJuncSL*SAjunc/distJuncSL*1e-10;% [L/msec] = 6.6240e-013
%J_na_slmyo = DnaSLcyto*SAsl/distSLcyto*1e-10;  % [L/msec] = 2.5618e-011
% tau's from c-code, not used here
J_ca_juncsl =1/1.2134e12; % [L/msec] = 8.2413e-13
J_ca_slmyo = 1/2.68510e11; % [L/msec] = 3.2743e-12
J_na_juncsl = 1/(1.6382e12/3*100); % [L/msec] = 6.1043e-13
J_na_slmyo = 1/(1.8308e10/3*100);  % [L/msec] = 5.4621e-11

% Fractional currents in compartments
Fjunc = 0.11;   
Fsl = 1-Fjunc;
Fjunc_CaL = 0.9; 
Fsl_CaL = 1-Fjunc_CaL;

% Fixed ion concentrations     
Cli = 15;   % Intracellular Cl  [mM]
Clo = 150;  % Extracellular Cl  [mM]
Ko = 5.4;   % Extracellular K   [mM]
Nao = 140;  % Extracellular Na  [mM]
Cao = 1.8;  % Extracellular Ca  [mM]
Mgi = 1;    % Intracellular Mg  [mM]

% Reversal potential that does not change with time
ecl = (1/FoRT)*log(Cli/Clo);            % [mV]

% Na transport parameters
      
%%
GNa=23*(1-0.1*AF);  % [mS/uF]
GNaB = 0.597e-3;    % [mS/uF] 
IbarNaK = 1.26;     % [uA/uF]
KmNaip = 11*(1-0.25*ISO);         % [mM]11
KmKo =1.5;         % [mM]1.5
Q10NaK = 1.63;  
Q10KmNai = 1.39;

%% K current parameters
pNaK = 0.01833;      
gkp = 0.002;

% Cl current parameters
GClCa =0.0548;   % [mS/uF]
GClB = 9e-3;        % [mS/uF]
KdClCa = 100e-3;    % [mM]

% I_Ca parameters
pNa = (1+0.5*ISO)*(1-0.5*AF)*0.75e-8;       % [cm/sec]
pCa = (1+0.5*ISO)*(1-0.5*AF)*2.7e-4;       % [cm/sec]
pK = (1+0.5*ISO)*(1-0.5*AF)*1.35e-7;        % [cm/sec]
Q10CaL = 1.8;       

%% Ca transport parameters
IbarNCX = (1+0.4*AF)*3.15;      % [uA/uF]5.5 before - 9 in rabbit
KmCai = 3.59e-3;    % [mM]
KmCao = 1.3;        % [mM]
KmNai = 12.29;      % [mM]
KmNao = 87.5;       % [mM]
ksat = 0.27;        % [none]  
nu = 0.35;          % [none]
Kdact =0.384e-3;   % [mM] 0.256 rabbit
Q10NCX = 1.57;      % [none]
IbarSLCaP =  0.0471; % IbarSLCaP FEI changed [uA/uF](2.2 umol/L cytosol/sec) jeff 0.093 [uA/uF]
KmPCa =0.5e-3;     % [mM] 
GCaB = 6.0643e-4;    % [uA/uF] 3
Q10SLCaP = 2.35;    % [none]

% SR flux parameters
Q10SRCaP = 2.6;          % [none]
Vmax_SRCaP = 5.3114e-3;  % [mM/msec] (286 umol/L cytosol/sec)
Kmf = (2.5-1.25*ISO)*0.246e-3;          % [mM] default
Kmr = 1.7;               % [mM]L cytosol
hillSRCaP = 1.787;       % [mM]
ks = 25;                 % [1/ms]      
koCa = 10+20*AF+10*ISO*(1-AF);               % [mM^-2 1/ms]   %default 10   modified 20
kom = 0.06;              % [1/ms]     
kiCa = 0.5;              % [1/mM/ms]
kim = 0.005;             % [1/ms]
ec50SR = 0.45;           % [mM]

% Buffering parameters
% koff: [1/s] = 1e-3*[1/ms];  kon: [1/uM/s] = [1/mM/ms]
Bmax_Naj = 7.561;       % [mM] % Na buffering
Bmax_Nasl = 1.65;       % [mM]
koff_na = 1e-3;         % [1/ms]
kon_na = 0.1e-3;        % [1/mM/ms]
Bmax_TnClow = 70e-3;    % [mM]                      % TnC low affinity
koff_tncl = (1+0.5*ISO)*19.6e-3;    % [1/ms] 
kon_tncl = 32.7;        % [1/mM/ms]
Bmax_TnChigh = 140e-3;  % [mM]                      % TnC high affinity 
koff_tnchca = 0.032e-3; % [1/ms] 
kon_tnchca = 2.37;      % [1/mM/ms]
koff_tnchmg = 3.33e-3;  % [1/ms] 
kon_tnchmg = 3e-3;      % [1/mM/ms]
Bmax_CaM = 24e-3;       % [mM] **? about setting to 0 in c-code**   % CaM buffering
koff_cam = 238e-3;      % [1/ms] 
kon_cam = 34;           % [1/mM/ms]
Bmax_myosin = 140e-3;   % [mM]                      % Myosin buffering
koff_myoca = 0.46e-3;   % [1/ms]
kon_myoca = 13.8;       % [1/mM/ms]
koff_myomg = 0.057e-3;  % [1/ms]
kon_myomg = 0.0157;     % [1/mM/ms]
Bmax_SR = 19*.9e-3;     % [mM] (Bers text says 47e-3) 19e-3
koff_sr = 60e-3;        % [1/ms]
kon_sr = 100;           % [1/mM/ms]
Bmax_SLlowsl = 37.4e-3*Vmyo/Vsl;        % [mM]    % SL buffering
Bmax_SLlowj = 4.6e-3*Vmyo/Vjunc*0.1;    % [mM]    %Fei *0.1!!! junction reduction factor
koff_sll = 1300e-3;     % [1/ms]
kon_sll = 100;          % [1/mM/ms]
Bmax_SLhighsl = 13.4e-3*Vmyo/Vsl;       % [mM] 
Bmax_SLhighj = 1.65e-3*Vmyo/Vjunc*0.1;  % [mM] %Fei *0.1!!! junction reduction factor
koff_slh = 30e-3;       % [1/ms]
kon_slh = 100;          % [1/mM/ms]
Bmax_Csqn = 140e-3*Vmyo/Vsr;            % [mM] % Bmax_Csqn = 2.6;      % Csqn buffering
koff_csqn = 65;         % [1/ms] 
kon_csqn = 100;         % [1/mM/ms] 

gkr =0.035*sqrt(Ko/5.4);

gks_bar = 0.0035 ;
gks_junc=1*(1+1*AF+2*ISO)*gks_bar*1;
gks_sl=1*(1+1*AF+2*ISO)*gks_bar*1; %FRA

GtoFast=(1.0-0.7*AF)*0.165*1.0; %nS/pF maleckar; %human atrium
Gkur = 1*(1.0-0.5*AF)*(1+2*ISO)* 0.045*(1+0.2*RA); %nS/pF maleckar 0.045


% baselineparameters = [GNa ; GNaB ; gkr ; gks_bar ; GtoFast ; ...
%   Gkur ; IbarNaK ; gkp ; GClCa ; ...
%   GClB ; pNa ; pCa ; pK ; IbarNCX ; Vmax_SRCaP] ;
%       
% n_parameters = length(baselineparameters) ;
% variations = 5 ;
% actionpotentials = cell(2*variations,1) ;
% sigma = 0.05 ;     % standard deviation of variation of parameters
% allparameters = zeros(n_parameters,variations) ;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% Step 2:  Define simulation, stimulus, and recording parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


stiminterval = 1000 ;      % interval bewteen stimuli, ms
stimdelay = 100 ;
stimdur = 5 ;
stim_amp = 12.5 ;

n_stimuli = 1 ;
tend = n_stimuli*stiminterval ;              % end of simulation, ms

stim_starts = stimdelay + stiminterval*(0:n_stimuli-1)  ;
stim_ends = stim_starts + stimdur ;

simints = 3*n_stimuli ;
for i=1:n_stimuli
  intervals(3*i-2,:) = [stiminterval*(i-1),stim_starts(i)] ;
  intervals(3*i-1,:) = [stim_starts(i),stim_ends(i)] ;
  intervals(3*i,:) = [stim_ends(i),stiminterval*i] ;
end
intervals(end,:) = [stim_ends(end),tend] ;

Istim = zeros(simints,1) ;
stimindices = 3*(1:n_stimuli) - 1 ;
Istim(stimindices) = -stim_amp ;

% Initial conditions
% Obtained after steady-state pacing at 1 Hz

figure
handle = gcf ;
hold on


mo=1.405627e-3;
ho= 9.867005e-1;
jo=9.915620e-1; 
do=7.175662e-6; 
fo=1.000681; 
fcaBjo=2.421991e-2;
fcaBslo=1.452605e-2;
xtoso=4.051574e-3;
ytoso=9.945511e-1; 
xtofo=4.051574e-3; 
ytofo= 9.945511e-1; 
xkro=8.641386e-3; 
xkso= 5.412034e-3;
RyRro=8.884332e-1;
RyRoo=8.156628e-7; 
RyRio=1.024274e-7; 
NaBjo=3.539892;
NaBslo=7.720854e-1; 
TnCLo=8.773191e-3; 
TnCHco=1.078283e-1; 
TnCHmo=1.524002e-2; 
CaMo=2.911916e-4; 
Myoco=1.298754e-3; 
Myomo=1.381982e-1;
SRBo=2.143165e-3; 
SLLjo=9.566355e-3; 
SLLslo=1.110363e-1; 
SLHjo=7.347888e-3; 
SLHslo=7.297378e-2; 
Csqnbo= 1.242988;
Ca_sro=0.1e-1; %5.545201e-1; 
Najo=9.136;%8.80329; 
Naslo=9.136;%8.80733; 
Naio=9.136;%8.80853; 
Kio=120; 
Cajo=1.737475e-4; 
Caslo= 1.031812e-4; 
Caio=8.597401e-5; 
Vmo=-8.09763e+1; 
rtoso=0.9946; 
ICajuncinto=1; 
ICaslinto=0;
C1o=0.0015;       % [] 
C2o=0.0244;       % [] 
C3o=0.1494;       % [] 
C4o=0.4071;       % [] 
C5o=0.4161;       % [] 
C7o=0.0001;       % [] 
C8o=0.0006;       % [] 
C9o=0.0008;       % [] 
C10o=0;           % [] 
C11o=0;           % [] 
C12o=0;           % [] 
C13o=0;           % [] 
C14o=0;           % [] 
C15o=0;           % [] 
O1o=0;            % [] 
O2o=0;            % [] 
C6o=1-(C1o+C2o+C3o+C4o+C5o+C7o+C8o+C9o+C10o+C11o+C12o+C13o+C14o+C15o+O1o+O2o);       % []


%initial conditions for IKur
rkuro = 0;
skuro = 1.0;

% Gating variables      
%   1       2       3       4       5       6       7       8       9       10      11      12      13
%%   m       h       j       d       f       fcaBj   fcaBsl   xtos    ytos    xtof    ytof    xkr     xks   
%y10=[1.2e-3;0.99;   0.99;   0.0;    1.0;    0.0141; 0.0141;     0;      1;      0.0;    1.0;    0.0;    6e-3;];
y10=[mo; ho; jo; do; fo; fcaBjo; fcaBslo; xtoso; ytoso; xtofo; ytofo; xkro; xkso;];   
% RyR and Buffering variables
%   14      15      16      17      18      19      20      21      22      23      24
%%   RyRr    RyRo    RyRi    NaBj    NaBsl   TnCL    TnCHc   TnCHm   CaM     Myoc    Myom  
y20=[RyRro; RyRoo; RyRio; NaBjo; NaBslo; TnCLo; TnCHco; TnCHmo; CaMo; Myoco; Myomo;];           
%y20=[1;     0;      0;      1.8;   0.8;    0.012;   0.112;  0.01;   0.4e-3; 1.9e-3; 0.135;];
% More buffering variables
%   25      26      27      28      29      30
%%   SRB     SLLj   SLLsl    SLHj    SLHsl  Csqnb
y30=[SRBo; SLLjo; SLLslo; SLHjo; SLHslo; Csqnbo];
%y30=[3.3e-3; 0.012; 0.012; 0.13;  0.13;  1.5;];
%   Intracellular concentrations/ Membrane voltage
%    31      32      33      34      35      36      37     38     39    40   41
%%    Ca_sr   Naj     Nasl    Nai     Ki      Caj    Casl    Cai   Vm  rtos ?
y40=[Ca_sro; Najo; Naslo; Naio; Kio; Cajo; Caslo; Caio; Vmo; rtoso; 1]; 
y50=[C1o; C2o; C3o; C4o; C5o; C6o; C7o; C8o; C9o; C10o; C11o; C12o; C13o; C14o; C15o; O1o];
% y40=[0.9;    8.8;    8.8;    8.8;    135;    0.1e-3; 0.1e-3; 0.1e-3; -88;  0.89; 0;          0;];
% y50=[UIC3o; UIC2o; UIFo; UIM1o; UC3o; UC2o; UC1o; UOo; UIM2o; LC3o; LC2o; LC1o; LOo ];    

% SVP put initial variables for IKur; 11/12/09
%(y(58); y(59));
y60=[rkuro; skuro];

% EG put initial variables for INaL; 
%(y(60); y(61));
y70=[1; 0; 0];

% Put everything together
%y0  = [y10;y20;y30;y40;y50]    
y0  = [y10;y20;y30;y40;y50;y60;y70];

statevar_i = y0 ;

options = odeset('RelTol',1e-5,'MaxStep',1,'Stats','on'); 
    t = 0 ;
    statevars = statevar_i' ;
    for i=1:simints
      [post,posstatevars] = ode15s(@dydt_grandi,intervals(i,:),statevar_i,options,Istim(i)) ;
      t = [t;post(2:end)] ;
      statevars = [statevars;posstatevars(2:end,:)] ;
      statevar_i = posstatevars(end,:) ;
    end % for

outputcell = num2cell(statevars,1) ;
nvariables = length(outputcell) ;

V = statevars(:,39) ;

plot(t,V) 




