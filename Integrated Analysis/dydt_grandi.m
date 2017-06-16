function deriv = dydt_grandi(t,y,I_app)

global epi AF ISO RA

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


% Nernst Potentials
ena_junc = (1/FoRT)*log(Nao/y(32));     % [mV]
ena_sl = (1/FoRT)*log(Nao/y(33));       % [mV]
ek = (1/FoRT)*log(Ko/y(35));	        % [mV]
eca_junc = (1/FoRT/2)*log(Cao/y(36));   % [mV]
eca_sl = (1/FoRT/2)*log(Cao/y(37));     % [mV]

%% Membrane Currents
mss = 1 / ((1 + exp( -(56.86 + y(39)) / 9.03 ))^2);
taum = 0.1292 * exp(-((y(39)+45.79)/15.54)^2) + 0.06487 * exp(-((y(39)-4.823)/51.12)^2);                 
 
ah = (y(39) >= -40) * (0)... 
   + (y(39) < -40) * (0.057 * exp( -(y(39) + 80) / 6.8 )); 
bh = (y(39) >= -40) * (0.77 / (0.13*(1 + exp( -(y(39) + 10.66) / 11.1 )))) ...
   + (y(39) < -40) * ((2.7 * exp( 0.079 * y(39)) + 3.1*10^5 * exp(0.3485 * y(39)))); 
tauh = 1 / (ah + bh); 
hss = 1 / ((1 + exp( (y(39) + 71.55)/7.43 ))^2);
 
aj = (y(39) >= -40) * (0) ...
    +(y(39) < -40) * (((-2.5428 * 10^4*exp(0.2444*y(39)) - 6.948*10^-6 * exp(-0.04391*y(39))) * (y(39) + 37.78)) / ...
                     (1 + exp( 0.311 * (y(39) + 79.23) )));
bj = (y(39) >= -40) * ((0.6 * exp( 0.057 * y(39))) / (1 + exp( -0.1 * (y(39) + 32) ))) ...
   + (y(39) < -40) * ((0.02424 * exp( -0.01052 * y(39) )) / (1 + exp( -0.1378 * (y(39) + 40.14) ))); 
tauj = 1 / (aj + bj);
jss = 1 / ((1 + exp( (y(39) + 71.55)/7.43 ))^2);         
 
ydot(1) = (mss - y(1)) / taum;
ydot(2) = (hss - y(2)) / tauh;
ydot(3) = (jss - y(3)) / tauj;
    
I_Na_junc = Fjunc*GNa*y(1)^3*y(2)*y(3)*(y(39)-ena_junc);
I_Na_sl = Fsl*GNa*y(1)^3*y(2)*y(3)*(y(39)-ena_sl);
I_Na = I_Na_junc+I_Na_sl;


% Late I_Na
GNaL=0.0025*AF;
aml = 0.32*(y(39)+47.13)/(1-exp(-0.1*(y(39)+47.13)));
bml = 0.08*exp(-y(39)/11);
hlinf = 1/(1+exp((y(39)+91)/6.1));
tauhl=600;
ydot(60) = aml*(1-y(60))-bml*y(60);
ydot(61) = (hlinf-y(61))/tauhl;

I_NaL_junc = Fjunc*GNaL*y(60)^3*y(61)*(y(39)-ena_junc);
I_NaL_sl = Fsl*GNaL*y(60)^3*y(61)*(y(39)-ena_sl);
I_NaL = I_NaL_junc + I_NaL_sl;
if t<9050
    ydot(62)=0;
else
    ydot(62)=I_NaL;
end
% I_nabk: Na Background Current
I_nabk_junc = Fjunc*GNaB*(y(39)-ena_junc);
I_nabk_sl = Fsl*GNaB*(y(39)-ena_sl);
I_nabk = I_nabk_junc+I_nabk_sl;

% I_nak: Na/K Pump Current
sigma = (exp(Nao/67.3)-1)/7;
fnak = 1/(1+0.1245*exp(-0.1*y(39)*FoRT)+0.0365*sigma*exp(-y(39)*FoRT));
I_nak_junc = 1*Fjunc*IbarNaK*fnak*Ko /(1+(KmNaip/y(32))^4) /(Ko+KmKo);
I_nak_sl = 1*Fsl*IbarNaK*fnak*Ko /(1+(KmNaip/y(33))^4) /(Ko+KmKo);
I_nak = I_nak_junc+I_nak_sl;

%% I_kr: Rapidly Activating K Current
gkr =0.035*sqrt(Ko/5.4);
xrss = 1/(1+exp(-(y(39)+10)/5));
tauxr = 550/(1+exp((-22-y(39))/9))*6/(1+exp((y(39)-(-11))/9))+230/(1+exp((y(39)-(-40))/20));
ydot(12) = (xrss-y(12))/tauxr;
rkr = 1/(1+exp((y(39)+74)/24));
I_kr = gkr*y(12)*rkr*(y(39)-ek);


%% I_ks: Slowly Activating K Current
markov_iks=0;
% pcaks_junc = -log10(y(36))+3.0; 
% pcaks_sl = -log10(y(37))+3.0;  
% gks_junc = 0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_junc)/0.6)));
% gks_sl = 0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_sl)/0.6)));     

eks = (1/FoRT)*log((Ko+pNaK*Nao)/(y(35)+pNaK*y(34)));

if markov_iks==0;
gks_junc=1*(1+1*AF+2*ISO)*0.0035*1;
gks_sl=1*(1+1*AF+2*ISO)*0.0035*1; %FRA
xsss = 1 / (1+exp(-(y(39)+40*ISO + 3.8)/14.25)); % fitting Fra
tauxs=990.1/(1+exp(-(y(39)+40*ISO+2.436)/14.12));
ydot(13) = (xsss-y(13))/tauxs;
I_ks_junc = Fjunc*gks_junc*y(13)^2*(y(39)-eks);
I_ks_sl = Fsl*gks_sl*y(13)^2*(y(39)-eks);                                                                                                                                   
I_ks = I_ks_junc+I_ks_sl;
else
    gks_junc=1*0.0065;
    gks_sl=1*0.0065; %FRA
    alpha=3.98e-4*exp(3.61e-1*y(39)*FoRT);
    beta=5.74e-5*exp(-9.23e-2*y(39)*FoRT);
    gamma=3.41e-3*exp(8.68e-1*y(39)*FoRT);
    delta=1.2e-3*exp(-3.3e-1*y(39)*FoRT);
    teta=6.47e-3;
    eta=1.25e-2*exp(-4.81e-1*y(39)*FoRT);
    psi=6.33e-3*exp(1.27*y(39)*FoRT);
    omega=4.91e-3*exp(-6.79e-1*y(39)*FoRT);
    
ydot(42)=-4*alpha*y(42)+beta*y(43);
ydot(43)=4*alpha*y(42)-(beta+gamma+3*alpha)*y(43)+2*beta*y(44);
ydot(44)=3*alpha*y(43)-(2*beta+2*gamma+2*alpha)*y(44)+3*beta*y(45);
ydot(45)=2*alpha*y(44)-(3*beta+3*gamma+alpha)*y(45)+4*beta*y(46);
ydot(46)=1*alpha*y(44)-(4*beta+4*gamma)*y(46)+delta*y(50);    
ydot(47)=gamma*y(43)-(delta+3*alpha)*y(47)+beta*y(48);   
ydot(48)=2*gamma*y(44)+3*alpha*y(47)-(delta+beta+2*alpha+gamma)*y(48)+2*beta*y(49)+2*delta*y(51);
ydot(49)=3*gamma*y(45)+2*alpha*y(48)-(delta+2*beta+1*alpha+2*gamma)*y(49)+3*beta*y(50)+2*delta*y(52);
ydot(50)=4*gamma*y(46)+1*alpha*y(49)-(delta+3*beta+0*alpha+3*gamma)*y(50)+2*delta*y(53);
ydot(51)=1*gamma*y(48)-(2*delta+2*alpha)*y(51)+beta*y(52);  
ydot(52)=2*gamma*y(49)+2*alpha*y(51)-(2*delta+beta+1*alpha+gamma)*y(52)+2*beta*y(53)+3*delta*y(54);
ydot(53)=3*gamma*y(50)+1*alpha*y(52)-(2*delta+2*beta+2*gamma)*y(53)+3*delta*y(55);
ydot(54)=1*gamma*y(52)-(3*delta+1*alpha)*y(54)+beta*y(55);  
ydot(55)=2*gamma*y(53)+1*alpha*y(54)-(3*delta+1*beta+1*gamma)*y(55)+4*delta*y(56);
ydot(56)=1*gamma*y(55)-(4*delta+teta)*y(56)+eta*y(57);
O2=1-(y(42)+y(43)+y(44)+y(45)+y(46)+y(47)+y(49)+y(48)+y(50)+y(51)+y(52)+y(53)+y(54)+y(55)+y(56)+y(57));
ydot(57)=1*teta*y(56)-(eta+psi)*y(57)+omega*O2;
I_ks_junc = Fjunc*gks_junc*(y(57)+O2)*(y(39)-eks);
I_ks_sl = Fsl*gks_sl*(y(57)+O2)*(y(39)-eks);                                                                                                                                   
I_ks = I_ks_junc+I_ks_sl;
end
%I_kp: Plateau K current
kp_kp = 1/(1+exp(7.488-y(39)/5.98));
I_kp_junc = Fjunc*gkp*kp_kp*(y(39)-ek);
I_kp_sl = Fsl*gkp*kp_kp*(y(39)-ek);
I_kp = I_kp_junc+I_kp_sl;

%% I_to: Transient Outward K Current (slow and fast components)
% modified for human myocytes

GtoFast=(1.0-0.7*AF)*0.165*1.0; %nS/pF maleckar; %human atrium

%11/12/09; changed Itof to that from maleckar/giles/2009; removed I_tos
%atrium
%equations for activation; 
xtoss = ( (1)./ ( 1 + exp( -(y(39)+1.0)/11.0 ) ) );
tauxtof = 3.5*exp(-((y(39)/30.0)^2.0))+1.5;

%equations for inactivation;
ytoss = ( (1.0)./ ( 1 + exp( (y(39)+40.5)/11.5) ) ) ;
tauytof =25.635*exp(-(((y(39)+52.45)/15.8827)^2.0))+24.14;%14.14

ydot(10) = (xtoss-y(10))/tauxtof;
ydot(11) = (ytoss-y(11))/tauytof;
I_tof = 1.0*GtoFast*y(10)*y(11)*(y(39)-ek);

I_to = 1*I_tof;

%% I_kur: Ultra rapid delayed rectifier Outward K Current
%Equation for IKur; from Maleckar et al. 2009 - EG
%atrium
%equations for activation;
Gkur = 1*(1.0-0.5*AF)*(1+2*ISO)* 0.045*(1+0.2*RA); %nS/pF maleckar 0.045
xkurss = ( (1)./ ( 1 + exp( (y(39)+6)/-8.6 ) ) );
tauxkur = 9/(1+exp((y(39)+5)/12.0))+0.5;

%equations for inactivation;
ykurss = ( (1)./ ( 1 + exp( (y(39)+7.5)/10 ) ) );
tauykur = 590/(1+exp((y(39)+60)/10.0))+3050;

ydot(58) = (xkurss-y(58))/tauxkur;
ydot(59) = (ykurss-y(59))/tauykur;
I_kur = 1*Gkur*y(58)*y(59)*(y(39)-ek);

%% I_ki: Time-Independent K Current
aki = 1.02/(1+exp(0.2385*(y(39)-ek-59.215)));
bki =(0.49124*exp(0.08032*(y(39)+5.476-ek)) + exp(0.06175*(y(39)-ek-594.31))) /(1 + exp(-0.5143*(y(39)-ek+4.753)));
kiss = aki/(aki+bki);

%I_ki =1* 0.35*sqrt(Ko/5.4)*kiss*(y(39)-ek);
%SVP 11/11/09
%multiplieD IK1 by 0.15 to scale it to single cell isolated atrial cell
%resting potential
I_ki =(1+1*AF)*0.0525*sqrt(Ko/5.4)*kiss*(y(39)-ek);

% I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/y(36))*(y(39)-ecl);
I_ClCa_sl = Fsl*GClCa/(1+KdClCa/y(37))*(y(39)-ecl);
I_ClCa = I_ClCa_junc+I_ClCa_sl;
I_Clbk = GClB*(y(39)-ecl);

GClCFTR=0;%4.9e-3*ISO;     % [mS/uF]
I_ClCFTR = GClCFTR*(y(39)-ecl);

%% I_Ca: L-type Calcium Current
dss = 1/(1+exp(-(y(39)+3*ISO+9)/6)); %in Maleckar v1/2=-9 S=6 (mV); Courtemanche v1/2=-9 S=5.8 (mV)
taud = 1*dss*(1-exp(-(y(39)+3*ISO+9)/6))/(0.035*(y(39)+3*ISO+9)); 
fss = 1/(1+exp((y(39)+3*ISO+30)/7))+0.2/(1+exp((50-y(39)-3*ISO)/20)); % in Maleckar v1/2=-27.4 S=7.1 (mV); Courtemanche v1/2=-28 S=6.9 (mV)
tauf = 1/(0.0197*exp( -(0.0337*(y(39)+3*ISO+25))^2 )+0.02);
ydot(4) = (dss-y(4))/taud;
ydot(5) = (fss-y(5))/tauf;
ydot(6) = 1.7*y(36)*(1-y(6))-1*11.9e-3*y(6); % fCa_junc   koff!!!!!!!!
ydot(7) = 1.7*y(37)*(1-y(7))-1*11.9e-3*y(7); % fCa_sl
fcaCaMSL= 0.1/(1+(0.01/y(37)));
fcaCaj= 0.1/(1+(0.01/y(36)));
fcaCaMSL=0;
fcaCaj= 0;
ibarca_j = pCa*4*(y(39)*Frdy*FoRT) * (0.341*y(36)*exp(2*y(39)*FoRT)-0.341*Cao) /(exp(2*y(39)*FoRT)-1);
ibarca_sl = pCa*4*(y(39)*Frdy*FoRT) * (0.341*y(37)*exp(2*y(39)*FoRT)-0.341*Cao) /(exp(2*y(39)*FoRT)-1);
ibark = pK*(y(39)*Frdy*FoRT)*(0.75*y(35)*exp(y(39)*FoRT)-0.75*Ko) /(exp(y(39)*FoRT)-1);
ibarna_j = pNa*(y(39)*Frdy*FoRT) *(0.75*y(32)*exp(y(39)*FoRT)-0.75*Nao)  /(exp(y(39)*FoRT)-1);
ibarna_sl = pNa*(y(39)*Frdy*FoRT) *(0.75*y(33)*exp(y(39)*FoRT)-0.75*Nao)  /(exp(y(39)*FoRT)-1);
I_Ca_junc = (Fjunc_CaL*ibarca_j*y(4)*y(5)*((1-y(6))+fcaCaj)*Q10CaL^Qpow)*0.45;
I_Ca_sl = (Fsl_CaL*ibarca_sl*y(4)*y(5)*((1-y(7))+fcaCaMSL)*Q10CaL^Qpow)*0.45;
I_Ca = I_Ca_junc+I_Ca_sl;
I_CaK = (ibark*y(4)*y(5)*(Fjunc_CaL*(fcaCaj+(1-y(6)))+Fsl_CaL*(fcaCaMSL+(1-y(7))))*Q10CaL^Qpow)*0.45;
I_CaNa_junc = (Fjunc_CaL*ibarna_j*y(4)*y(5)*((1-y(6))+fcaCaj)*Q10CaL^Qpow)*0.45;
I_CaNa_sl = (Fsl_CaL*ibarna_sl*y(4)*y(5)*((1-y(7))+fcaCaMSL)*Q10CaL^Qpow)*.45;
I_CaNa = I_CaNa_junc+I_CaNa_sl;
I_Catot = I_Ca+I_CaK+I_CaNa;

% I_ncx: Na/Ca Exchanger flux
Ka_junc = 1/(1+(Kdact/y(36))^2);
Ka_sl = 1/(1+(Kdact/y(37))^2);
s1_junc = exp(nu*y(39)*FoRT)*y(32)^3*Cao;
s1_sl = exp(nu*y(39)*FoRT)*y(33)^3*Cao;
s2_junc = exp((nu-1)*y(39)*FoRT)*Nao^3*y(36);
s3_junc = KmCai*Nao^3*(1+(y(32)/KmNai)^3) + KmNao^3*y(36)*(1+y(36)/KmCai)+KmCao*y(32)^3+y(32)^3*Cao+Nao^3*y(36);
s2_sl = exp((nu-1)*y(39)*FoRT)*Nao^3*y(37);
s3_sl = KmCai*Nao^3*(1+(y(33)/KmNai)^3) + KmNao^3*y(37)*(1+y(37)/KmCai)+KmCao*y(33)^3+y(33)^3*Cao+Nao^3*y(37);


I_ncx_junc = Fjunc*IbarNCX*Q10NCX^Qpow*Ka_junc*(s1_junc-s2_junc)/s3_junc/(1+ksat*exp((nu-1)*y(39)*FoRT));
I_ncx_sl = Fsl*IbarNCX*Q10NCX^Qpow*Ka_sl*(s1_sl-s2_sl)/s3_sl/(1+ksat*exp((nu-1)*y(39)*FoRT));
I_ncx = I_ncx_junc+I_ncx_sl;

% I_pca: Sarcolemmal Ca Pump Current
I_pca_junc = Fjunc*Q10SLCaP^Qpow*IbarSLCaP*y(36)^1.6/(KmPCa^1.6+y(36)^1.6);
I_pca_sl = Fsl*Q10SLCaP^Qpow*IbarSLCaP*y(37)^1.6/(KmPCa^1.6+y(37)^1.6);
I_pca = I_pca_junc+I_pca_sl;

% I_cabk: Ca Background Current
I_cabk_junc = Fjunc*GCaB*(y(39)-eca_junc);
I_cabk_sl = Fsl*GCaB*(y(39)-eca_sl);
I_cabk = I_cabk_junc+I_cabk_sl;

%% SR fluxes: Calcium Release, SR Ca pump, SR Ca leak
MaxSR = 15; MinSR = 1;
kCaSR = MaxSR - (MaxSR-MinSR)/(1+(ec50SR/y(31))^2.5);
koSRCa = (1)*koCa/kCaSR;%
kiSRCa = kiCa*kCaSR;
RI = 1-y(14)-y(15)-y(16);
ydot(14) = (kim*RI-kiSRCa*y(36)*y(14))-(koSRCa*y(36)^2*y(14)-kom*y(15));   % R
ydot(15) = (koSRCa*y(36)^2*y(14)-kom*y(15))-(kiSRCa*y(36)*y(15)-kim*y(16));% O
ydot(16) = (kiSRCa*y(36)*y(15)-kim*y(16))-(kom*y(16)-koSRCa*y(36)^2*RI);   % I
J_SRCarel = ks*y(15)*(y(31)-y(36));          % [mM/ms]

J_serca = 1.0*Q10SRCaP^Qpow*Vmax_SRCaP*((y(38)/Kmf)^hillSRCaP-(y(31)/Kmr)^hillSRCaP)...
    /(1+(y(38)/Kmf)^hillSRCaP+(y(31)/Kmr)^hillSRCaP);
J_SRleak = (1)*(1.0+0.25*AF)*5.348e-6*(y(31)-y(36));           %   [mM/ms]


%% Sodium and Calcium Buffering
ydot(17) = kon_na*y(32)*(Bmax_Naj-y(17))-koff_na*y(17);        % NaBj      [mM/ms]
ydot(18) = kon_na*y(33)*(Bmax_Nasl-y(18))-koff_na*y(18);       % NaBsl     [mM/ms]

% Cytosolic Ca Buffers
ydot(19) = kon_tncl*y(38)*(Bmax_TnClow-y(19))-koff_tncl*y(19);            % TnCL      [mM/ms]
ydot(20) = kon_tnchca*y(38)*(Bmax_TnChigh-y(20)-y(21))-koff_tnchca*y(20); % TnCHc     [mM/ms]
ydot(21) = kon_tnchmg*Mgi*(Bmax_TnChigh-y(20)-y(21))-koff_tnchmg*y(21);   % TnCHm     [mM/ms]
ydot(22) = kon_cam*y(38)*(Bmax_CaM-y(22))-koff_cam*y(22);                 % CaM       [mM/ms]
ydot(23) = kon_myoca*y(38)*(Bmax_myosin-y(23)-y(24))-koff_myoca*y(23);    % Myosin_ca [mM/ms]
ydot(24) = kon_myomg*Mgi*(Bmax_myosin-y(23)-y(24))-koff_myomg*y(24);      % Myosin_mg [mM/ms]
ydot(25) = kon_sr*y(38)*(Bmax_SR-y(25))-koff_sr*y(25);                    % SRB       [mM/ms]
J_CaB_cytosol = sum(ydot(19:25));

% Junctional and SL Ca Buffers
ydot(26) = kon_sll*y(36)*(Bmax_SLlowj-y(26))-koff_sll*y(26);       % SLLj      [mM/ms]
ydot(27) = kon_sll*y(37)*(Bmax_SLlowsl-y(27))-koff_sll*y(27);      % SLLsl     [mM/ms]
ydot(28) = kon_slh*y(36)*(Bmax_SLhighj-y(28))-koff_slh*y(28);      % SLHj      [mM/ms]
ydot(29) = kon_slh*y(37)*(Bmax_SLhighsl-y(29))-koff_slh*y(29);     % SLHsl     [mM/ms]
J_CaB_junction = ydot(26)+ydot(28);
J_CaB_sl = ydot(27)+ydot(29);

%% Ion concentrations
% SR Ca Concentrations
ydot(30) = kon_csqn*y(31)*(Bmax_Csqn-y(30))-koff_csqn*y(30);       % Csqn      [mM/ms]
ydot(31) = J_serca-(J_SRleak*Vmyo/Vsr+J_SRCarel)-ydot(30);         % Ca_sr     [mM/ms] %Ratio 3 leak current
% ydot(31)=0;

% Sodium Concentrations
I_Na_tot_junc = I_Na_junc+I_nabk_junc+3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc+I_NaL_junc;   % [uA/uF]
I_Na_tot_sl = I_Na_sl+I_nabk_sl+3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl+I_NaL_sl;   % [uA/uF]
I_Na_tot_sl2 = 3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl;   % [uA/uF]
I_Na_tot_junc2 = 3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc;   % [uA/uF]

ydot(32) = -I_Na_tot_junc*Cmem/(Vjunc*Frdy)+J_na_juncsl/Vjunc*(y(33)-y(32))-ydot(17);
ydot(33) = -I_Na_tot_sl*Cmem/(Vsl*Frdy)+J_na_juncsl/Vsl*(y(32)-y(33))...
   +J_na_slmyo/Vsl*(y(34)-y(33))-ydot(18);
%FluxNaSL=ydot(33);
% ydot(32) = 0;
% ydot(33) = 0;
ydot(34) = J_na_slmyo/Vmyo*(y(33)-y(34));             % [mM/msec] 
% ydot(34)=0;

% Potassium Concentration
I_K_tot = I_to+I_kr+I_ks+I_ki-2*I_nak+I_CaK+I_kp+I_kur;     % [uA/uF] %SVP: added IKur
% ydot(35) = 0; %-I_K_tot*Cmem/(Vmyo*Frdy);           % [mM/msec]
ydot(35) =0; % -I_K_tot*Cmem/(Vmyo*Frdy);

% Calcium Concentrations
I_Ca_tot_junc = I_Ca_junc+I_cabk_junc+I_pca_junc-2*I_ncx_junc;                   % [uA/uF]
I_Ca_tot_sl = I_Ca_sl+I_cabk_sl+I_pca_sl-2*I_ncx_sl;            % [uA/uF]
ydot(36) = -I_Ca_tot_junc*Cmem/(Vjunc*2*Frdy)+J_ca_juncsl/Vjunc*(y(37)-y(36))...
    -J_CaB_junction+(J_SRCarel)*Vsr/Vjunc+J_SRleak*Vmyo/Vjunc;  % Ca_j
ydot(37) = -I_Ca_tot_sl*Cmem/(Vsl*2*Frdy)+J_ca_juncsl/Vsl*(y(36)-y(37))...
    + J_ca_slmyo/Vsl*(y(38)-y(37))-J_CaB_sl;   % Ca_sl
% ydot(36)=0;
% ydot(37)=0;
% ydot(38) = -J_serca*Vsr/Vmyo-J_CaB_cytosol;%+J_ca_slmyo/Vmyo*(y(37)-y(38));    % [mM/msec]
ydot(38) = -J_serca*Vsr/Vmyo-J_CaB_cytosol +J_ca_slmyo/Vmyo*(y(37)-y(38));
% ydot(38)=0;

 

%% Membrane Potential
%%
I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;          % [uA/uF]
I_Cl_tot = I_ClCa+I_Clbk++I_ClCFTR;                        % [uA/uF]
I_Ca_tot = I_Ca_tot_junc+I_Ca_tot_sl;
I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot;

% ydot(39) = -(I_tot-I_app);

ydot(39) = -(I_tot+I_app) ;

deriv = ydot' ;

return






