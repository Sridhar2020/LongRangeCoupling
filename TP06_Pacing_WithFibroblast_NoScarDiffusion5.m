%%%Author S.Sridhar (dharmails@gmail.com)
%%%Affl: University of Sheffield, UK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%This program generates simulates the electrical activity for
%%%2D left ventricular tissue for the human heart using the
%%%TNNP-TP06 model and MacCannell active fibroblast model.
%%%The tissue comprises of a central non-conducting scar region
%%%with uncoupled myocytes and fibroblasts attached randomly
%%%to the myocytes in the scar region.Some of these fibroblasts are also
%%%be randomly connected to myocytes surrounding the scar region.
%%%Pacing waves are generated from one edge of the domain and 
%%%the resulting interaction between the fibrotic tissue and the pacing 
%%%wave is investigated.

%%%The input tissue parameters are set to Shallow restitution curve
%%%The resting membrane potential for the fibroblast is set to V_FR = -24.5 mV
%%%Input files 
%%% a) The tissue initial conditions - Shallow_TP06_Size400_D0_1_FromT6610ms_Upto7200ms.mat UU Nt   
%%% b) The adjacency link matrix with lambda = 50, np = 30K - 
%%%    LinkMatrix_Size400_ScarRad80_and_BorderRad100_Lam50_Np30000_NewCase5.mat

%%%Output files
%%% a) PacingOutput_Period300_Gs1_Range10_NewCase3.dat 
%%% The output is the voltage at all points on the grid every 20 ms. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

tic
for i2= 1:1
    Per = 300; %Pacing period
    for j2=1:1
        %Loading link adjacency matrix for FDist = 10, np 30000 
        load LinkMatrix_Size400_ScarRad80_and_BorderRad100_Lam50_Np30000.mat adj np FDist
        %Loading initial conditions for TP06 myocytes
        load Shallow_TP06_Size400_D0_1_FromT6610ms_Upto7200ms.mat UU Nt
        

        threshold = -73.0; % myocyte Voltage threshold 
        
        %Simulation parameters
        N=400;%lattice size
        N2=N*N;
        mmm=0;
        dx=0.25;dt=0.01;
        
       %Duration of simulation
        Ntold = Nt;% Copy older Nt to Ntold
        clear Nt;
        Nt= Ntold+600000; %Duration of simulation
        %n=[N 1:N-1]; s=[2:N 1]; %Periodic BC
        n=[2 1:N-1]; s=[2:N N-1];e=n;w=s; % No flux
        
        %Myocyte and fibroblast coupling parameters
        Gs= 1.0;  %M-F Coupling strength
        MCap=185; %Myocyte capacitance
        FCap=50; %Fibroblast capacitance
        NUMFIB=8; %Number of fibroblasts in a fibroblast unit
        D = 0.1; %Diffusion constant
        MyoT = NUMFIB*(Gs/MCap);
        FibT = (Gs/FCap);

        %Setting scar boundary diffusion values
        x_cent = N/2; y_cent = N/2;
        radius = 80; % Radius of scar region
        for i=1:N
            for j=1:N
                     diffcoff(i,j)=D.*(0.5+0.5*sign((i-x_cent).^2+(j-y_cent).^2- eps - (radius.^2)));               
            end
        end
        
    	%Output file 
        fname1 = sprintf('PacingOutput_Period%d_Gs%f_Range%d_New.dat',Per,Gs,FDist);
        f2=fopen(fname1,'w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Myocyte Parameters
% Terms for Solution of Conductance and Reversal Potential
Rgas = 8314.472;      % Universal Gas Constant (J/kmol*K)
Frdy = 96485.3415;  % Faraday's Constant (C/mol)
Temp = 310.0;    % Temperature (K) 37C
RTonF = (Rgas * Temp) / Frdy;

% Cellular capacitance
CAPACITANCE = 0.185;

% External concentrations
Ko = 5.4;
KoNorm = 5.4;
Cao = 2.0;
Nao = 140.0;
Nao3 = 2744000.0;

% Intracellular volumes
Vc = 0.016404;
Vsr = 0.001094;
Vss = 0.00005468;

% Calcium buffering dynamics
Bufc = 0.2;
Kbufc = 0.001;
Bufsr = 10.0;
Kbufsr = 0.3;
Bufss = 0.4;
Kbufss = 0.00025;

% Intracellular calcium flux dynamics
Vmaxup = 0.006375;
Kup = 0.00025;
Vrel = 0.102; % 40.8;
k1bar = 0.15;
k2bar = 0.045;
k3 = 0.060;
k4 = 0.005; % 0.000015;
EC = 1.5;
maxsr = 2.5;
minsr = 1.0;
Vleak = 0.00036;
Vxfer = 0.0038;

%RestitutionVariables
% Parameters for IKr
GkrPar1 = 0.134;
GkrPar2 = 0.153;
GkrPar3 = 0.172;
GkrPar4 = 0.172;
GkrR2 = 0.134 * 1.75;
GkrR1 = 0.134 * 0.8;

% Parameters for Iks
pKNa = 0.03;
GksEpi = 0.392;
GksEndo = 0.392;
GksMcell = 0.098;
GksPar1 = 0.270;
GksPar2 = 0.392;
GksPar3 = 0.441;
GksPar4 = 0.441;
GksR2 = 0.270 * 1.75;
GksR1 = 0.270 * 0.8;	

%RestitutionVariables
% Parameters for Ik1
GK1=5.405;
% Parameters for Ito
GtoEpi=0.294;
GtoEndo=0.073;
GtoMcell=0.294;
%	Parameters for INa
GNa=14.838; % nS/PF
INa_Vshift = 0.0; 
%Parameters for IbNa
GbNa=0.00029;
%Parameters for INaK
KmK=1.0;
KmNa=40.0;
knak=2.724;
%Parameters for ICaL
GCaL=0.00003980;
GCaL_atp = 1.0; 
GCaL_pH = 1.0; %
%	Parameters for IbCa
GbCa=0.000592;
%Parameters for INaCa
naca_pH = 1.0; 
knaca=1000;
KmNai=87.5;
KmCa=1.38;
ksat=0.1;
nn=0.35;
%Parameters for IpCa
 GpCaPar1=0.0619;
 GpCaPar2=0.1238;
 GpCaPar3=0.3714;
 GpCaPar4=0.8666;
 KpCa=0.0005;
%Parameters for IpK;
 GpKPar1=0.0730;
 GpKPar2=0.0146;
 GpKPar3=0.0073;
 GpKPar4=0.00219;
 inverseVcF2=1.0/(2.0*Vc*Frdy);
 inverseVcF=1.0/(Vc*Frdy);
 inversevssF2=1.0/(2.0*Vss*Frdy);


tau_f_multiplier = 0.6; % par for shallow restitution

myo_current =zeros(1,N2) ; %myocyte ionic current term initialising

%%%%%%%%%%%%%%%%%%%% Fibroblast Dynamics%%%%%%%%%%%%%%%%%%
     %Parameters for fibroblasts
FIB_GKV = 0.155; %//0.25; // units: ns/pF
FIB_GK1 = 0.4822;% //units: ns/pF
FIB_I_NAK = 2.002;% //units pA/pF
FIB_G_BNA = 0.0095; %//units ns/pF
fib_Ek = -87.0; %//units mV
fib_B = -200.0; %//units mV
V_rev = -150.0; %//units mV
FIB_K_mK = 1.0; %//units mmol/L
FIB_K_mna = 11.0; %//units mmol/L
FIB_Ko = 5.3581; %//units mmol/L
FIB_Nai = 8.5547; %//units mmol/L
FIB_CAP = 0.050; %//nF //units pF 
FIB_Nao = 130.0;
% /* Fibroblast sodium reversal potentials */
 Ena_f = RTonF*log(FIB_Nao/FIB_Nai);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %num_states1= 20;
        num_states2 = 3; 
        UU1= ones(np,num_states2);
        ve = UU(:,1);
        U1(:,1) = -24.5;
        U1(:,2) = 0;
        U1(:,3) = 1; 
        Iion= zeros(1,N*N);
        Ifib= zeros(np,1);
        
%Time loop
for ti=Ntold+1:Nt
          
     U=UU; U1=UU1;UUold = UU(:,1);
     %Myocyte Dynamics
     for i1 = 1:N2 
        
        %Stimulating tissue to generate pacing waves
         if (mod(ti,(Per/dt)) > 0) && (mod(ti,(Per/dt))< 101) && (i1 <= 10*N)
           stimCurrent1 = -52.0;
        else
           stimCurrent1 = 0.0;
         end


     VmoRTonF = U(i1,1)/RTonF;

    %Reversal potentials
    Ena = RTonF*log(Nao/U(i1,19));
    Ek = RTonF*(log((Ko/U(i1,20))));
    Eks = RTonF*(log((Ko+pKNa*Nao)/(U(i1,20)+pKNa*U(i1,19))));
    Eca = 0.5*RTonF*(log((Cao/U(i1,18))));
      
   %%%%%%%%Ion channel dynamics%%%%%%%%%
   % Inward current iNa
    alpha_m = 1.0/(1.0+exp((-60.0-U(i1,1))/5.0));
    beta_m = 0.1/(1.0+exp((U(i1,1)+35.0)/5.0))+0.10/(1.0+exp((U(i1,1)-50.0)/200.0));
    tau_m = alpha_m * beta_m;
    m_inf = 1.0/((1.0+exp((-56.86-U(i1,1))/9.03))*(1.0+exp((-56.86-U(i1,1))/9.03)));
       if (U(i1,1) >= -40.0)
	      alpha_h = 0.0;
	      beta_h = (0.77/(0.13*(1.0 + exp(-(U(i1,1) + 10.66)/11.1))));
       else
	    alpha_h = (0.057 * exp(-(U(i1,1) + 80.0)/6.8));
	    beta_h = (2.7 * exp(0.079 * U(i1,1))+(3.1e5) * exp(0.3485 * U(i1,1)));
       end
	tau_h = 1.0/(alpha_h + beta_h);
	h_inf = 1.0/((1.0 + exp((U(i1,1) + 71.55)/7.43))*(1.0 + exp((U(i1,1) + 71.55)/7.43)));

	 
	      if (U(i1,1) >= -40.0)
                 alpha_j = 0;
	         beta_j = ((0.6 * exp((0.057) * U(i1,1)))./(1.0 + exp(-0.1 * (U(i1,1) + 32.0))));           
              else
	         alpha_j = (((-2.5428e4)*exp(0.2444*U(i1,1))-(6.948e-6)*exp(-0.04391*U(i1,1)))*(U(i1,1)+37.78)./(1.0+exp(0.311*(U(i1,1)+79.23))));
		  beta_j = (0.02424*exp(-0.01052*U(i1,1))/(1.0+exp(-0.1378*(U(i1,1)+40.14))));
       end

	tau_j = 1.0/(alpha_j + beta_j);
	j_inf = h_inf;
       	  
 
  	U(i1,2) = m_inf - ( m_inf - U(i1,2) ) * exp( -dt / tau_m );
  	U(i1,3) = h_inf - ( h_inf - U(i1,3) ) * exp( -dt / tau_h );
	U(i1,4) = j_inf - ( j_inf - U(i1,4) ) * exp( -dt / tau_j );

  	INa = GNa*U(i1,2)*U(i1,2)*U(i1,2)*U(i1,3)*U(i1,4)*(U(i1,1)-Ena);

% Currents in Ca channels 
    d_inf = 1.0/(1.0+exp((-8.0-U(i1,1))/7.5));
    ad = 1.4/(1.0+exp((-35.0-U(i1,1))/13.0))+0.25;
    bd=1.4/(1.0+exp((U(i1,1)+5.0)/5.0));
    cd1=1.0/(1.0+exp((50.0-U(i1,1))/20.0));
    tau_d = ad * bd + cd1;

    f_inf = 1.0/(1.0+exp((U(i1,1)+20.0)/7.0));
    af=1102.5*exp(-(U(i1,1)+27.0)*(U(i1,1)+27.0)/225.0);
	bf=200.0/(1.0+exp((13.0-U(i1,1))/10.0));
	cf=(180.0/(1.0+exp((U(i1,1)+30.0)/10.0)))+20.0;
   if (U(i1,1) < 0)
    tau_f = af + bf + cf;
   else
    tau_f= tau_f_multiplier*(af + bf + cf);
   end
   
  
    f2_inf = 0.67/(1.0+exp((U(i1,1)+35.0)/7.0))+0.33;
    af2=562.0*exp(-(U(i1,1)+27.0)*(U(i1,1)+27.0)/240.0);
    bf2=31.0/(1.0+exp((25.0-U(i1,1))/10.0));
    cf2=80.0/(1.0+exp((U(i1,1)+30.0)/10.0));
    tau_f2 = af2+bf2+cf2;
  

	fCass_inf = 0.6/(1.0+(U(i1,16)/0.05)*(U(i1,16)/0.05))+0.4;
	tau_fCass = 80.0/(1.0+(U(i1,16)/0.05)*(U(i1,16)/0.05))+2.0;

        U(i1,7) = d_inf - (d_inf - U(i1,7)) * exp( -dt / tau_d );
  	U(i1,8) = f_inf - (f_inf - U(i1,8)) * exp( -dt / tau_f );
  	U(i1,9) = f2_inf - (f2_inf - U(i1,9)) * exp( -dt / tau_f2 );
  	U(i1,10) = fCass_inf - (fCass_inf - U(i1,10)) * exp( -dt / tau_fCass );

  	ICaL = GCaL_pH*GCaL_atp*GCaL*U(i1,7)*U(i1,8)*U(i1,9)*U(i1,10)*4.0*(U(i1,1)-15.0)*(Frdy/RTonF)*(0.25*exp(2.0*(U(i1,1)-15.0)/RTonF)*U(i1,16)-Cao)/(exp(2.0*(U(i1,1)-15.0)/RTonF)-1.0);

% Rapidly inactivating K current 
  	
    xr1_inf = 1.0/(1.0+exp((-26.0-U(i1,1))/7.0));
    axr1 = 450.0/(1.0+exp((-45.0-U(i1,1))/10.0));
    bxr1 = 6.0/(1.0+exp((U(i1,1)+30.0)/11.5));
    tau_xr1 = axr1 * bxr1;

    xr2_inf = 1.0/(1.0+exp((U(i1,1)+88.0)/24.0));
    axr2 = 3.0/(1.0+exp((-60.0-U(i1,1))/20.0));
    bxr2 = 1.12/(1.0+exp((U(i1,1)-60.0)/20.0));
    tau_xr2 = axr2 * bxr2;
     
 
    Gkr = GkrPar1;%for shallow restitution
    U(i1,11) = xr1_inf - (xr1_inf - U(i1,11)) * exp( -dt / tau_xr1 );
    U(i1,12) = xr2_inf - (xr2_inf - U(i1,12)) * exp( -dt / tau_xr2 );
    IKr = Gkr*sqrt(Ko/5.4)*U(i1,11)*U(i1,12)*(U(i1,1)-Ek);

% Slowly inactivating K current
        xs_inf = 1.0/(1.0+exp((-5.0-U(i1,1))/14.0));
	axs = (1400.0/(sqrt(1.0+exp((5.0-U(i1,1))/6.0))));
	bxs = (1.0/(1.0+exp((U(i1,1)-35.0)/15.0)));
	tau_xs = axs * bxs + 80.0;
     	

	U(i1,13) = xs_inf - (xs_inf - U(i1,13)) * exp( -dt / tau_xs );

	Gks = GksPar1; %for shallow restitution
	IKs = Gks*U(i1,13)*U(i1,13)*(U(i1,1)-Eks);

  % Time independent K current 
  	Ak1 = 0.1/(1.0+exp(0.06*(U(i1,1)-Ek-200.0)));
	Bk1 = (3.0*exp(0.0002*(U(i1,1)-Ek+100.0))+exp(0.1*(U(i1,1)-Ek-10.0)))/(1.0+exp(-0.5*(U(i1,1)-Ek)));
	rec_iK1 = Ak1/(Ak1+Bk1);
	IK1 = GK1*rec_iK1*(U(i1,1) - Ek);

  %Plateau K current
  	GpK=GpKPar1; % for shallow restitution
	rec_ipK = 1.0/(1.0+exp((25.0-U(i1,1))/5.98));
	IpK=GpK*rec_ipK*(U(i1,1)-Ek);


 %transient outward current
    r_inf = 1.0/(1.0+exp((20.0-U(i1,1))/6.0));
	s_inf = 1.0/(1.0+exp((U(i1,1)+20.0)/5.0));
	tau_r = 9.5*exp(-(U(i1,1)+40.0)*(U(i1,1)+40.0)/1800.0)+0.8;
	tau_s = 85.0*exp(-(U(i1,1)+45.0)*(U(i1,1)+45.0)/320.0)+5.0/(1.0+exp((U(i1,1)-20.0)/5.0))+3.0;
	       	
	Gto = GtoEpi;

	U(i1,6) = s_inf - (s_inf - U(i1,6)) * exp(-dt / tau_s);
	U(i1,5) = r_inf - (r_inf - U(i1,5)) * exp(-dt / tau_r);
	Ito = Gto*U(i1,5)*U(i1,6)*(U(i1,1)-Ek);


  %Na Ca exchanger
	naca1 = knaca*(1.0/(KmNai*KmNai*KmNai+Nao3))*(1.0/(KmCa+Cao));
	naca2 = (1.0/(1.0+ksat*exp((nn-1.0)*VmoRTonF)));
	naca3 = (exp(nn*VmoRTonF)*U(i1,19)*U(i1,19)*U(i1,19)*Cao-exp((nn-1.0)*VmoRTonF)*Nao3*U(i1,18)*2.5);
	INaCa = naca_pH * naca1 * naca2 * naca3;

	
  %Background Na current
	IbNa=GbNa*(U(i1,1)-Ena);

  %iNaK
	rec_iNaK = (1.0/(1.0+0.1245*exp(-0.1*VmoRTonF)+0.0353*exp(-VmoRTonF)));
	INaK=knak*(Ko/(Ko+KmK))*(U(i1,19)/(U(i1,19)+KmNa))*rec_iNaK;

  %Plateau Ca current
 	GpCa=GpCaPar1; %for shallow restitution
  	IpCa=GpCa*U(i1,18)/(KpCa+U(i1,18));

  % Background Ca current
	IbCa=GbCa*(U(i1,1)-Eca);


   %intracellular ion concentrations 
	kCaSR = maxsr-((maxsr-minsr)/(1.0+(EC/U(i1,17))*(EC/U(i1,17))));
	k1 = k1bar/kCaSR;
	k2 = k2bar*kCaSR;
	dRR = k4 * (1.0-U(i1,14)) - k2*U(i1,16)*U(i1,14);
	U(i1,14) = U(i1,14) + dt*dRR;
	U(i1,15) = k1*U(i1,16)*U(i1,16)*U(i1,14)/(k3+k1*U(i1,16)*U(i1,16));
	Irel = Vrel*U(i1,15)*(U(i1,17)-U(i1,16));
	Ileak = Vleak*(U(i1,17)-U(i1,18));
	Iup = Vmaxup/(1.0+((Kup*Kup)/(U(i1,18)*U(i1,18))));
	Ixfer = Vxfer*(U(i1,16) - U(i1,18));

	CaCSQN = Bufsr*U(i1,17)/(U(i1,17)+Kbufsr);
	dCaSR = dt*(Iup-Irel-Ileak);
	bjsr = Bufsr-CaCSQN-dCaSR-U(i1,17)+Kbufsr;
	cjsr = Kbufsr*(CaCSQN+dCaSR+U(i1,17));
	U(i1,17) = (sqrt(bjsr*bjsr+4.0*cjsr)-bjsr)/2.0;

	CaSSBuf=Bufss*U(i1,16)/(U(i1,16)+Kbufss);
	dCaSS = dt*(-Ixfer*(Vc/Vss)+Irel*(Vsr/Vss)+(-ICaL*inversevssF2*CAPACITANCE));
	bcss = Bufss-CaSSBuf-dCaSS-U(i1,16)+Kbufss;
	ccss = Kbufss*(CaSSBuf+dCaSS+U(i1,16));
	U(i1,16) = (sqrt(bcss*bcss+4.0*ccss)-bcss)/2.0;

	CaBuf = Bufc*U(i1,18)/(U(i1,18)+Kbufc);
	dCai = dt*((-(IbCa+IpCa-2.0*INaCa)*inverseVcF2*CAPACITANCE)-(Iup-Ileak)*(Vsr/Vc)+Ixfer);
	bc = Bufc-CaBuf-dCai-U(i1,18)+Kbufc;
	cc = Kbufc*(CaBuf+dCai+U(i1,18));
	U(i1,18) = (sqrt(bc*bc+4.0*cc)-bc)/2.0;

	dNai = -(INa+IbNa+3.0*INaK+3.0*INaCa)*inverseVcF*CAPACITANCE;
	U(i1,19) = U(i1,19) +  dt*dNai;

	dKi= -(stimCurrent1+IK1+Ito+IKr+IKs-2.0*INaK+IpK)*inverseVcF*CAPACITANCE;
	U(i1,20) = U(i1,20) + dt*dKi;

        myo_current(i1)=IKr + IKs + IK1 + Ito + INa + IbNa + ICaL + IbCa + INaK + INaCa + IpCa + IpK + stimCurrent1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
     UU=U;
     clear U;
     Iion = reshape(myo_current,N,N);
     ve =  reshape(UU(:,1),N,N);  
    
 %%%%%%%%%%%%%%%%%Fibroblast dynamics%%%%%%%%%%%%%%%%%%
  for i1 = 1:np
   % /*inward rectifying potassium*/ 
    fib_alpha_K1 = 0.1/(1+exp(0.06*(U1(i1,1)-fib_Ek-200)));
    fib_beta_K1= 3*exp(0.0002*(U1(i1,1)-fib_Ek+100))+exp(0.1*(U1(i1,1)-fib_Ek-10))/(1+exp(-0.5*(U1(i1,1)-fib_Ek)));
    fib_rec_IK1= fib_alpha_K1/(fib_alpha_K1 + fib_beta_K1);
    fib_IK1 = FIB_GK1*fib_rec_IK1*(U1(i1,1)-fib_Ek);

    %/*Na-K pump in the fibroblast*/
     fib_INaK = FIB_I_NAK*(FIB_Ko/(FIB_Ko+FIB_K_mK))*(FIB_Nai.^(1.5))/((FIB_Nai.^(1.5))+(FIB_K_mna.^(1.5)))*((U1(i1,1)-V_rev)/(U1(i1,1) - fib_B));  

   %/*inward Kv current*/
     fib_r_bar = 1/(1 + exp(-(U1(i1,1) + 20.0 -30)/11));
     fib_s_bar = 1/(1 + exp((U1(i1,1) + 23.0 -30 )/7));
     tr1 = (U1(i1,1) + 20 -30)/25.9;
     tr2 = (U1(i1,1) + 23 -30)/22.7;  
     fib_tau_r = 20.3 + 138*exp(-(tr1*tr1));
     fib_tau_s = 1574 + 5268*exp(-(tr2*tr2));	    
     U1(i1,2) = fib_r_bar - (fib_r_bar-U1(i1,2))*exp(-dt/fib_tau_r);
     U1(i1,3) = fib_s_bar - (fib_s_bar-U1(i1,3))*exp(-dt/fib_tau_s); 
     fib_IKv = FIB_GKV*U1(i1,2)*U1(i1,3)*(U1(i1,1)-fib_Ek);       
   
 % /*Background sodium current*/
     fib_IbNa = FIB_G_BNA*(U1(i1,1) - Ena_f);    
     Ifib = fib_IKv + fib_IK1 + fib_INaK + fib_IbNa;

  end
     UU1=U1;clear U1;
     vp=UU1(:,1);
     spat = diffcoff.*(((0.5+ 0.5*sign(diffcoff(n,:)-eps)).*(ve(n,:)-ve))+((0.5+ 0.5*sign(diffcoff(:,n)-eps)).*(ve(:,n)-ve))+((0.5+ 0.5*sign(diffcoff(:,s)-eps)).*(ve(:,s)-ve)) + ((0.5+ 0.5*sign(diffcoff(s,:)-eps)).*(ve(s,:)-ve)));
     ve=ve+dt*(-Iion + MyoT*(reshape((adj*vp),N,N)-reshape(sum(adj,2),N,N).*ve)+ spat/dx^2);
     vp=vp+dt*(-Ifib'+FibT*((reshape(ve,1,N*N)*adj)'- (sum(adj)'.*vp)));     
     
     UU(:,1) = reshape(ve,1,N*N);
     UU1(:,1)=vp;

   %Recording voltage at all points every 20 ms.
   if mod(ti,2000)==0
    fprintf(f2,'%f\n',ve); 
   end
end
  fclose(f2);

    end
end

toc



