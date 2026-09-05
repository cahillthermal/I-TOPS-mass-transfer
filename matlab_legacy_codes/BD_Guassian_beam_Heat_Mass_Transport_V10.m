%%% FDPBD modeling and fitting on a layered structure, including: 
%%% (1)Processing data (phase correction)
%%% (2)Calculating temperature excursion
%%% (3)Calculating beam deflection
%%% (4)Fitting the beam deflection data
%%% (5)Saving the results
%%% The first layer and the last layer are always thermally thick
%%% The heating layer is the top surface of 'nheat'
%%% The output temperature is the top surface of 'ntemp'

clc;clear;
close all;tic;
%% Program Choices and Input Files
%(1)Import Data for Comparision with simulation? 0 for no, 1 for yes
Importdata=1;
filename_Exp='091917_2_R23_Scan_Vacuum(0.04psi)_PMMA(114nm)_Al_SiO2_5X_pump3mW_probe6mW_TC3s_12dB_Sens200uV_PSD920mV_Wait15s_dewell15s_12um';

%(2)Phase correction? 0 for no, 1 for yes; %Amplitude correction? 0 for no, 1
%for yes (only do amplitude correction if phase correction is allowed)
Phasecorrection=1; Amplitudecorrection=1; Plot_phasecorrection=0;
filename_Cali='091917_6_Cali_R23_Before';

%(3)Data value moving downwards?(Fix the issue of labview recording for data collected before 10/18/2017) 0 for no, 1 for yes  
Fix_labview=1;

%(4)Calculate the temperature rise (K) from multiple layers due to heat transport? 0 for no, 1 for yes
Calculate_temperature_rise=0;     

%(5)Calculate the probe beam deflection? 0 for no, 1 for yes (if nheat=2, then do not calculate mass transport)
Simulate_Beam_Deflection=1;    

%(6)Fit the deflection data? 0 for no, 1 for yes (only valid if importdata==1, automatically updating fitting results)
Fit_data=0;  

%(7)Plot heat- and mass-induced probe beam deflection, respectively? Plot overall beam deflection? 
Plot_deflection_heat=0; Plot_deflection_mass=0; Plot_all=1;

%(10)Save the data/imported/simulated/fitted?
Save_data=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Input Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Input Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% All Parameters that need to change

%%%-------------------Material Parameter(Heat)-----------------------%%%
%%%(Air kth=0.026 W/m/K, Cv,p=0.0012 J/cm^3/K, n=1, dn/dT=-9.1e-7/K)
%%%(Water vapor kth=0.019 W/m/K, Cv,p= 40e-6 J/cm^3/K, n=1, dn/dT=-2.1e-8/K)
%%%(PMMA kth=0.19 W/m/K, Cv=1.60 J/cm^3/K,n=1.49,dn/dT=-0.78e-4/K,aCET=33e-6,v=0.37,Y=2.5Gpa)
%%%(Al kth=160 W/m/K, Cv=2.42 J/cm^3/K,nAl=1.96+7.10*1i, dn/dT=calculated, aCET=23e-6/K, v=0.33, Y=69Gpa)
%%%(Fused silica kth=1.32 W/m/K, Cv=1.63 J/cm^3/K, n=1.4536, dn/dT=1e-5, aCET=0.54e-6/K,v=0.17,Y=73GPa)


%%% Note: Most left is infinite air or water vapor
%%%       Most right infinite substrate
%%%       Heating layer is the upper surface of the nheat layer
kth=[0.019 0.19 160 1.32];                      % W/m/k, thermal conductivity of substrate---right layer is the substrate
Cv=[40e-6 1.60 2.42 1.63]*1e6;                  % J/m^3/K, heat capacity of the substrate---right layer is the substrate
t=[1e9 100 80 1e9]*1e-9;                        % m, thickness of the sample
n_refrac=[1 1.49 (1.95+7.10*1i) 1.45];          % 1, refractive index of material

dRdT=1e-4;                % 1/K, reflectance change with T, used for getting dndT_Al (obtained from FDTR)
dMdT=1/2*dRdT/0.87;       % 1/K, E-field amplitude change with T, used for getting dndT_Al
dfydT_Vac=[3e-5];         % rad/K, used for getting dndT_Al (obtained by fitting Al/SiO2 sample in vaccum)
nAl=1.96+7.10*1i;         % Basic values used to derive the dfyH2O_dT or the other layer
dndT_Al=0.5*(nAl^2-1)*(dMdT+1i*dfydT_Vac); % Shall use this to fit later on

dndT=[-2.1e-8 -0.78e-4 dndT_Al NaN];       % /K, thermo-optic coefficient of materials;
aCET=[NaN 33 23 0.54]*1e-6;                % /K, thermal expansion coefficient
%dndT=[-0e-8 -0e-4 0 NaN];                 % /K, thermo-optic coefficient of materials;
%aCET=[NaN 0 0 0]*1e-6;                    % /K, thermal expansion coefficient

v=[NaN 0.37 0.33 0.17];                    % 1, possion's ratio
Y=[NaN 2.5 69 73]*1e9;                     % Pa, Young's modulus
Vl=[NaN NaN NaN 5970];                     % m/s, longitudinal speed of sound of amorppous SiO2
Vt=[NaN NaN NaN 3770];                     % m/s, transverse speed of sound of amorphous SiO2
eta=ones(1,numel(kth));                    % Isotropic layers, eta=kx/ky;

%%%-------------------Material Parameter(Mass)-----------------------%%%
%Mass transport fitting need 
Dm=1e-12;                  % m2/s Water diffusion coefficient in polymer 
dndT_mass=0.5e-4;          % /K, This is the dn/dC*dC/dT       (dn/dC negative dC/dT negative)
aCET_mass=-0.25e-4;        % /K, this is the (dL/L/dC)*(dC/dT) (dL/L/dC positive, dC/dT negative)
%dndT_mass=0e-4;            % /K, This is the dn/dC*dC/dT
%aCET_mass=-0e-4;          % /K, this is the (dL/L/dC)*(dC/dT)

%Xguess=[dndT_mass*1e4,Dm*1e12]; % Initial guess for the fitting solution; dndT_mass[1e-4],Dm[1e-12] easy for fitting
Xguess=[dndT_Al];


%%%-------------------Experimental Parameters-----------------------%%%
f=logspace(log10(1),log10(100e3),100);       % Heating frequency used for simulation and range of fitting
nheat=3;                                     % The layer that generates heat,designed to be the upper interface of the transducer(nheat)
ntemp=3;                                     % The layer where the top surface temperature is calculated
 
%typically nheat>2(the top surface of the first layer can not be reached)
%if nheat=2,then topsurface of the second layer is heating.
%if nheat=length(t),then the topsurface of the last layer is heating(first layer and last layer always has infinite thickness), however, 
%nheat can not be 1, because the top surface of the 1st layer(which is infinite) can not be the heating layer

P=3.0e-3;                   % W, incident pump laser power measured at the back of the objective lens
Namada=783e-9;              % m, laser wavelength
W0=11.3e-6*(Namada/783e-9); % m,size of the laser beam 2.7um (20X) 10.6um(5X) 1.0um (50X)
r0=3e-6;                   % m, this is the offset of the pump and probe beam if doing deflection measurement 
F=40e-3;                    % m, focal length of the objective lens F=200e-3/magnification;

V0_multimeter=830e-3;       % V, usually shown as mV(use for converting to the actual signal)
G2=1;                       % second stage gain of the PSD (1,or 3, or 10, use for converting to the actual signal)

k=((0.01/W0):(0.01/W0):(5/W0))'; % Spatial frequency vector[Note needs be a collumn vector with dimension [k,1] ]

%%%------Derivative Parameter (heat and mass share)-------------%%%
R=BD_Multilayer_reflection_V10( t,n_refrac,0,Namada );    % Overall reflectivity  
Trans=0.90*0.92;            % Transmission of the 5X objective lens 0.9, 20X only 0.7, the chamber glass 0.92
A=(4/pi*P)*(1-R)*1.08*Trans;% When using Al as the absorption layer
radii=Namada*F/pi/W0;       % Free space laser diameter, usually measured in front of the PSD
radii=1.55e-3/2;            % Free space laser diameter measured by beam profiler

%%%------Derivative Parameter (do not change)-------------%%%
D=kth./Cv;                  % Thermal diffusivity
omega=2*pi*f;               % Angular frequency for heating 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------- CALCULATE STEADY STATE HEATING--------------
A_SS = (2*A+2*A)*pi/4;      % Incident light intensity,Watts, A is the intensity for pump, probe intensity set as 2A; All modes contribute to the ss-heating so factor 4/pi cancelled
X_tempL=-0e-9;              % The depth of the temperature sensing, negative means above the ntemp, positive means below the ntemp
                            % if negative, abs(X_tempL)<= t(ntemp-1); if positive,abs(X_tempL)<=t(ntemp)
dT_SS = BD_SS_Heating(k,0.01,kth,Cv,t,eta,nheat,ntemp,X_tempL,A_SS,W0) %max steady state temperature rise

%----------------CALCULATE TEMPERATURE RISE----------------
if Calculate_temperature_rise==1  
X_tempL=-0e-9;  % The depth of the temperature sensing, negative means above the ntemp, positive means below the ntemp
[Vin_T,Vout_T] = BD_Temperature_rise(k,f,kth,Cv,t,eta,nheat,ntemp,X_tempL,A,W0);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import experimental data, perform coherent correction and phase correction
%%%-----------------Import Experimental Data-----------------------%%%
if Importdata==1
filename1 = sprintf('%s.txt',filename_Exp);
file_exp  = fopen(filename1,'r');
ExpData   = fscanf(file_exp,'%12e %12e %12e %12e %12e %12e %12e',[7,5000]);
f_exp = ExpData(1,:)';
t_exp = ExpData(2,:)';
Vin_exp = ExpData(3,:)';
Vout_exp= ExpData(4,:)';
V0_exp=ExpData(6,:)';
Complex_exp=Vin_exp+1i*Vout_exp;
R_exp=abs(Complex_exp);
Phase_exp=angle(Complex_exp)/pi*180; % in degrees
end 

%%%---------------Calibration Phase Correction----------------------%%%
if Phasecorrection==1
%-----------------load Phase Calibration data--------------%
filename3 = sprintf('%s.txt',filename_Cali);
file_cali  = fopen(filename3,'r');
CaliData   = fscanf(file_cali,'%12e %12e %12e %12e %12e %12e %12e',[7,5000]);
f_cali = CaliData(1,:)';
Vin_cali = CaliData(3,:)';
Vout_cali = CaliData(4,:)';
Complex_cali = Vin_cali+1i*Vout_cali; % Vout=Rsin(theta) Vin=Rcos(theta)
R_cali=abs(Complex_cali);
Phase_cali=angle(Complex_cali)/pi*180;
%-----------------Plot Phase Calibration Curve--------------%
if Plot_phasecorrection==1
figure(7);
semilogx(f_cali,R_cali,'ok');
xlabel('frequency (Hz)','FontSize',15)
ylabel('Amplitude(uV)','FontSize',15)
title('Calibration Amplitude','FontSize',15)
set(gca,'FontSize',15)
figure(8);
semilogx(f_cali,Phase_cali,'or');
xlabel('frequency (Hz)','FontSize',15)
ylabel('Phase(degree)','FontSize',15)
title('Calibration Phase','FontSize',15)
set(gca,'FontSize',15)
end 
%------------------Phase Calibration process------------------%
Theta=interp1(f_cali,Phase_cali,f_exp,'linear'); % Intercept if frequency for cali and exp is not the same;
Vin_exp1=Vin_exp.*cos(Theta/180*pi)+Vout_exp.*sin(Theta/180*pi);
Vout_exp1=-Vin_exp.*sin(Theta/180*pi)+Vout_exp.*cos(Theta/180*pi);

if Amplitudecorrection==1
R_cali_normal=R_cali/mean(R_cali);               %  Normalize the signal to the average signals
R_cali_normal_int=interp1(f_cali,R_cali_normal,f_exp,'linear');
Vin_exp2=Vin_exp1./R_cali_normal_int;
Vout_exp2=Vout_exp1./R_cali_normal_int;
end 

Vin_exp=Vin_exp2;
Vout_exp=Vout_exp2;
Complex_exp=Vin_exp+1i*Vout_exp;
R_exp=abs(Complex_exp);
Phase_exp=angle(Complex_exp)/pi*180;
end


%%%-----Extract the range of data for fitting(reverse the sequence of the data)------%%%
if Importdata==1
f_min=min(f); 
f_max=max(f);
[f_data,R_data] = BD_extract_interior(f_exp',R_exp',f_min,f_max); % need to change f_data and other data from collumn vector to row vector
[f_data,Phase_data] = BD_extract_interior(f_exp',Phase_exp',f_min,f_max);
[f_data,Vin_data] = BD_extract_interior(f_exp',Vin_exp',f_min,f_max);
[f_data,Vout_data] = BD_extract_interior(f_exp',Vout_exp',f_min,f_max);

if Fix_labview==1                   % Fix the labview issue (only needed for data before 10/18/2017)
points=length(f_data);              % add to fix the labview of frequency recording issues
f_data=f_data(2:(points));          % add to fix the labview of frequency recording issues
R_data=R_data(1:(points-1));        % add to fix the labview of frequency recording issues
Phase_data=Phase_data(1:(points-1));% add to fix the labview of frequency recording issues
Vin_data=Vin_data(1:(points-1));    % add to fix the labview of frequency recording issues
Vout_data=Vout_data(1:(points-1));  % add to fix the labview of frequency recording issues
end                                 % Fix the labview ends


V_data=Vin_data+1i*Vout_data;
else
f_data=f';
V_data=zeros(length(f),1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation for the beam deflection

%-------------Calculate Beam Deflection Caused by Heat and Mass Flow--------------
if Simulate_Beam_Deflection==1
[Simu_signal] = BD_Beam_deflection_V10(k,f,kth,Cv,t,eta,n_refrac,dndT,aCET,v,Y,Vl,Vt,nheat,Namada,F,W0,G2,radii,A,r0,V0_multimeter,dndT_mass,Dm,aCET_mass);
Simu_R=abs(Simu_signal); Simu_Phase=angle(Simu_signal)/pi*180; Simu_Vin=real(Simu_signal);Simu_Vout=imag(Simu_signal);

Mass_signal=Simu_signal*0;
%[Mass_signal] = BD_Beam_deflection_V10(k,f,kth,Cv,t,eta,n_refrac,zeros(1,length(dndT)),zeros(1,length(aCET)),v,Y,Vl,Vt,nheat,Namada,F,W0,G2,radii,A,r0,V0_multimeter,dndT_mass,Dm,aCET_mass);
Simu_R_Mass=abs(Mass_signal); Simu_Phase_Mass=angle(Mass_signal)/pi*180; Simu_Vin_Mass=real(Mass_signal);Simu_Vout_Mass=imag(Mass_signal);

Heat_signal=Simu_signal*0;
%[Heat_signal] = BD_Beam_deflection_V10(k,f,kth,Cv,t,eta,n_refrac,dndT,aCET,v,Y,Vl,Vt,nheat,Namada,F,W0,G2,radii,A,r0,V0_multimeter,0,Dm,0);
Simu_R_Heat=abs(Heat_signal); Simu_Phase_Heat=angle(Heat_signal)/pi*180; Simu_Vin_Heat=real(Heat_signal);Simu_Vout_Heat=imag(Heat_signal);

%---------Plot simulation results------%
%'Importdata', 'plot_deflection_heat', 'plot_deflection_mass', 'plot_all'
% are all '0' or '1'
BD_plot_all(Importdata,f_data,V_data,f,Simu_signal,Heat_signal,Mass_signal,Plot_deflection_heat,Plot_deflection_mass,1,2,3,4,Plot_all);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Fitting
%-------------------Fit Experimental data and plot results------------------
if Importdata==1
if Fit_data==1
 
Xsol=fminsearch(@(X) BD_Beam_deflection_fit_V9(X,Vin_data,Vout_data,f_data,k,kth,Cv,t,eta,n_refrac,dndT,aCET,v,Y,Vl,Vt,nheat,Namada,F,W0,A,r0,V0_multimeter,dndT_mass,Dm),Xguess,optimset('TolFun',1e-2,'TolX',1e-2))
fprintf('Data fit completed\n')    
[Z,Heat_fit,Mass_fit] = BD_Beam_deflection_fit_V9(Xsol,Vin_data,Vout_data,f_data,k,kth,Cv,t,eta,n_refrac,dndT,aCET,v,Y,Vl,Vt,nheat,Namada,F,W0,A,r0,V0_multimeter,dndT_mass,Dm);
Z
R_fit=abs(Heat_fit+Mass_fit);
Phase_fit=angle(Heat_fit+Mass_fit)/pi*180;
Vin_fit=real(Heat_fit+Mass_fit);
Vout_fit=imag(Heat_fit+Mass_fit);

%-----------------------Plot Fitting Results------------------------------
BD_plot_all(importdata,f_data,V_data,f_data,Heat_fit,Mass_fit,plot_deflection_heat,plot_deflection_mass,21,22,23,24);
    
end % for Fit_data=1
end % for importdata=1

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Saving all data
if Save_data==1
% If fit, then save fit&exp. data, if not fit, then save simulation&exp.data
  if importdata==1     
     if Fit_data==1
        CalcRatioFileName = sprintf('%s%s.txt','Fit_Exp_',filename_Exp);
        fid2 = fopen(CalcRatioFileName,'w'); 
        for j=1:length(f_data)
        fprintf(fid2,'\n %6f %6f %6f %6f %6f',f_data(j),Vin_data(j),Vout_data(j),Vin_fit(j),Vout_fit(j)); 
        end
        fclose(fid2);
     elseif Simulate_Beam_Deflection==1 
        CalcRatioFileName = sprintf('%s%s.txt','Sim_Exp_',filename_Exp);
        fid2 = fopen(CalcRatioFileName,'w'); 
        for j=1:length(f_data)
        fprintf(fid2,'\n %6f %6f %6f %6f %6f',f_data(j),Vin_data(j),Vout_data(j),Simu_Vin(j),Simu_Vout(j)); 
        end
        fclose(fid2);    
     end % for Fit_data=1
  elseif Simulate_Beam_Deflection==1
     CalcRatioFileName = sprintf('%s%s.txt','Sim_',filename_Exp);
     fid2 = fopen(CalcRatioFileName,'w'); 
     for j=1:length(f_data)
     fprintf(fid2,'\n %6f %6f %6f',f(j),Simu_Vin(j),Simu_Vout(j)); 
     end
     fclose(fid2);        
  end  % for importdata=1

end;   % for Save_data=1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc;  

