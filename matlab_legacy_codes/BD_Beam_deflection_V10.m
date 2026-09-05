function [Signal_all] = BD_Beam_deflection_V10(k,f,kth,Cv,t,eta,n_refrac,dndT,aCET,v,Y,Vl,Vt,nheat,Namada,F,W0,Gain2,radii,A,r0,V0_multimeter,dndT_mass,Dm,aCET_mass)
% This function calculated Beam Deflections caused by the heat field and
% water concentration field 
% Z0 from the added phase of fresnel reflection by using multilayer optics
%(including effects from refractive index change in the polymer film, 
% expansion of the polymer film and interface effects)
% Z1 from the thermal expansion of thin films (Al transducer and thin film below)
% Z2 from the elastic deformation of the substrate due to lateral stresses
% in the thin films above
% Z3 from the thermal expansion of the substrate
% Z4 from the thermo-optic effect in the water vapor 

%% Parameters used in the calculation
D=kth./Cv;                       % Thermal diffusivity
omega=2*pi*f;                    % Angular frequency for heating 
%% Calulation of the temperature field 'G_nheat' (Tp/P in the paper) 
Nlayers=length(kth);             % # of layers
Nfreq=length(f);                 % # of frequency f
Nint=length(k);                  % # of k vector k, spatial frequency
kvect=k*ones(1,Nfreq);           % change k vector to matrix, with number of frequency [k,f]
kvect2=kvect.^2;                 % matrix [k,f]
kterm2=4*pi^2*kvect2;            % matrix [k,f]
ii=sqrt(-1); 
[G_nheat,TempBplus_nheat,TempBminus_nheat,TempAplus_nheat,TempAminus_nheat] = BD_Spatial_frequency_TEMP_V10(k,f,kth,Cv,t,eta,nheat,nheat,0);% matrix [k,f]


%% Calculation of the equivalent surface displacements due to Z0 to Z4
%  all [k,f] matrix
deta_phi0=zeros(Nint,Nfreq);  % Actual Phase
deta_phi1=zeros(Nint,Nfreq); 
deta_phi2=zeros(Nint,Nfreq); 
deta_phi3=zeros(Nint,Nfreq); 
deta_phi4=zeros(Nint,Nfreq); 
deta_Z0=zeros(Nint,Nfreq);    % Actual Thickness Change
deta_Z1=zeros(Nint,Nfreq); 
deta_Z2=zeros(Nint,Nfreq); 
deta_Z3=zeros(Nint,Nfreq); 
deta_Z4=zeros(Nint,Nfreq); 
deta_Z0m=zeros(Nint,Nfreq);   % Measurable Thickness Change (when comes out of materials, light change angle by n/1 times)
deta_Z1m=zeros(Nint,Nfreq); 
deta_Z2m=zeros(Nint,Nfreq); 
deta_Z3m=zeros(Nint,Nfreq); 
deta_Z4m=zeros(Nint,Nfreq); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Special Note: 
%A: The sign of deta_phi and deta_Z, decreasing the optical path
%is positive, increasing the optical path is negative
%B: Pk,Sk,Zk are calculated separated and times together in the convolution 
%Here Zk (deta_phi,deta_Z,deta_Zm)does not contains the heating power of
%the pump beam Pk, which were added afterwards
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%----(0)Fresnel reflection(Al and PMMA surface, dn/dT and aphaCET inside PMMA by heat and mass)----%
[ dR,dphi] = BD_Multilayer_reflection_heat_mass_V10(k,f,kth,Cv,t,eta,n_refrac,v,dndT,dndT_mass,aCET,aCET_mass,Dm,Namada,nheat);
deta_phi0 = -dphi;  % I assume that all top layers have the same temperature/increase the optical path so negative nominally
deta_Z0 = -dphi*Namada/(4*pi*n_refrac(1));
deta_Z0m = -dphi*Namada/(4*pi*n_refrac(1))*n_refrac(1);

%----(1)Thermal expansion of thin films(nheating layer only)-----%
if nheat<Nlayers
for i=nheat:1:Nlayers-1
q2=ones(Nint,1)*(ii*omega./D(i));
un=sqrt(4*pi^2*eta(i)*kvect2+q2); 
[Gi,TempBplusi,TempBminusi,TempAplusi,TempAminusi]=BD_Spatial_frequency_TEMP_V10(k,f,kth,Cv,t,eta,nheat,i,0);

deta_phi1 = deta_phi1+4*pi*n_refrac(1)/Namada*(1+v(i))/(1-v(i))*aCET(i)./un.*((exp(un*t(i))-1).*TempBplusi-(exp(-un*t(i))-1).*TempBminusi);
deta_Z1 = deta_Z1+(1+v(i))/(1-v(i))*aCET(i)./un.*((exp(un*t(i))-1).*TempBplusi-(exp(-un*t(i))-1).*TempBminusi);
deta_Z1m = deta_Z1m+n_refrac(1)*(1+v(i))/(1-v(i))*aCET(i)./un.*((exp(un*t(i))-1).*TempBplusi-(exp(-un*t(i))-1).*TempBminusi);
end
%The Al contribute positive signal due to expansion(norminally);measurable
%expansion of PMMA needs to time its refractive index;Other substrates time
%the refractive index of n(1), thinking from the actual effects on the
%optical path
end % This ends the calculation of (1)layer expansion

%----(2)Thermal expansion of substrate due to stresses in the thin films-----%
gT=zeros(Nint,Nfreq);    
for i=2:1:Nlayers-1 % stress due to thermal expansion
q2=ones(Nint,1)*(ii*omega./D(i));
un=sqrt(4*pi^2*eta(i)*kvect2+q2); 
[Gi,TempBplusi,TempBminusi,TempAplusi,TempAminusi]=BD_Spatial_frequency_TEMP_V10(k,f,kth,Cv,t,eta,nheat,i,0);
gT=gT+Y(i)/(1-v(i))*aCET(i)./un.*((exp(un*t(i))-1).*TempBplusi-(exp(-un*t(i))-1).*TempBminusi);
end

%%%%added stress due to the mass diffusion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q2mass=ones(Nint,1)*(ii*omega./Dm);                                       %
unmass=sqrt(4*pi^2*eta(1)*kvect2+q2mass);                                 %
l=t(nheat-1);                                                             %
gTM=gT+Y(2)/(1-v(2))*aCET_mass.*G_nheat./unmass.*sinh(unmass*l)./cosh(unmass*l);
%%%%added stress due to the mass diffusion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


kvect=k*ones(1,Nfreq);           % change k vector to matrix, with number of frequency [k,f]
kvect2=kvect.^2;                 % matrix [k,f]
kterm2=4*pi^2*kvect2;            % matrix [k,f]
W2=ones(Nint,1)*omega.^2;
apha=sqrt(4*pi^2*eta(Nlayers)*kvect2-W2/Vl(Nlayers)^2);
beta=sqrt(4*pi^2*eta(Nlayers)*kvect2-W2/Vt(Nlayers)^2);
deta_Z2_part=kterm2.*(beta.^2+kterm2-2*apha.*beta)./((beta.^2+4*pi^2*kvect2).^2-4*apha.*beta.*4*pi^2.*kvect2); %This term negative
%----add this low frequency limits to avoid the instability----
for i=1:1:Nint
    for j=1:1:Nfreq
    if W2(i,j)*10000<Vt(Nlayers)^2*kterm2(i,j)
    deta_Z2_part(i,j)=1/2*Vt(Nlayers)^2/(Vt(Nlayers)^2-Vl(Nlayers)^2);% This is an approximation to the original solution(negative)
    end
    end
end
deta_Z2=-deta_Z2_part*2*(1+v(Nlayers))/Y(Nlayers).*gTM;% Expansion of the substrate upwards so contribute positively
deta_Z2m=deta_Z2*n_refrac(1);
deta_phi2=deta_Z2*(4*pi*n_refrac(1))/Namada;


%---------(3)Thermal expansion of the substrate---------%
q2=ones(Nint,1)*(ii*omega./D(Nlayers));
un=sqrt(4*pi^2*eta(Nlayers)*kvect2+q2); 
[GN,TempBplusN,TempBminusN,TempAplusN,TempAminusN]=BD_Spatial_frequency_TEMP_V10(k,f,kth,Cv,t,eta,nheat,Nlayers,0);
W2=ones(Nint,1)*omega.^2;
apha=sqrt(4*pi^2*eta(Nlayers)*kvect2-W2/Vl(Nlayers)^2);
beta=sqrt(4*pi^2*eta(Nlayers)*kvect2-W2/Vt(Nlayers)^2);
xi=sqrt(apha.^2+q2);
deta_Z3_part=(beta.^4-16*pi^4*kvect.^4)./((beta.^2+4*pi^2*kvect.^2).^2-16*apha.*beta*pi^2.*kvect.^2); %This term is positive

%----add this to avoid the instability----
for i=1:1:Nint
    for j=1:1:Nfreq
    if W2(i,j)*10000<Vt(Nlayers)^2*kterm2(i,j)
    deta_Z3_part(i,j)=Vl(Nlayers)^2/(Vl(Nlayers)^2-Vt(Nlayers)^2);% This is an approximation to the original solution
    end
    end
end

deta_Z3=deta_Z3_part*(1+v(Nlayers))/(1-v(Nlayers))*aCET(Nlayers).*(un./q2).*(1-apha./xi).*TempBminusN; %(1-apha./xi) is positive, so overall this term contribute positively

deta_Z3m=deta_Z3*n_refrac(1); % Substrate expansion is positive
deta_phi3=deta_Z3*(4*pi*n_refrac(1))/Namada;



%-(4)Thermo-optical effect in air(due to dn/dT integration layer one)----------%
q2=ones(Nint,1)*(ii*omega./D(1));
un=sqrt(4*pi^2*eta(1)*kvect2+q2); 
[G2,TempBplus2,TempBminus2,TempAplus2,TempAminus2]=BD_Spatial_frequency_TEMP_V10(k,f,kth,Cv,t,eta,nheat,2,0);
deta_phi4 = -4*pi/Namada*dndT(1)*TempAplus2./un; % For layer 1, norminally negative since it increases the optical path norminally
deta_Z4 = deta_phi4*Namada/(4*pi*n_refrac(1));
deta_Z4m = deta_Z4*n_refrac(1);



%---------------Total Values----------------%
deta_phitot = deta_phi0+deta_phi1+deta_phi2+deta_phi3+deta_phi4; % Actual phase, Z4 is negnigible
deta_Ztot = deta_Z0+deta_Z1+deta_Z2+deta_Z3+deta_Z4;             % Actual displacement
deta_Ztotm = deta_Z0m+deta_Z1m+deta_Z2m+deta_Z3m+deta_Z4m;       % Measurable displacement




%% Calculation of Beam Deflection Angles base on Z0 to Z5
Pk=A*exp(-pi^2*k.^2*W0^2/2);    % [k,1] Vector, Hankel of the pump beam,including the profile and magnitude A 
Sk=1.00*exp(-pi^2*k.^2*W0^2/2); % [k,1] The old sk with pump and probe beam overlap

%%%Fishers beam offset simplified version%%%
%Sk=BD_TDTR_GetSHankel(k,r0,W0,1);% [k,1]  %
%%%Fishers beam offset simplified version%%%

Def_ave=zeros(1,length(omega));
for i=1:1:length(omega)   
Def=8*pi^2*Pk.*deta_Ztotm(:,i).*Sk.*besselj(1,2*pi*k*r0).*k.^2; % 8pi^2 is the D0 in Xuan's thesis
Def_ave(i)=trapz(k,Def); % this is the average beam deflection angle in radii 
end;

%----------------Add the perfect model of beam offset------------------%
perfect=0;
if perfect==1
   % first calculate the Sk by integrating Sr
   r=(0.001*W0):(0.001*W0):(5*W0);   %   [1,Nr] vector
   PreSr=2/(pi*W0^2)*exp(-2*(r.^2+r0^2)/W0^2);
   Nrmax = 80; % here only add 40 terms
   PolySr=1+zeros(1,length(r));
       for n=1:1:Nrmax
           PolySr=PolySr+1/factorial(n)^2*(4*r0^2*r.^2/W0^4).^n;
       end
   
   Sr=PreSr.*PolySr;            %  [1,Nr] vector
   
   sk=zeros(length(k),1);
   for j=1:1:length(k)
   Sr_integ=Sr.*besselj(1,2*pi*k(j)*r0).*r; % [1,Nr] vector
   sk(j)=trapz(r,Sr_integ);     %   [Nk,1] vector
   end
   
   
   Def_ave=zeros(1,length(omega));
   for i=1:1:length(omega)
   Def=16*pi^3*Pk.*deta_Ztotm(:,i).*sk.*k.^2; 
   Def_ave(i)=trapz(k,Def);     %   [1,Nf] vector
   end
   
end
%-------------------Perfect model of beam offset ends------------------%

%% Convert to actual signals measured by lock-in
Def_Prefactor=F*V0_multimeter*Gain2/(0.65*radii*2^0.5)*1e6; % This is to convert to the real signal(uV) using 3mW gain for the PSD

Signal_all=Def_ave*Def_Prefactor;     
end %end of the function
  

