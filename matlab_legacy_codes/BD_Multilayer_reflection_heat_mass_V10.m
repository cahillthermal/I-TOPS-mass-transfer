function [ dR_heat_mass,dphi_heat_mass] = BD_Multilayer_reflection_heat_mass_V10(k,f,kth,Cv,t,eta,n_refrac,v,dndT,dndT_mass,aCET,aCET_mass,Dm,Namada,nheat)
%%% Calculation of the enhanced magnitude and added phase of the reflectance by multilayer optics 
%%% Include effects from the temperature field and concentration field
%%% k(Nk,1), f(1,Nf), Gnheat(Nk,Nf), output are also matrix [Nk,Nf]
%%% good for the case air/Al/SiO2 and vapor/polymer/Al/SiO2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%-All Derivative Parameters-%%%%%%
Nf=length(f);                    % number of the heating frequency
OmegaH=2*pi*f;                   % Hz, Angular frequency for heating [1,Nf]
Nk=length(k);                    % number of the spatial frequency
kvect=k*ones(1,Nf);              % change k vector to matrix, with number of frequency [Nk,Nf] 
kvect2=kvect.^2;                 % matrix [k,f]
ii=sqrt(-1); 
C0=3e8;                          % m/s, speed of light
OmegaL=2*pi*C0/Namada;           % Hz, light optical frequency [1,1]
[G_nheat,TempBplus_nheat,TempBminus_nheat,TempAplus_nheat,TempAminus_nheat] = BD_Spatial_frequency_TEMP_V10(k,f,kth,Cv,t,eta,nheat,nheat,0);% matrix [k,f]

%%%%%%-Descretize the PMMA layer into multiple layers-%%%%%%
l=t(nheat-1);                    % m, thickness of the polymer layer 
Omegamax=2*pi*max(f);            % Hz, maximum heating frequency
l_diff=(Dm/Omegamax)^0.5;        % m, mass diffusion length
N=round(l/(l_diff/5))+1;         % descretize the PMMA layer into N new layers, at least having one layer

Nnew=length(t)-1+N;              % Total number of the new layers

%%%%%%-Define thickness and material parameters in the new Layer System-%%%%%%
t_new=zeros(1,Nnew);
n_refrac_new=zeros(1,Nnew);
v_new=zeros(1,Nnew);
aCET_new=zeros(1,Nnew);
dndT_new=zeros(1,Nnew);

if nheat==2 % Then the layer system does not change
%%%%%%-Add the special case of air/Al/SiO2 (where no polymer exist)-%%%%%%
t_new=t;
n_refrac_new=n_refrac;
v_new=v;
aCET_new=aCET; 
dndT_new=dndT;
else
%%%%%%-when the polymer exist,descretize it-%%%%%%    
for i=1:1:Nnew
if i<(nheat-1)                      % above the polymer layer
t_new(i)=t(i);
n_refrac_new(i)=n_refrac(i);
v_new(i)=v(i);
aCET_new(i)=aCET(i);
dndT_new(i)=dndT(i);
elseif (nheat-1)<=i&&i<=(nheat-2+N) % total N layers, replacing the (nheat-1) layer
t_new(i)=t(nheat-1)/N;
n_refrac_new(i)=n_refrac(nheat-1);
v_new(i)=v(nheat-1);
aCET_new(i)=aCET(nheat-1);
dndT_new(i)=dndT(nheat-1);
elseif i>=nheat-1+N                 % below the polymer layer
t_new(i)=t(i-N+1);
n_refrac_new(i)=n_refrac(i-N+1);
aCET_new(i)=aCET(i-N+1);
dndT_new(i)=dndT(i-N+1);
end
end

   
end % end for if nheat==2 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calulation of the original reflectance  
Nlayers=length(t_new);           % # of new layers
Fplus=ones(Nk,Nf);          % matrix of [Nk,Nf] start from the last layer (bottom layer)
Fminus=zeros(Nk,Nf);        % matrix of [Nk,Nf] start from the last layer (bottom layer)_no reflection from the outlet
n_refrac_Nlayersplus1=1;         % assume the air layer below the last layer has refractive index of 1
q2=ones(Nk,Nf)*(n_refrac_Nlayersplus1*OmegaL/C0)^2;
vm=sqrt(4*pi^2*kvect2-q2);
gammam=vm;
% Every time gets the light at the top surface of layer m, but the iteration 
% starts from the top of the layer N+1; when goes to the top of the layer 2 
% and convert to the bottom of the layer 1 

    for m=(Nlayers+1):-1:3
        q2=ones(Nk,Nf)*(n_refrac_new(m-1)*OmegaL/C0)^2;
        vmminus=sqrt(4*pi^2*kvect2-q2);
        gammamminus=vmminus;
        AA=gammamminus+gammam;
        BB=gammamminus-gammam;
        temp1=AA.*Fplus+BB.*Fminus;
        temp2=BB.*Fplus+AA.*Fminus;
        expterm=exp(vmminus*t_new(m-1));
        Fplus=(0.5./(gammamminus.*expterm)).*temp1; 
        Fminus=0.5./(gammamminus).*expterm.*temp2;
        % These next 3 lines fix the issue if the layer m-1 is optically opaque;
        L_penetration=Namada/(2*pi*imag(n_refrac_new(m-1)));      %penetration depth
        penetration_logic=logical(t_new(m-1)/L_penetration>10);   %if pentration is smaller than layer...set to semi-inf
        Fplus(penetration_logic) = 1;
        Fminus(penetration_logic) = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vm=vmminus;
        gammam=gammamminus;
    end

% Convert to the bottom of the layer 1
q2=ones(Nk,Nf)*(n_refrac_new(1)*OmegaL/C0)^2;
vmminus=sqrt(4*pi^2*kvect2-q2);
gammamminus=vmminus;
AA=gammamminus+gammam;
BB=gammamminus-gammam;
temp1=AA.*Fplus+BB.*Fminus;
temp2=BB.*Fplus+AA.*Fminus;
Cplus=0.5./gammamminus.*temp1;  
Cminus=0.5./gammamminus.*temp2;    
    
Eix=Cplus;   %incident wave in layer 1
Erx=Cminus;  %reflected wave in layer 1
%---reflectance and phase---%
r=Erx./Eix;
R=abs(Erx./Eix).^2;
phi=angle(r);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calulation of the reflectance in temperature and water concentration
%% "Real" and "Imaginary" part regarding the heating frequency

%%%%%%-Material and Experiment Parameter-%%%%%%
qmass2=ones(Nk,1)*(ii*OmegaH/Dm); % note this is at the heat frequency, not the optical frequency
Namada_mass=sqrt(4*pi^2*kvect2+qmass2); % mass transport related

%%%%%%%%%%%%%%%%%%%%Calculation of the 'real' part%%%%%%%%%%%%%%%%%%%
%--------Updated thickness and material parameters 'real'--------%
tT=zeros(Nk,Nf,Nlayers);         % the new thickness at Tp        matrix arrays [k,f,Nlayers]
n_refracT=zeros(Nk,Nf,Nlayers);  % the new refractive index at Tp matrix arrays [k,f,Nlayers]

deltaT=real(G_nheat);         % get the real value of Tp


if nheat==2 % For Al/SiO2 sample, then the layer system does not change 
%%%%%%-Add the special case of air/Al/SiO2 (where no polymer exist)-%%%%%%
  for i=1:1:Nlayers   % Nlayers=Nnew only when nheat~=2
  tT(:,:,i)=t_new(i); % ignore the expansion for this case cause has no effect on optical matrix
  n_refracT(:,:,i)=n_refrac_new(i)+dndT_new(i)*deltaT;
     if isnan(dndT_new(i))==1 % if NaN exists
     n_refracT(:,:,i)=n_refrac_new(i);
     end
  end
else
%%%%%%-when the polymer exist,change n and t by deltaT deltaC-%%%%%% 
for i=1:1:Nnew
if i<(nheat-1)
tT(:,:,i)=t_new(i);
n_refracT(:,:,i)=n_refrac_new(i)+dndT_new(i)*deltaT;
   if isnan(dndT_new(i))==1 % if NaN exists
   n_refracT(:,:,i)=n_refrac_new(i);
   end

elseif (nheat-1)<=i&&i<=(nheat-2+N) % total N layers, replacing the nheat-1 layer
Z=(N-(i-(nheat-2))+1/2)*t_new(i);   % i-(nheat-2) is layer 1, corresponds to the coordinate Z=N-1+0.5; 
deltaC=real(G_nheat.*cosh(Namada_mass*Z)./cosh(Namada_mass*l)); % real part of the concentration
tT(:,:,i)=t_new(i)+aCET_new(i)*t_new(i)*(1+v_new(i))/(1-v_new(i))*deltaT+aCET_mass*t_new(i)*(1+v_new(i))/(1-v_new(i))*deltaC;
n_refracT(:,:,i)=n_refrac_new(i)+dndT_new(i)*deltaT+dndT_mass*deltaC;

elseif i>=nheat-1+N
tT(:,:,i)=t_new(i);
n_refracT(:,:,i)=n_refrac_new(i)+dndT_new(i)*deltaT;
   if isnan(dndT_new(i))==1 % if NaN exists
   n_refracT(:,:,i)=n_refrac_new(i);
   end
end % end the if loop
end % end the for loop
end % end of the if nheat==2


%--------Updated reflectance by transfer matrix 'real' --------%
Nlayers=length(t_new);           % # of new layers
Fplus=ones(Nk,Nf);          % matrix of [Nk,Nf] start from the last layer (bottom layer)
Fminus=zeros(Nk,Nf);        % matrix of [Nk,Nf] start from the last layer (bottom layer)_no reflection from the outlet
n_refrac_Nlayersplus1T=ones(Nk,Nf);% assume the air layer below the last layer has refractive index of 1
q2=(n_refrac_Nlayersplus1T*OmegaL/C0).^2;
vm=sqrt(4*pi^2*kvect2-q2);
gammam=vm;
% Every time gets the light at the top surface of layer m, but the iteration 
% starts from the top of the layer N+1; when goes to the top of the layer 2 
% and convert to the bottom of the layer 1 

    for m=(Nlayers+1):-1:3
        q2=(n_refracT(:,:,m-1)*OmegaL/C0).^2;
        vmminus=sqrt(4*pi^2*kvect2-q2);
        gammamminus=vmminus;
        AA=gammamminus+gammam;
        BB=gammamminus-gammam;
        temp1=AA.*Fplus+BB.*Fminus;
        temp2=BB.*Fplus+AA.*Fminus;
        expterm=exp(vmminus.*tT(:,:,m-1));
        Fplus=(0.5./(gammamminus.*expterm)).*temp1; 
        Fminus=0.5./(gammamminus).*expterm.*temp2;
        % These next 3 lines fix the issue if the layer m-1 is optically opaque;
        if (m-1)~=1
        L_penetration=Namada./(2*pi*imag(n_refracT(:,:,m-1)));      %penetration depth
        penetration_logic=logical(tT(:,:,m-1)./L_penetration>10);   %if pentration is smaller than layer...set to semi-inf
        Fplus(penetration_logic) = 1;
        Fminus(penetration_logic) = 0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vm=vmminus;
        gammam=gammamminus;
    end

% Convert to the bottom of the layer 1
q2=(n_refracT(1)*OmegaL/C0).^2;
vmminus=sqrt(4*pi^2*kvect2-q2);
gammamminus=vmminus;
AA=gammamminus+gammam;
BB=gammamminus-gammam;
temp1=AA.*Fplus+BB.*Fminus;
temp2=BB.*Fplus+AA.*Fminus;
Cplus=0.5./gammamminus.*temp1;  
Cminus=0.5./gammamminus.*temp2;    
    
EixT=Cplus;   %incident wave in layer 1
ErxT=Cminus;  %reflected wave in layer 1

%--------Updated reflectance and phase 'real'--------%
rT=ErxT./EixT;
RT=abs(ErxT./EixT).^2;
phiT=angle(rT);
dphi_real=phiT-phi;      %(Nk,Nf) % This is the real part of dphi
dphi_real2=imag(rT./r);  % another definition
dR_real=RT-R;            %(Nk,Nf) % This is the real part of dR



%%%%%%%%%%%%%%%%%%%%Calculation of the 'imag' part%%%%%%%%%%%%%%%%%%%
%--------Updated thickness and material parameters 'imag'--------%
tT=zeros(Nk,Nf,Nlayers);         % the new thickness at Tp        matrix arrays [k,f,Nlayers]
n_refracT=zeros(Nk,Nf,Nlayers);  % the new refractive index at Tp matrix arrays [k,f,Nlayers]

deltaT=imag(G_nheat);         % get the imaginary value of Tp


if nheat==2 % For Al/SiO2 sample, then the layer system does not change 
%%%%%%-Add the special case of air/Al/SiO2 (where no polymer exist)-%%%%%%
  for i=1:1:Nlayers   % Nlayers=Nnew only when nheat~=2
  tT(:,:,i)=t_new(i); % ignore the expansion for this case cause has no effect on optical matrix
  n_refracT(:,:,i)=n_refrac_new(i)+dndT_new(i)*deltaT;
    if isnan(dndT_new(i))==1 % if NaN exists
    n_refracT(:,:,i)=n_refrac_new(i);
    end
  end
else
%%%%%%-when the polymer exist,change n and t by deltaT deltaC-%%%%%% 
for i=1:1:Nnew
if i<(nheat-1)
tT(:,:,i)=t_new(i);
n_refracT(:,:,i)=n_refrac_new(i)+dndT_new(i)*deltaT;
   if isnan(dndT_new(i))==1 % if NaN exists
   n_refracT(:,:,i)=n_refrac_new(i);
   end

elseif (nheat-1)<=i&&i<=(nheat-2+N) % total N layers, replacing the nheat-1 layer
Z=(N-(i-(nheat-2))+1/2)*t_new(i);   % i-(nheat-2) is layer 1, corresponds to the coordinate Z=N-1+0.5; 
deltaC=imag(G_nheat.*cosh(Namada_mass*Z)./cosh(Namada_mass*l)); % real part of the concentration
tT(:,:,i)=t_new(i)+aCET_new(i)*t_new(i)*(1+v_new(i))/(1-v_new(i))*deltaT+aCET_mass*t_new(i)*(1+v_new(i))/(1-v_new(i))*deltaC;
n_refracT(:,:,i)=n_refrac_new(i)+dndT_new(i)*deltaT+dndT_mass*deltaC;

elseif i>=nheat-1+N
tT(:,:,i)=t_new(i);
n_refracT(:,:,i)=n_refrac_new(i)+dndT_new(i)*deltaT;
   if isnan(dndT_new(i))==1 % if NaN exists
   n_refracT(:,:,i)=n_refrac_new(i);
   end
end % end the if loop
end % end the for loop
end % end of the if nheat==2


%--------Updated reflectance by transfer matrix 'imag' --------%
Nlayers=length(t_new);           % # of new layers
Fplus=ones(Nk,Nf);          % matrix of [Nk,Nf] start from the last layer (bottom layer)
Fminus=zeros(Nk,Nf);        % matrix of [Nk,Nf] start from the last layer (bottom layer)_no reflection from the outlet
n_refrac_Nlayersplus1T=ones(Nk,Nf);% assume the air layer below the last layer has refractive index of 1
q2=(n_refrac_Nlayersplus1T*OmegaL/C0).^2;
vm=sqrt(4*pi^2*kvect2-q2);
gammam=vm;
% Every time gets the light at the top surface of layer m, but the iteration 
% starts from the top of the layer N+1; when goes to the top of the layer 2 
% and convert to the bottom of the layer 1 

    for m=(Nlayers+1):-1:3
        q2=(n_refracT(:,:,m-1)*OmegaL/C0).^2;
        vmminus=sqrt(4*pi^2*kvect2-q2);
        gammamminus=vmminus;
        AA=gammamminus+gammam;
        BB=gammamminus-gammam;
        temp1=AA.*Fplus+BB.*Fminus;
        temp2=BB.*Fplus+AA.*Fminus;
        expterm=exp(vmminus.*tT(:,:,m-1));
        Fplus=(0.5./(gammamminus.*expterm)).*temp1; 
        Fminus=0.5./(gammamminus).*expterm.*temp2;
        % These next 3 lines fix the issue if the layer m-1 is optically opaque;
        if (m-1)~=1
        L_penetration=Namada./(2*pi*imag(n_refracT(:,:,m-1)));      %penetration depth
        penetration_logic=logical(tT(:,:,m-1)./L_penetration>10);   %if pentration is smaller than layer...set to semi-inf
        Fplus(penetration_logic) = 1;
        Fminus(penetration_logic) = 0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vm=vmminus;
        gammam=gammamminus;
    end

% Convert to the bottom of the layer 1
q2=(n_refracT(1)*OmegaL/C0).^2;
vmminus=sqrt(4*pi^2*kvect2-q2);
gammamminus=vmminus;
AA=gammamminus+gammam;
BB=gammamminus-gammam;
temp1=AA.*Fplus+BB.*Fminus;
temp2=BB.*Fplus+AA.*Fminus;
Cplus=0.5./gammamminus.*temp1;  
Cminus=0.5./gammamminus.*temp2;    
    
EixT=Cplus;   %incident wave in layer 1
ErxT=Cminus;  %reflected wave in layer 1

%--------Updated reflectance and phase 'imag'--------%
rT=ErxT./EixT;
RT=abs(ErxT./EixT).^2;
phiT=angle(rT);
dphi_imag=phiT-phi;      %(Nk,Nf) % This is the real part of dphi
dphi_imag2=imag(rT./r);  % another definition
dR_imag=RT-R;            %(Nk,Nf) % This is the real part of dR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%The Final dphidT and dRdT(Complex Number)%%%%%%%%%%%%%%%%%
dphi_heat_mass=dphi_real+dphi_imag*ii; %(Nk,Nf)
dR_heat_mass=dR_real+dR_imag*ii;       %(Nf,Nf)
end % end of program 
