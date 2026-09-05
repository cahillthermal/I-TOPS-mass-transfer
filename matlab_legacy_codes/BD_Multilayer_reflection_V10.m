function [ Reflec ] = BD_Multilayer_reflection_V10( t,n_refrac,theta,Namada )
%%% This is for the calculation of the reflectance and transmitance of a TM
%%% incident power asssume as 1, and output power reflectance
%%% operate in the time and real space domain
theta(1)=theta/180*pi;     % Change unit from theta to rad, incident angle

%%%------Derivative Parameter (heat and mass share)-------------%%%
for i=1:1:length(n_refrac)
theta(i)=asin(sin(theta(1))*n_refrac(1)/n_refrac(i));
end

C0=3e8;                                  % m/s, speed of light
OmegaL=2*pi*C0/Namada;                   % Hz, light frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Calulation of the trasfer matrix 'G' 
ii=sqrt(-1);                     % Define imagination number i
Nlayers=length(t);               % # of layers

Fplus=ones(1,1);                 % matrix of [1,1] 1 start from the last layer (bottom layer)
Fminus=zeros(1,1);               % matrix of [1,1] 0 start from the last layer (bottom layer)_no reflection from the outlet
n_refrac_Nlayersplus1=1;         % assume the air layer below the last layer has refractive index of 1
theta_Nlayersplus1=asin(sin(theta(1))*n_refrac(1)/n_refrac_Nlayersplus1); 
gammam=n_refrac_Nlayersplus1/cos(theta_Nlayersplus1);  % value of the imaginary air layer below the last layer
% Every time gets the E-field of the top surface of the layer m, 
% but the iteration starts from the layer after the last layer N+1, 
% ends at the second layer(m=3 to obtain the value of the layer 2)
% finally reformulate the top surface of layer 2 to bottom of the layer 1

    for m=(Nlayers+1):-1:3
        vmminus=ii*OmegaL*n_refrac(m-1)*cos(theta(m-1))/C0;
        gammamminus=n_refrac(m-1)/cos(theta(m-1));
        AA=gammamminus+gammam;
        BB=gammamminus-gammam;
        temp1=AA.*Fplus+BB.*Fminus;
        temp2=BB.*Fplus+AA.*Fminus;
        expterm=exp(vmminus*t(m-1));
        Fplus=(0.5/(gammamminus.*expterm))*temp1; 
        Fminus=0.5/(gammamminus)*expterm*temp2;
        % These next 3 lines fix a numerical stability issue if one of the
        % layers is very thick or resistive;
        L_penetration=Namada/(2*pi*imag(n_refrac(m-1))*cos(theta(m-1))); %penetration depth
        penetration_logic=logical(t(m-1)/L_penetration>10);        %if pentration is smaller than layer...set to semi-inf
        Fplus(penetration_logic) = 1;
        Fminus(penetration_logic) = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vm=vmminus;
        gammam=gammamminus;
    end

% Convert to the bottom of the layer 1
vmminus=ii*OmegaL*n_refrac(1)*cos(theta(1))/C0;
gammamminus=n_refrac(1)/cos(theta(1));
AA=gammamminus+gammam;
BB=gammamminus-gammam;
temp1=AA.*Fplus+BB.*Fminus;
temp2=BB.*Fplus+AA.*Fminus;
Cplus=0.5/gammamminus*temp1; 
Cminus=0.5/gammamminus*temp2;    
    
Eix=Cplus;   %incident wave
Erx=Cminus;  %reflected wave
Etx=1;       %transmitted wave

%% Calculation of the reflectance and transmitance
Reflec=abs(Erx/Eix)^2;
Trasmi=abs(Etx/Eix)^2;

end







