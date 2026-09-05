%Computes spatial frequency domain temperature response to 
%a unit heating source
%%% Bidirectional model,total Nlayers, heating layer is at the 'nheat' top
%%% surface;uplayer is from the left and always put an infinite layer as
%%% the first layer(both most left and most right boundary condition assume
%%% infinite thickness);in k ans f space, upward and downward thermal
%%% resistance are connected in parallel


%Definitions
%f: excitation frequency (Hz), ROW vector
%lambda: vector of thermal conductivities, 
%lambda(1)=top surface,(W/m-K)
%C: vector of volumetric specific heat (J/m3-K)
%t: thicknesses of each layer (layer N will NOT be used, semiinfinite)
%G: temperature at X_tempL of ntemp
%TempBplus: B(ntemp,t)plus       TempBminus: B(ntemp,t)minus
%TempAplus: A(ntemp-1,b)plus     TempAplus:  A(ntemp-1,b)minus
%G(ntemp)=TempBplus+TempBminus=TempAplus+TempAminus
function [G,TempBplus,TempBminus,TempAplus,TempAminus] = BD_Spatial_frequency_TEMP_V10(k,f,kth,Cv,t,eta,nheat,ntemp,X_tempL)
D=kth./Cv;                       % Thermal diffusivity
omega=2*pi*f;                    % Angular frequency for heating 
ii=sqrt(-1);                     % Define imagination number i
Nlayers=length(kth);             % # of layers
Nfreq=length(f);                 % # of frequency f
Nint=length(k);                  % # of k vector k, spatial frequency
%k is a COLUMN vector (actually a matrix that changes down the rows)
%f is a ROW vector
%D=kth./Cv; %already defined previously  %alpha=lambda./C; %here alpha is D---diffusivity
%omega=2*pi*f; %already defined previously
kvect=k*ones(1,Nfreq);           % change k vector to matrix, with number of frequency [k,f]
kvect2=kvect.^2;                 % matrix [k,f]
kterm2=4*pi^2*kvect2;            % matrix [k,f]

%%%-----------------------apply transfer matrix for up layers-----------------------%%%
q2=ones(Nint,1)*(ii*omega./D(1));
un=sqrt(4*pi^2*eta(1)*kvect2+q2);
gamman=kth(1)*un;
Aplus=ones(Nint,Nfreq);      % Aplus,Aminus represent the values for the bottom surface of the layer 
Aminus=zeros(Nint,Nfreq);
if nheat == 1
   gammanplus = NaN;
   Aplus = NaN;
   Aminus = NaN;
   fprintf('Program does not work for n_heat = 1 because the top surface of the first layer can not be the heating source! Use "n_heat>=2!!')
elseif nheat==2  % which means bottom of the top surface layer is heating
   %fprintf('nheat=2 ') 
   Aplus=ones(Nint,Nfreq);   % set Aplus value as 1(bottom of the top layer)
   Aminus=zeros(Nint,Nfreq); % set Aminus value as 0(bottom of the top layer)
else   
    for n=1:1:nheat-2 % n cycle down to nheat-2, got Aplus,Aminus values for the bottom of the 'nheat-1' layer  
        q2=ones(Nint,1)*(ii*omega./D(n+1));
        unplus=sqrt(eta(n+1)*kterm2+q2);
        gammanplus=kth(n+1)*unplus;
        CC=gammanplus+gamman;
        DD=gammanplus-gamman;
        temp1=CC.*Aplus+DD.*Aminus;
        temp2=DD.*Aplus+CC.*Aminus;
        expterm=exp(unplus*t(n+1));
        Aplus=0.5./(gammanplus).*expterm.*temp1;         %matrix [k,f]
        Aminus=(0.5./(gammanplus.*expterm)).*temp2;      %matrix [k,f]
        
        % To prevent NaN if pentration is smaller than layer...set to semi-inf
        % Added by Jungwoo
        penetration_logic=logical(t(n+1)*abs(unplus)>100); 
        Aplus(penetration_logic)=1;
        Aminus(penetration_logic)=0;
        %%%%%%%%%%%%%%%%%%%%
        %un=unplus;
        gamman=gammanplus;
    end
end
gamman_nheat_minus1=gamman; %put the gamman value of the layer above the heating layer


%%%----------------------------apply transfer matrix for down layers----------------------%%%
q2=ones(Nint,1)*(ii*omega./D(Nlayers)); % change frequency to matrix, with number of k vector [k,f], start with values of the last layer
un=sqrt(4*pi^2*eta(Nlayers)*kvect2+q2); % value of the last layer
gamman=kth(Nlayers)*un;                 % value of the last layer
% the three lines are for calculating gamman(Nlayers)

Bplus=zeros(Nint,Nfreq);                % matrix of [k,f] values=0 start from the last layer (bottom layer)
Bminus=ones(Nint,Nfreq);                % matrix of [k,f] values=1 start from the last layer (bottom layer)
                                        % Bplus,Bminus represent the values for the top surface of the layer 

if nheat==Nlayers        % which means top of the bottom surface layer is heating
   %fprintf('ntemp=Nlayers ')
   Bplus=zeros(Nint,Nfreq);
   Bminus=ones(Nint,Nfreq);
else
    for n=Nlayers:-1:(nheat+1) % n cycle up to nheat+1, got Bplus,Bminus values for the top of the 'nheat' layer  
        q2=ones(Nint,1)*(ii*omega./D(n-1));
        unminus=sqrt(eta(n-1)*kterm2+q2);
        gammanminus=kth(n-1)*unminus;
        AA=gammanminus+gamman;
        BB=gammanminus-gamman;
        temp1=AA.*Bplus+BB.*Bminus;
        temp2=BB.*Bplus+AA.*Bminus;
        expterm=exp(unminus*t(n-1));
        Bplus=(0.5./(gammanminus.*expterm)).*temp1; % used divided so the expterm changes sign
        Bminus=0.5./(gammanminus).*expterm.*temp2;
        % These next 3 lines fix a numerical stability issue if one of the
        % layers is very thick or resistive;
        penetration_logic=logical(t(n-1)*abs(unminus)>100);  %if pentration is smaller than layer...set to semi-inf
        Bplus(penetration_logic)=0;
        Bminus(penetration_logic)=1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %un=unminus;
        gamman=gammanminus;
    end
end
gamman_nheat=gamman; %put the gamman value of the heating layer        
Gup=(Aplus+Aminus)./(Aminus-Aplus)./gamman_nheat_minus1;
Gdown=(Bplus+Bminus)./(Bminus-Bplus)./gamman_nheat;

%G=(Bplus+Bminus)./(Bminus-Bplus)./gamman; 
%The layer G(k) is a matrix [k,f];old version for one directional heat transport 

G_nheat=1./(-1./Gup+1./Gdown);            
%Old version for bi-directional heat flow, compute the temperature at the heating surface

%% Calculation of A1plus(first layer bottom surface) and BNminus(last layer top surface)
A1plus = Gup.*Gdown./(Gup-Gdown)./(Aplus+Aminus);  % Temp of the bottom of first layer, note A1minus is 0
BNminus = Gup.*Gdown./(Gup-Gdown)./(Bplus+Bminus); % Temp of the top of last layer, note BNplus is 0
%need to verify

%% Get Aplus(or Bplus) and Aminus(or Bminus) for the temperature sensing layer 
if ntemp == 1
   fprintf('Program does not work for n_heat or n_temp = 1! Use "n_heat and n_temp>=2!!')
   
elseif ntemp == nheat
TempAplus=Aplus.*A1plus;     %set the T+ for the bottom of layer nheat-1
TempAminus=Aminus.*A1plus;   %set the T- for the bottom of layer nheat-1
TempBplus=Bplus.*BNminus;    %set the T+ for the top of layer nheat
TempBminus=Bminus.*BNminus;  %set the T- for the top of layer nheat


elseif ntemp == 2  
   %fprintf('ntemp=2 ')
   TempAplus=A1plus;             % set the Tplus for the bottom of n_temp-1 (layer 1) as A1plus
   TempAminus=zeros(Nint,Nfreq); % set the Tminus for the bottom of n_temp-1(layer 1) as 0
   q2=ones(Nint,1)*(ii*omega./D(1));
   un=sqrt(4*pi^2*eta(1)*kvect2+q2);
   gamman1=kth(1)*un;
   q2=ones(Nint,1)*(ii*omega./D(2));
   un=sqrt(4*pi^2*eta(2)*kvect2+q2);
   gamman2=kth(2)*un;
   AA=gamman2+gamman1;
   BB=gamman2-gamman1;
   temp1=AA.*TempAplus+BB.*TempAminus;
   temp2=BB.*TempAplus+AA.*TempAminus;
   TempBplus=0.5./gamman2.*temp1;           % set the Tplus for the top of ntemp(layer2)
   TempBminus=0.5./gamman2.*temp2;          % set the Tminus for the top of ntem(layer2)
   
elseif ntemp == Nlayers
   %fprintf('ntemp=Nlayers ')
   TempBplus=zeros(Nint,Nfreq);             % set the Tplus for the top of ntemp(layerN)as 0
   TempBminus=BNminus;                      % set the Tminus for the top of ntem(layerN)as BNminus
   q2=ones(Nint,1)*(ii*omega./D(Nlayers));
   un=sqrt(4*pi^2*eta(Nlayers)*kvect2+q2);
   gammanN=kth(Nlayers)*un;
   q2=ones(Nint,1)*(ii*omega./D(Nlayers-1));
   un=sqrt(4*pi^2*eta(Nlayers-1)*kvect2+q2);
   gammanNminus=kth(Nlayers-1)*un;
   AA=gammanNminus+gammanN;
   BB=gammanNminus-gammanN;
   temp1=AA.*TempBplus+BB.*TempBminus;
   temp2=BB.*TempBplus+AA.*TempBminus;
   TempAplus=0.5./gammanNminus.*temp1;      % set the Tplus for the bottom of ntemp(layerN-1)
   TempAminus=0.5./gammanNminus.*temp2;     % set the Tminus for the bottom of ntemp(layerN-1)   
   
elseif ntemp < nheat
%fprintf('ntemp<nheat ')
ATplus=ones(Nint,Nfreq);       
ATminus=zeros(Nint,Nfreq);
q2=ones(Nint,1)*(ii*omega./D(1));
un=sqrt(4*pi^2*eta(1)*kvect2+q2);
gamman=kth(1)*un;
    for n = 1:1:(ntemp-2)
        q2=ones(Nint,1)*(ii*omega./D(n+1));
        unplus=sqrt(eta(n+1)*kterm2+q2);
        gammanplus=kth(n+1)*unplus;
        AA=gammanplus+gamman;
        BB=gammanplus-gamman;
        temp1=AA.*ATplus+BB.*ATminus;
        temp2=BB.*ATplus+AA.*ATminus;
        expterm=exp(unplus*t(n+1));
        ATplus=0.5./(gammanplus).*expterm.*temp1;
        ATminus=(0.5./(gammanplus.*expterm)).*temp2;
        
        % To prevent NaN if pentration is smaller than layer...set to semi-inf
        % Added by Jungwoo
        penetration_logic=logical(t(n+1)*abs(unplus)>100); 
        ATplus(penetration_logic)=1;
        ATminus(penetration_logic)=0;
        %%%%%%%%%%%%%%%%%%%%
        %un=unplus;
        gamman=gammanplus;
    end        
 TempAplus = ATplus.*A1plus;   %set the T+ for the bottom of layer ntemp-1
 TempAminus = ATminus.*A1plus; %set the T- for the bottom of layer ntemp-1   
 q2=ones(Nint,1)*(ii*omega./D(ntemp));
 unplus=sqrt(eta(ntemp)*kterm2+q2);
 gammanplus=kth(ntemp)*unplus;
 AA=gammanplus+gamman;
 BB=gammanplus-gamman;
 temp1=AA.*TempAplus+BB.*TempAminus;
 temp2=BB.*TempAplus+AA.*TempAminus;
 TempBplus=0.5./gammanplus.*temp1;     %set the T+ for the top of layer ntemp
 TempBminus=0.5./gammanplus.*temp2;    %set the T- for the top of layer ntemp

 

elseif ntemp > nheat    
%fprintf('ntemp>=nheat ')
BTplus=zeros(Nint,Nfreq);
BTminus=ones(Nint,Nfreq);
q2=ones(Nint,1)*(ii*omega./D(Nlayers)); 
un=sqrt(4*pi^2*eta(Nlayers)*kvect2+q2); 
gamman=kth(Nlayers)*un;                 
    for n=Nlayers:-1:(ntemp+1) % n cycle up to nheat+1, got Bplus,Bminus values for the top of the 'nheat' layer  
        q2=ones(Nint,1)*(ii*omega./D(n-1));
        unminus=sqrt(eta(n-1)*kterm2+q2);
        gammanminus=kth(n-1)*unminus;
        AA=gammanminus+gamman;
        BB=gammanminus-gamman;
        temp1=AA.*BTplus+BB.*BTminus;
        temp2=BB.*BTplus+AA.*BTminus;
        expterm=exp(unminus*t(n-1));
        BTplus=(0.5./(gammanminus.*expterm)).*temp1; % used divided so the expterm changes sign
        BTminus=0.5./(gammanminus).*expterm.*temp2;
        % These next 3 lines fix a numerical stability issue if one of the
        % layers is very thick or resistive;
        penetration_logic=logical(t(n-1)*abs(unminus)>100);  %if pentration is smaller than layer...set to semi-inf
        BTplus(penetration_logic)=0;
        BTminus(penetration_logic)=1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %un=unminus;
        gamman=gammanminus;
    end 
 TempBplus = BTplus.*BNminus;   %set the T+ for the top of layer ntemp
 TempBminus = BTminus.*BNminus; %set the T- for the top of layer ntemp
 q2=ones(Nint,1)*(ii*omega./D(ntemp-1));
 unminus=sqrt(eta(ntemp-1)*kterm2+q2);
 gammanminus=kth(ntemp-1)*unminus;
 AA=gammanminus+gamman;
 BB=gammanminus-gamman;
 temp1=AA.*TempBplus+BB.*TempBminus;
 temp2=BB.*TempBplus+AA.*TempBminus;
 TempAplus=0.5./gammanminus.*temp1;     %set the T+ for the top of layer ntemp
 TempAminus=0.5./gammanminus.*temp2;    %set the T- for the top of layer ntemp
 
end

if X_tempL<=0 % sensing position above the top surface(x=0)of the ntemp layer 
q2 = ones(Nint,1)*(ii*omega./D(ntemp-1));
unminus = sqrt(4*pi^2*eta(ntemp-1)*kvect2+q2);
G = TempAplus.*exp(unminus*X_tempL) + TempAminus.*exp(-unminus*X_tempL); %The layer G(k)
elseif X_tempL>0  % sensing position below the top surface(x=0) of the ntemp layer  
q2 = ones(Nint,1)*(ii*omega./D(ntemp));
un = sqrt(4*pi^2*eta(ntemp)*kvect2+q2);
G = TempBplus.*exp(un*X_tempL) + TempBminus.*exp(-un*X_tempL); %The layer G(k) matrix[k,f]
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%For testing purpose only%%%%%%%% 
%q2 = ones(Nint,1)*(ii*omega./D(ntemp-1));
%unminus = sqrt(4*pi^2*eta(ntemp-1)*kvect2+q2);
%gammanminus=kth(ntemp-1)*unminus;
%q2 = ones(Nint,1)*(ii*omega./D(ntemp));
%un = sqrt(4*pi^2*eta(ntemp)*kvect2+q2);
%gamman=kth(ntemp)*un;
%G=(gamman.*(TempBminus-TempBplus)-gammanminus.*(TempAminus-TempAplus));
%should be 0 if nheat~=ntemp, should be 1 if nheat=ntemp

end

