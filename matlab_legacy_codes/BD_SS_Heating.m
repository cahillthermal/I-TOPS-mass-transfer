function [ T_ss ] = BD_SS_Heating(k,f,kth,Cv,t,eta,nheat,ntemp,X_tempL,A_pump,r_pump )
%   Detailed explanation goes here
%   Used to calculate the steady-state heating
[G_ntemp,TempBplus,TempBminus,TempAplus,TempAminus]=BD_Spatial_frequency_TEMP_V10(k,f,kth,Cv,t,eta,nheat,ntemp,X_tempL);

T_signal=zeros(1,length(f));
for i=1:1:length(f)
T1=2*pi*A_pump.*G_ntemp(:,i).*exp(-pi^2*k.^2*r_pump^2).*k; % This is to do an average over r
T_signal(i)=trapz(k,T1); % this is the average temperature over the beam size
end;

T_ss = real(T_signal);
end

