function [Vin_T,Vout_T] = BD_Temperature_rise(k,f,kth,Cv,t,eta,nheat,ntemp,X_tempL,A_pump,r_pump)
% Detailed explanation goes here
% Calculate the inphase and out-of-phase temperature rise
[G_ntemp,TempBplus,TempBminus,TempAplus,TempAminus]=BD_Spatial_frequency_TEMP_V10(k,f,kth,Cv,t,eta,nheat,ntemp,X_tempL);

T_signal=zeros(1,length(f));
for i=1:1:length(f)
T1=2*pi*A_pump.*G_ntemp(:,i).*exp(-pi^2*k.^2*r_pump^2).*k; % This is to do an average over r
T_signal(i)=trapz(k,T1); % this is the average temperature over the beam size
end;

Vin_T = real(T_signal);
Vout_T = imag(T_signal);
R_T = abs(T_signal);
Phase_T = angle(T_signal)/pi*180;

% Plot the calculated data
figure(51);
set(gcf,'units','normalized','position',[0.04 0.22 0.44 0.53])
semilogx(f,Vin_T,'--b');
xlabel('frequency (Hz)','FontSize',12)
ylabel('Vin(K)','FontSize',12)
set(gca,'FontSize',12)
figure(52);
set(gcf,'units','normalized','position',[0.53 0.22 0.44 0.53])
semilogx(f,Vout_T,'--r');
xlabel('frequency (Hz)','FontSize',12)
ylabel('Vout(K)','FontSize',12)
set(gca,'FontSize',12) 

%figure(53);
%set(gcf,'units','normalized','position',[0.04 0.22 0.44 0.53])
%semilogx(f,R_T,'--k');
%xlabel('frequency (Hz)','FontSize',12)
%ylabel('Amp(K)','FontSize',12)
%set(gca,'FontSize',12)
end

