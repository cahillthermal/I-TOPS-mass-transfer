function [ Z,Heat_fit,Mass_fit ] = BD_Beam_deflection_fit_V9(X,Vin_data,Vout_data,f_data,k,kth,Cv,t,eta,n_refrac,dndT,aCET,v,Y,Vl,Vt,nheat,Namada,F,W0,A,r0,V0_multimeter,dndT_mass,Dm)
%Summary of this function goes here
%Fit the Beam Deflection data Obtained
%dndT_mass=X(1)*1e-4;
%Dm=X(2)*1e-12;
dndT(2)=X(1);
[Heat_fit,Mass_fit] = BD_Beam_deflection_V9(k,f_data,kth,Cv,t,eta,n_refrac,dndT,aCET,v,Y,Vl,Vt,nheat,Namada,F,W0,A,r0,V0_multimeter,dndT_mass,Dm);
Vin_fit=real(Heat_fit+Mass_fit);
Vout_fit=imag(Heat_fit+Mass_fit);
Phase_fit=angle(Heat_fit+Mass_fit)/pi*180;
Phase_data=angle(Vin_data+1i*Vout_data)/pi*180;
%%%%%See the optimization process%%%%%
%Plot Vin
figure(23);
set(gcf,'units','normalized','position',[0.04 0.22 0.44 0.53])
semilogx(f_data,Vin_data,'ob',f_data,Vin_fit,'-b'); 
xlabel('frequency (Hz)','FontSize',12)
ylabel('Vin(\muV)','FontSize',12)
set(gca,'FontSize',12)
%Plot Vout
figure(24);
set(gcf,'units','normalized','position',[0.53 0.22 0.44 0.53])
semilogx(f_data,Vout_data,'or',f_data,Vout_fit,'-r'); 
xlabel('frequency (Hz)','FontSize',12)
ylabel('Vout(\muV)','FontSize',12)
set(gca,'FontSize',12)

pause(0.1)
X
%res=(Vin_fit-Vin_data).^2+10*(Vout_fit-Vout_data).^2;
res=(Phase_fit-Phase_data).^2;
Z=sum(res*1e6);%This uses uV as the unit
end

