function [ output_args ] = BD_plot_all(Importdata,f_data,V_data,f,All_signal,Heat_signal,Mass_signal,plot_deflection_heat,plot_deflection_mass,P_amplitude,P_phase,P_inphase,P_outofphase,plot_all)
%This plots the amplitude, phase, in-phase and out-of-phase
Simu_R=abs(All_signal);
Simu_Phase=angle(All_signal)/pi*180;
Simu_Vin=real(All_signal);
Simu_Vout=imag(All_signal);
R_data=abs(V_data);
Phase_data=angle(V_data)/pi*180;
Vin_data=real(V_data);
Vout_data=imag(V_data);

%-------------Plot Amplitude-----------%
figure(P_amplitude);set(gcf,'units','normalized','position',[0.04 0.38 0.44 0.53])
if plot_all==1
semilogx(f,Simu_R,'-k');
hold on;
end
if plot_deflection_heat==1
semilogx(f,abs(Heat_signal),'--k');
hold on
end
if plot_deflection_mass==1
semilogx(f,abs(Mass_signal),':k');
hold on
end
if Importdata==1
semilogx(f_data,R_data,'ok');
hold on;
end 
xlabel('frequency (Hz)','FontSize',12)
    ylabel('Amplitude(\muV)','FontSize',12)
    %title('Amplitude','FontSize',12)
    set(gca,'FontSize',12)
%------------Plot Phase-----------%
figure(P_phase);set(gcf,'units','normalized','position',[0.53 0.38 0.44 0.53])
if plot_all==1
semilogx(f,Simu_Phase,'-g');
hold on;
end
if plot_deflection_heat==1
semilogx(f,angle(Heat_signal)/pi*180,'--g');
hold on
end
if plot_deflection_mass==1
semilogx(f,angle(Mass_signal)/pi*180,':g');
hold on
end
if Importdata==1
semilogx(f_data,Phase_data,'og');
end 

xlabel('frequency (Hz)','FontSize',12)
    ylabel('Phase(degree)','FontSize',12)
    %title('Phase','FontSize',12)
    set(gca,'FontSize',12)
    %axis([min(f) max(f) -75 25])
%--------------Plot In Phase-------------%
figure(P_inphase);set(gcf,'units','normalized','position',[0.04 0.06 0.44 0.53])
if plot_all==1
semilogx(f,Simu_Vin,'-b');
hold on
end
if plot_deflection_heat==1   
semilogx(f,real(Heat_signal),'--b');
hold on
end
if plot_deflection_mass==1
semilogx(f,real(Mass_signal),':b');
hold on;
end
if Importdata==1
semilogx(f_data,Vin_data,'ob');
end 

    xlabel('frequency (Hz)','FontSize',12)
    ylabel('Vin(\muV)','FontSize',12)
    %title('Vin','FontSize',12)
    set(gca,'FontSize',12)
    %axis([min(f) max(f) 0 150])

%-------------Plot Out of Phase-----------%
figure(P_outofphase);set(gcf,'units','normalized','position',[0.53 0.06 0.44 0.53])
if plot_all==1
semilogx(f,Simu_Vout,'-r');
hold on
end
if plot_deflection_heat==1
semilogx(f,imag(Heat_signal),'--r');
hold on
end
if plot_deflection_mass==1
semilogx(f,imag(Mass_signal),':r');
hold on
end
if Importdata==1
semilogx(f_data,Vout_data,'or');
end 

xlabel('frequency (Hz)','FontSize',12)
    ylabel('Vout(\muV)','FontSize',12)
    %title('Vout','FontSize',12)
    set(gca,'FontSize',12)  
end



