function [im,SNR]=Thoracic_image_with_SPICEIT(path,deskpath,N,skipcurr,skipvolt,...
    Circuit_Electrodes,electrode_deviation,electrode_size,sig_frequency,fs,...
    Lbits,NofSinePeriods,Vref,Noise,Algorithm,prior,lambda,getcell)
%%%%%%acquire thoracic images
%%% Inputs:
%%%-path: the path to EIDORS
%%%-deskpath: the path to the transient SPICE measurement files (folder
%%%Measurements)
%%%-N: the number of electrodes used for the measurements
%%%-skipcurr: current-skip protocol (if available)
%%%-skipvolt: voltage-skip protocol (if available)
%%%-Circuit_Electrodes: 'Passive' for passive electrode configuration,
%%%'Active_Rec' for active electrode configuration
%%%-electrode_deviation: 'Symmetric' for no electrode deviation,
%%%'Asymmetric' for minor electrode deviation
%%%-electrode_size: '005' or '003' (5 or 3cm) radius circular electrodes
%%%-sig_frequency: input signal frequency 15000 or 100000 (in Hz)
%%%-fs: ADC sampling frequency (in Hz)
%%%-Lbits: ADC resolution
%%%-NofSinePeriods: No of periods to be sampled (integer, 1-4) 
%%%-Noise: white noise signal (in Volts), as computed from SPICE simulation
%%%-Algorithm: 'direct' or 'one_step' for single-step approach, 
%%% 'TV' or 'Total_Variation' or 'Pdipm' or 'pdipm' or 'PDIPM'
%%% for Total Variation 
%%% 'iterationalGN' or 'mulistepGN' or 'iterationGN' or 'GN2' for iterative
%%% (non-linear) Gauss-Newton approach
%%%-prior: 1,'NOSER','Noser' or 'noser' for NOSER,
%%% 2,'LAPLACE','Laplace' or'laplace' for Laplace
%%% 3,'Gaussian','Gauss','GAUSSIAN' or 'gaussian' for Gaussian HPF
%%% 4,'TV' or 'Total_Variation' for Total Variation
%%% lambda: hyperparameter of the inverse problem

%%%Outputs:
%%%-im: reconstructed image struct (EIDORS format)

%%%============================================================================
%%%============================================================================
situation1='deflated';
situation2='inflated';
if skipvolt==0&&skipcurr==0
    skip='';
else
    skip=['skip',num2str(skipcurr),'_',num2str(skipvolt),'_'];
end

measdest_deflated=['EIT_2019_',num2str(round(sig_frequency/1000)),...
    'k_',skip,'Mirrored_Current_Pt_Bright_VAR' situation1];
measdest_inflated=['EIT_2019_',num2str(round(sig_frequency/1000)),...
    'k_',skip,'Mirrored_Current_Pt_Bright_VAR' situation2];

%%%%deflated case measurements
[Amplitudesdef,Realdef,Imagdef,Phasesdef,VoltageCell1]=Get_Spice_Measurements(N,skipcurr,...
    skipvolt,sig_frequency,Lbits,fs,measdest_deflated,electrode_deviation,...
    Circuit_Electrodes,electrode_size,deskpath,Noise,NofSinePeriods,Vref,'Thoracic',getcell,0);
%%%%inflated case measurements
[Amplitudesinf,Realinf,Imaginf,Phasesinf,VoltageCell2]=Get_Spice_Measurements(N,skipcurr,...
    skipvolt,sig_frequency,Lbits,fs,measdest_inflated,electrode_deviation,...
    Circuit_Electrodes,electrode_size,deskpath,Noise,NofSinePeriods,Vref,'Thoracic',getcell,0);
%%%%set reconstruction model
target='d2T3';
IMDL=set_inverse_model(N,skipcurr,skipvolt,Algorithm,prior,lambda,target);
IMDL.fwd_model.electrode([9:-1:1,16:-1:10])=IMDL.fwd_model.electrode;

Imreal=inv_solve(IMDL,Realinf/max(Amplitudesdef),Realdef/max(Amplitudesdef));
Maxscale=max(real(Imreal.elem_data));
%Minscale=min(real(Imreal.elem_data));
Imreal.elem_data=(Imreal.elem_data)/Maxscale*0.3;
im=Imreal;

if getcell==1
%%%%%%Calculate mean signal SNR
SNRmean1=Calculate_SNR(VoltageCell1,Lbits,Vref,NofSinePeriods,sig_frequency,fs,Noise);
SNRmean2=Calculate_SNR(VoltageCell2,Lbits,Vref,NofSinePeriods,sig_frequency,fs,Noise);
SNR=(SNRmean1+SNRmean2)/2;
else
    SNR=[];
end

figure
H1=show_fem(Imreal,1);
axis off
set(H1, 'edgecolor', 'none');
axis square
figure
plot(Amplitudesdef/max(Amplitudesdef),'LineWidth',2)
hold on
plot(Amplitudesinf/max(Amplitudesinf),'LineWidth',2)
L=legend({'Deflated Lungs','Inflated Lungs'},'interpreter','latex');
L.FontSize=12;
XL1=xlabel('Measurement Index','interpreter','latex');
XL1.FontSize=14;
YL1=ylabel('Amplitudes (Normalized)','interpreter','latex');
YL1.FontSize=14;
xlim([0 210])
ylim([0 1.3])


end