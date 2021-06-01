%%%%%%%%%
%%%%SPICEIT Phantom Tank testbench
%%%%%%%%%
%%%path to eidors ('...\eidors-v3.9-ng\eidors\')
path2EID='C:\Users\Chris\Desktop\PHD_files\eidors-v3.9-ng\eidors\';
%%%start EIDORS
startpath=[path2EID 'startup'];
run(startpath)

%%%path to SPICEIT
path2SPICEIT='D:\EIT\SPICEIT\';
%%%%path to measurement files (...'LT_SPICE\Measurements\')
deskpath=[path2SPICEIT 'LT_SPICE\Measurements\'];

%%%No of electrodes
N=16;

%%%%measurement pattern (skip-0, skip-0 for both 15k and 100k freqs
%%%%, or skip-3, skip-0 for only 15kHz freqs)
%%%skip-current
skipcurr=0;
%%%skip-voltage
skipvolt=0;

%%%%electrode configuration (ONLY 'Passive')
Circuit_Electrodes='Passive';

%%%electrode deviation (ONLY 'Symmetric' (no deviation))
electrode_deviation='Symmetric';

%%%electrode size (ONLY '005')
electrode_size='005';

%%%sinusoidal signal frequency (15000 or 100000 in Hz). %%%NOTE THAT 100KHz
%%%are not available for current skip-3, voltage skip-0 pattern
sig_frequency=15000;
%sig_frequency=100000;

%%%ADC sampling rate (MSPS)
fs=1e+05;

%%%ADC resolution
Lbits=16;

%%%Number of periods sampled per measurement (1-4)
NofSinePeriods=4;

%%%ADC Voltage Reference
Vref=3.3;

%%%white noise resulting from LT SPICE simulation (Vrms)
Noise=1e-03; 

%%%reconstruction algorithm================================================

%%%-Algorithm: 'direct' or 'one_step' for single-step approach, 
%%% 'TV' or 'Total_Variation' or 'Pdipm' or 'pdipm' or 'PDIPM'
%%% for Total Variation 
%%% 'iterationalGN' or 'mulistepGN' or 'iterationGN' or 'GN2' for iterative
%%% (non-linear) Gauss-Newton approach

Algorithm='direct';

%%%reconstruction prior====================================================

%%%-prior: 1,'NOSER','Noser' or 'noser' for NOSER,
%%% 2,'LAPLACE','Laplace' or'laplace' for Laplace
%%% 3,'Gaussian','Gauss','GAUSSIAN' or 'gaussian' for Gaussian HPF
%%% 4,'TV' or 'Total_Variation' for Total Variation
%%% lambda: hyperparameter of the inverse problem

prior='Noser';

%%%%reconstruction hyperparameter
lambda=0.1;

%%%%%calculate SNR
getcell=1;

situation1='Homo';
situation2='inhomo';
%%%%
%%%tank case 1-3
tankcase=2;

measdest_homo=['Tank1Homo_Skip' num2str(skipcurr)];
measdest_inhomo=['Tank' num2str(tankcase) 'inhomo_Skip' num2str(skipcurr)];

%%%%deflated case measurements
[Amplitudeshomo,Realhomo,Imagehomo,Phaseshomo,VoltageCell1]=Get_Spice_Measurements(N,skipcurr,...
    skipvolt,sig_frequency,Lbits,fs,measdest_homo,electrode_deviation,...
    Circuit_Electrodes,electrode_size,deskpath,Vref,Noise,NofSinePeriods,Vref,'Tank',getcell,1);
%%%%inflated case measurements
[Amplitudesinh,Realinh,Imaginh,Phasesinh,VoltageCell2]=Get_Spice_Measurements(N,skipcurr,...
    skipvolt,sig_frequency,Lbits,fs,measdest_inhomo,electrode_deviation,...
    Circuit_Electrodes,electrode_size,deskpath,Vref,Noise,NofSinePeriods,Vref,'Tank',getcell,0);

%%%%set reconstruction model
target='d2c2';
IMDL=set_inverse_model(N,skipcurr,skipvolt,Algorithm,prior,lambda,target);

Imreal=inv_solve(IMDL,Realhomo/max(Realhomo),Realinh/max(Realhomo));
Maxscale=max(real(Imreal.elem_data));
%Minscale=min(real(Imreal.elem_data));
Imreal.elem_data=(Imreal.elem_data)/Maxscale*0.3;
im=Imreal;
figure
H1=show_fem(Imreal,1);
axis off
set(H1, 'edgecolor', 'none');
axis square
figure
plot(Amplitudeshomo/max(Amplitudeshomo),'LineWidth',2)
hold on
plot(Amplitudesinh/max(Amplitudesinh),'LineWidth',2)
L=legend({'Homogeneous','Inhomogeneous'},'interpreter','latex');
L.FontSize=12;
XL1=xlabel('Measurement Index','interpreter','latex');
XL1.FontSize=14;
YL1=ylabel('Amplitudes (Normalized)','interpreter','latex');
YL1.FontSize=14;
xlim([0 210])
ylim([0 1.3])


