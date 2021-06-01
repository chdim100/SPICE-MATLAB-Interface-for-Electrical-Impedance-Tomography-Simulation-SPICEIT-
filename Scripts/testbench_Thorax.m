%%%%%%%%%
%%%%SPICEIT_Thorax testbench
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

%%%%measurement pattern (only adjacent, i.e. skip-0, skip-0 is supported
%%%%for the thorax for the time)
%%%skip-current
skipcurr=0;
%%%skip-voltage
skipvolt=0;

%%%%electrode configuration ('Passive' or 'Active_Rec')
%Circuit_Electrodes='Passive';
Circuit_Electrodes='Active_Rec';

%%%electrode deviation ('Symmetric' or 'Asymmetric')
%electrode_deviation='Symmetric';
electrode_deviation='Asymmetric';

%%%electrode size ('005' or '003')
%electrode_size='003';
electrode_size='005';

%%%sinusoidal signal frequency (15000 or 100000 in Hz)
%sig_frequency=15000;
sig_frequency=100000;

%%%ADC sampling rate (SPS)
% fs=60000;
% fs=240000;
%fs=400000;
fs=1600000;

%%%ADC resolution
Lbits=16;

%%%Number of periods sampled per measurement (1-4)
NofSinePeriods=4;

%%%%ADC reference (V)
Vref=3.3;

%%%white noise resulting from LT SPICE simulation (Vrms)
Noise=6e-03; 

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

%%%%getSNR
getcell=1;


%%%%image reconstruction
[im,SNR]=Thoracic_image_with_SPICEIT(path,deskpath,N,skipcurr,skipvolt,...
    Circuit_Electrodes,electrode_deviation,electrode_size,sig_frequency,fs,...
    Lbits,NofSinePeriods,Vref,Noise,Algorithm,prior,lambda,getcell);

%%%%%load reference image for the corresponding frequency
load(['inv_ref_' num2str(round(sig_frequency/1000)) 'kHz.mat'])
%%%% compute image CC (R^2)
CCmat=corrcoef(inv_reference.elem_data,im.elem_data);
CC=CCmat(1,2);
SNR
CC

