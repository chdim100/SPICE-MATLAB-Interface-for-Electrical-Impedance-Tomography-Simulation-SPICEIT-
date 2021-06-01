function SNRmean=Calculate_SNR(VoltageCell,Lbits,Vref,NT,f,fs,Noise)
%%computes DU model SNR
%%%%Inputs:
%%%-VoltageCell: a cell array consisted of SPICE transient voltage measurements results
%%%-Lbits: ADC Resolution
%%%-Vref: ADC reference
%%%-NT: sine periods sampled
%%%-f: input signal frequency
%%%-fs: sampling frequency
%%%-Noise: Noise signal amplitude

%%%Outputs:
%%%-SNRmean: Mean signal SNR of the measurement frame
LM=length(VoltageCell);
A=zeros(LM,1);
for ii=1:LM
    A(ii)=max(VoltageCell{ii})-min(VoltageCell{ii});
end
LSB=Vref/(2^Lbits);
sn=Noise/sqrt(2);
N=NT*fs/f;
SNR=10*log10((mean(A).^2*N/2)/(LSB^2/12+sn^2));
SNRmean=mean(SNR);
end