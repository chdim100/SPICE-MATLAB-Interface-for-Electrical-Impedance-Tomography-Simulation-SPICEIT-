function [yf,b]=act_lowpass_filt(y,Fs,Windowlength,Order)

L=size(y,2);
K=size(y,1);
% t=linspace(0,L/Fs,L);
fc = 100;
Wn = (2/Fs)*fc;
b = fir1(Windowlength,Wn,'low',kaiser(Windowlength+1,Order));
delay=mean(grpdelay(b,Order,Fs));
% fvtool(b,1,'Fs',Fs)
% tt = t(1:end-round(delay*1.5));
%
for k=1:K
y1(k,:) = filter(b,1,y(k,:));
%yf(k,:) = y1(k,round(delay)+1:end);
yf(k,:) = y1(k,end-50:end);
end

end