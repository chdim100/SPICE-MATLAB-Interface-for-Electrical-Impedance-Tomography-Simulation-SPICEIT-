function [Amplitudes,Real,Imag,Phases,VoltageCell]=Get_Spice_Measurements(N,skipcurr,skipvolt,f,RESADC,fSample,subjectfile,Electrode_pos,Circuit_Electrodes,electrode_size,deskpath,Noise,NofSinePeriods,VADC,casetest,getcell,show_transient)

switch casetest
    case 'Thoracic'
        %confirm electrode size chars
        switch electrode_size
            case '005'
                sizel='';
            otherwise
                sizel='_Rel=003';
        end
        filemeas=[deskpath num2str(f/1000),'k\Pt_Bright\',Circuit_Electrodes,'\',subjectfile,...
            '_',Electrode_pos,'_',Circuit_Electrodes,'_','NonBalancedCaps','_with_skin',sizel,'.txt'];
    case 'Tank'
        filemeas=[deskpath 'TankConf\' subjectfile '.txt'];
    otherwise
        error('Under Update')
end
%==========================================================================
%1. Get Spice Data
%open Times and Transient Spice files
fileTimes=[deskpath 'Times\Times_',num2str(N),'_skip-',...
    num2str(skipcurr),'_skip-',num2str(skipvolt),'_',num2str(f/1000),'kHz.txt'];
Times=dlmread(fileTimes,'\t',0,0);

Values=dlmread(filemeas,'\t',1,0);
%Spice Transient Points
Spice_Transient_Voutput=Values(:,2);
Spice_Transient_Ioutput=Values(:,3);
Spice_Transient_Timepoints=Values(:,1);

%==========================================================================
%2. Initialize DAQ and some other specs
%2A. DAQ Specs
%ADC Voltage Swing
% VADC=3.3;
%LEAST SIGNIFICANT BITS
LSB=2;
%ADC ENOB
ACTUALRESADC=RESADC-LSB;
%Number of periods per window used for IQ
NT=NofSinePeriods;

%2B. Set the number of voltage measurements per current source electrode
%position
if skipcurr==skipvolt
    meas=N-3;
else
    meas=N-4;
end

%==========================================================================
%3.Place the Spice transient points in the corresponding measurement window
%3A. Initialize the Spice transient points per measurement matrices
point_per_measurement=round(fSample*(5/f))+1;
Vout_samples_permeasurement=zeros(N*meas,point_per_measurement);
Iout_samples_permeasurement=zeros(N*meas,point_per_measurement);

%run loop in measurements to propely complete those matrices
%set a counter (ii)
ii=1;
VoltageCell={};
for currentsource=1:N
    for voltagemeasurement=1:meas
        %the Times matrix is precomputed and indicates when to start and
        %stop each measuring window (avoiding transient effects)
        %Times is different for each frequency and measuring pattern
        
        %3.B. Set the start and stop times for this measuring window
        %(iith): When to start and stop getting Spice transient points
        %to complete the ith measure
        Tstart=Times(2*ii-1,1);
        Tstop=Times(2*ii,1);
        
        %3.C. find the corresponding indices of spice timepoints that are
        %included in the integral [Tstart, Tstop]
        indices_of_spicepoints=find(Spice_Transient_Timepoints>=Tstart&Spice_Transient_Timepoints<=Tstop);
        %add a last indice, +1 from the previous last
        last=indices_of_spicepoints(end);
        indices_of_spicepoints=[indices_of_spicepoints; last+1];
        
        %3.D.Get the corresponding values of Voltage (and current) outputs from
        %the Spice transients that are in [Tstart, Tstop]. Do same for
        %Spice transient Timepoints
        
        Voltage_Spice_transient_Window=Spice_Transient_Voutput(indices_of_spicepoints);
        Current_Spice_transient_Window=Spice_Transient_Ioutput(indices_of_spicepoints);
        Spice_Timepoints_Window=Spice_Transient_Timepoints(indices_of_spicepoints);
        
        %3.E. Add the precomputed noise to them:
        %Noise to the Voltage output
        Noisesignal=add_white_noise(length(Voltage_Spice_transient_Window),Noise);
        if getcell==1
            VoltageCell{end+1}=Voltage_Spice_transient_Window;
        end
        Noisy_Voltage_Spice_transient_Window=Voltage_Spice_transient_Window+Noisesignal;
        %Noise to the Current output (here we just add awgn of 70-90dB)
        Noisy_Current_Spice_transient_Window=awgn(Current_Spice_transient_Window,70);
        select_inj=[1,2,3,4,5,6,7,8];
        plottransient=1;
        
        %         for ils=1:length(select_inj)
        %         %Plot the window if selected so:
        %         if plottransient&&ii==select_inj(ils)
        %             currentsourceplotted=currentsource;
        %             voltagemeasurementplotted=voltagemeasurement;
        %             figure
        %             subplot(1,2,1)
        %             plot(Spice_Timepoints_Window,Noisy_Voltage_Spice_transient_Window)
        %             Title=title(['Voltage Transient Signal, Current:' num2str(currentsource) ...
        %                 ' Voltage:' num2str(voltagemeasurement)],'Interpreter','Latex');
        %             Title.FontSize=16;
        %             hold on
        %         end
        %         end
        
        %ADC
        
        %3.F. Since the LT Spice timeline does not have a constant step
        %we find the nearest LTspice timepoints to the timepoints we sample!
        %ideal ADC Samplepoints
        Ideal_ADCSample_Window_timepoints=Tstart:1/fSample:(point_per_measurement-1)/fSample+Tstart;
        %REPMAT ideal ADC Samplepoints and Spice Timepoints of this Window
        %in order to be able to compare them
        A = repmat(Ideal_ADCSample_Window_timepoints',[1 length(Spice_Timepoints_Window)]);
        B=repmat(Spice_Timepoints_Window,[1 length(Ideal_ADCSample_Window_timepoints)]);
        %find the indices of Spice Window Timepoints that are closest to ideal ADC Samplepoints
        [minValue,closestIndex] = min(abs(A'-B));
        %get the corresponding Spice Window Timepoints (Values in sec):
        closest_Spice_timepoint_Value=Spice_Timepoints_Window(closestIndex);
        %closestValue is the Spice Window timepoints we actually sample!!!!
        %closestIndex is the indices of Spice timepoints we sample!!!
        %so...
        
        %3.G. Get the Voltage and Current values that are to be sampled
        Noisy_Voltage_Spice_transient_Window=(10+0.05*randn(1))*(Noisy_Voltage_Spice_transient_Window-mean(Noisy_Voltage_Spice_transient_Window))+VADC/2;
        Vout_samples_permeasurement(ii,:)=Noisy_Voltage_Spice_transient_Window(closestIndex);
        Iout_samples_permeasurement(ii,:)=Noisy_Current_Spice_transient_Window(closestIndex);
        %and we got the samples!!!
        %plot them on the transient scheme if selected
        hh =  findobj('type','figure');
        nn = length(hh);
        if show_transient==1
            for ils=1:length(select_inj)
                if plottransient&&ii==select_inj(ils)
                    currentsourceplotted(ils)=currentsource;
                    voltagemeasurementplotted(ils)=voltagemeasurement;
                    figure(nn+1)
                    nno=nn+1;
                    subplot(1,2,1)
                    plot(Spice_Timepoints_Window,Noisy_Voltage_Spice_transient_Window)
                    Title=title(['Voltage Signal, Current:' num2str(currentsource) ...
                        ' Voltage:' num2str(voltagemeasurement)],'Interpreter','Latex');
                    Title.FontSize=16;
                    hold on
                    plot(closest_Spice_timepoint_Value,Vout_samples_permeasurement(ii,:),'r*','MarkerSize',10)
                    hold on
                    xlim([min(Spice_Timepoints_Window) max(Spice_Timepoints_Window)]);
                    YL1=ylabel('Voltage ($V$)','Interpreter','Latex');
                    YL1.FontSize=14;
                    XL1=xlabel('Time ($s$)','Interpreter','Latex');
                    XL1.FontSize=14;
                    LL=legend({'Voltage','Non-Quantized Samples'});
                    LL.FontSize=12;
                    grid on
                    att = gca;
                    att.FontSize = 12;
                end
            end
        end
        ii=ii+1;
    end
end

%==========================================================================
%4. ADC Sample and quantization

%4A. set each Voltage output window in the ADC voltage swing range
%(this can be done with the filter and the ADC driver)
% Vout_samples_permeasurement_ranged=((Vout_samples_permeasurement-...
%     mean(mean(Vout_samples_permeasurement))))*VADC+VADC/2;
Vout_samples_permeasurement_ranged=Vout_samples_permeasurement;
%check and rerange if any voltage is <0 (out of ADC's range):
while ~isempty(Vout_samples_permeasurement_ranged(Vout_samples_permeasurement_ranged<=0))
    Vout_samples_permeasurement_ranged=0.8*Vout_samples_permeasurement_ranged+0.3;
end
%current measurements are already in ADC voltage swing range!

%4.B. quantize sampled values
Vout_samples_permeasurement_quantized=round(Vout_samples_permeasurement_ranged/VADC*(2^RESADC-1));
Iout_samples_permeasurement_quantized=round(Iout_samples_permeasurement/VADC*(2^RESADC-1));

%4.C. create the sin input signal (plus 90 degrees)
Iout_sinus_samples_permeasurement_quantized=Iout_samples_permeasurement_quantized(:,round(fSample/f/4)+1:end);

%4.D. take the sampled and quantized signals to binary form
Vout_samples_permeasurement_quantized_binary=dec2bin(Vout_samples_permeasurement_quantized);
Iout_samples_permeasurement_quantized_binary=dec2bin(Iout_samples_permeasurement_quantized);
Iout_sinus_samples_permeasurement_quantized_binary=dec2bin(Iout_sinus_samples_permeasurement_quantized);

%4.E. throw out noisy LSB
Vout_samples_permeasurement_quantized_binary_reduced=Vout_samples_permeasurement_quantized_binary(:,1:end-LSB);
Iout_samples_permeasurement_quantized_binary_reduced=Iout_samples_permeasurement_quantized_binary(:,1:end-LSB);
Iout_sinus_samples_permeasurement_quantized_binary_reduced=Iout_sinus_samples_permeasurement_quantized_binary(:,1:end-LSB);

%4.G. take the signals in decimal form
%output voltage
Vout_samples_permeasurement_quantized_decimal_reduced=...
    reshape(bin2dec(Vout_samples_permeasurement_quantized_binary_reduced),...
    [N*meas size(Vout_samples_permeasurement_ranged,2)]);
%output current cos
Iout_samples_permeasurement_quantized_decimal_reduced=...
    reshape(bin2dec(Iout_samples_permeasurement_quantized_binary_reduced),...
    [N*meas size(Iout_samples_permeasurement_quantized,2)]);
%output current sin
Iout_sinus_samples_permeasurement_quantized_decimal_reduced=...
    reshape(bin2dec(Iout_sinus_samples_permeasurement_quantized_binary_reduced),...
    [N*meas size(Iout_sinus_samples_permeasurement_quantized,2)]);

%4.H. take the integer number of periods selected ([1,4])
Vout_samples_quantized_final=Vout_samples_permeasurement_quantized_decimal_reduced(:,1:round(NT*fSample/f)+1);
Iout_samples_quantized_final=Iout_samples_permeasurement_quantized_decimal_reduced(:,1:round(NT*fSample/f)+1);
Iout_sinus_samples_quantized_final=Iout_sinus_samples_permeasurement_quantized_decimal_reduced(:,1:round(NT*fSample/f)+1);

%plot if desired
if show_transient==1
    if plottransient
        for ils=1:length(select_inj)
            figure(nno-length(select_inj)+ils)
            subplot(1,2,2)
            plot(Vout_samples_quantized_final(select_inj(ils),:),'-s','MarkerSize',10,...
                'MarkerEdgeColor','red',...
                'MarkerFaceColor',[1 .6 .6])
            Title2=title(['ADC quantized Samples, Current:' num2str(currentsourceplotted(ils)) ...
                ' Voltage:' num2str(voltagemeasurementplotted(ils))],'Interpreter','Latex');
            Title2.FontSize=16;
            YL1=ylabel('ADC Value','Interpreter','Latex');
            YL1.FontSize=14;
            XL1=xlabel('No of Samples','Interpreter','Latex');
            XL1.FontSize=14;
            grid on
            att = gca;
            att.FontSize = 12;
        end
        
    end
end


%4.I. write samples at voltage form
Vout_samples_final=Vout_samples_quantized_final*VADC/(2^ACTUALRESADC-1);
Iout_samples_final=Iout_samples_quantized_final*VADC/(2^ACTUALRESADC-1);
Iout_sinus_samples_final=Iout_sinus_samples_quantized_final*VADC/(2^ACTUALRESADC-1);

%==========================================================================

%5. Interpolation and Multiplication
NN=1;
for mm=1:N*meas
    %for each measuring Window:
    
    %5.A. upsampling NN points
    if NN>1
        X=1:length(Vout_samples_final(1,:));
        NY1n=upsampling(Vout_samples_final(mm,:),X,NN);
        NY2n=upsampling(Iout_samples_final(mm,:),X,NN);
        NY3n=upsampling(Iout_sinus_samples_final(mm,:),X,NN);
    else
        NY1n=Vout_samples_final(mm,:);
        NY2n=Iout_samples_final(mm,:);
        NY3n=Iout_sinus_samples_final(mm,:);
    end
    
    %5.B. Multiply elementwise
    Cosf_expanded(mm,:)=(NY1n-mean(NY1n)).*(NY2n-mean(NY2n));
    Sinf_expanded(mm,:)=(NY1n-mean(NY1n)).*(NY3n-mean(NY2n));
end

%==========================================================================

%6. Act Lowpass matched digital filter, GET IQ

%6.A. Act matched filters

% [Cosfiltered,b]=act_lowpass_filt(Cosf_expanded,fSample,min(950,length(Cosf_expanded(1,:))),10);
% [Sinfiltered,b]=act_lowpass_filt(Sinf_expanded,fSample,min(950,length(Sinf_expanded(1,:))),10);

%6.B. Get IQ components
Cosf=mean(Cosf_expanded,2);
Sinf=mean(Sinf_expanded,2);
Amplitudes=sqrt(Cosf.^2+Sinf.^2);
Real=Cosf;
Imag=Sinf;
Phases=atan2(Sinf,Cosf);
end