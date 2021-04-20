function outputdata=P31Pipeline(rawdata,noisescan)
% Example V9999=P31Pipeline(rawdata,noisescan)
disp(strcat('Processing CSI data:',rawdata));
disp(strcat('Noise scan:',noisescan));

TE=0.56;% Echo time is hard coded
%Read seperate noise scan
outputdata.NoiseScan=NoiseCovarianceGeneration(noisescan);
outputdata.NoiseCov=outputdata.NoiseScan.noisecovariance;
%Read Philips .data file
currentfolder=cd;
outputdata.datapath=strcat(currentfolder,'\',rawdata);
[outputdata.data, outputdata.list]=csi_loadData(outputdata.datapath);

%% Parameters
outputdata.Parameters.gyromagneticratio=17.25144*10^6;
outputdata.Parameters.Tesla=7;
outputdata.Parameters.Freq=outputdata.Parameters.gyromagneticratio*outputdata.Parameters.Tesla;
outputdata.Parameters.BW=5000;disp(strcat('Spectral BW:',num2str(outputdata.Parameters.BW)));
outputdata.Parameters.ppmwindow=((outputdata.Parameters.Freq/outputdata.Parameters.BW)^-1)*10^6;
outputdata.Parameters.NP=outputdata.list.F_resolution;
outputdata.Parameters.time= linspace(0,outputdata.Parameters.NP/outputdata.Parameters.BW,outputdata.Parameters.NP);
outputdata.Parameters.zerofill= outputdata.Parameters.NP*2;

outputdata.xaxis=linspace(-outputdata.Parameters.ppmwindow/2,outputdata.Parameters.ppmwindow/2,outputdata.Parameters.NP);
outputdata.xaxiszerofill=linspace(-outputdata.Parameters.ppmwindow/2,outputdata.Parameters.ppmwindow/2,outputdata.Parameters.zerofill);
outputdata.Parameters.TE=TE*10^-3; % Echo time as an input. Reverse it later !!!
disp(strcat('Acquisition echo time:', num2str(outputdata.Parameters.TE*10^3),' ms'))
outputdata.Parameters.missingpoints=ceil(outputdata.Parameters.BW*outputdata.Parameters.TE);
outputdata.Parameters.Apodizationparameter=5; % Line broadening parameter(in Hertz)
outputdata.Parameters.apodfunc=exp(-pi*outputdata.Parameters.Apodizationparameter.*outputdata.Parameters.time).'; %Lorentzian apodization. This value could be decreased as linewidths are a bit broad now.
% Index
outputdata.Parameters.Index.channelIndex = find(contains(outputdata.data.labels,'chan'));
outputdata.Parameters.Index.averageIndex = find(contains(outputdata.data.labels,'aver'));
outputdata.Parameters.Index.kxIndex = find(contains(outputdata.data.labels,'kx'));
outputdata.Parameters.Index.kyIndex = find(contains(outputdata.data.labels,'ky'));
outputdata.Parameters.Index.kzIndex = find(contains(outputdata.data.labels,'kz'));
outputdata.Parameters.dims=size(outputdata.data.raw);
outputdata.Parameters.CSIdims=[outputdata.Parameters.dims(outputdata.Parameters.Index.kxIndex) outputdata.Parameters.dims(outputdata.Parameters.Index.kyIndex) outputdata.Parameters.dims(outputdata.Parameters.Index.kzIndex)];
outputdata.Parameters.FirstOrdPhaseFunct=exp(-1i * (2* pi * (outputdata.xaxis).'*(outputdata.Parameters.Freq/(10^6)) * outputdata.Parameters.TE));

% Average data
if outputdata.list.number_of_signal_averages==1
    outputdata.avgrawdata=squeeze((outputdata.data.raw)); disp('NSA=1')
else
    outputdata.avgrawdata=squeeze(mean(outputdata.data.raw,outputdata.Parameters.Index.averageIndex));
end

dataorder=1:ndims(outputdata.avgrawdata);
outputdata.avgrawdata=permute(outputdata.avgrawdata,dataorder);
% outputdata.avgrawdata=flip(outputdata.avgrawdata,5);
% Spatial filter-Hamming for 3D
outputdata.acqpattern=acquisitionpatterncheck(outputdata.data.raw);

if isempty(outputdata.Parameters.Index.channelIndex)==0
    outputdata.acqpattern=squeeze(outputdata.acqpattern(1,:,:,:)/max(outputdata.acqpattern,[],'all'));
end

[outputdata.idealhammingwindow,outputdata.Correctionmask,outputdata.avgrawdata]=Hammingfilterpostcorrection(outputdata.avgrawdata,outputdata.acqpattern,outputdata.Parameters);
outputdata.avgrawdata=permute(outputdata.avgrawdata,[1 2 4 3 5]);
outputdata.Parameters.CSIdims=[outputdata.Parameters.dims(outputdata.Parameters.Index.kyIndex) outputdata.Parameters.dims(outputdata.Parameters.Index.kxIndex) outputdata.Parameters.dims(outputdata.Parameters.Index.kzIndex)];

outputdata.avgrawdata=flip(outputdata.avgrawdata,3);disp('Fliped in AP')
outputdata.avgrawdata=flip(outputdata.avgrawdata,4);disp('Fliped in RL')

% FFT for each channel
outputdata.fftfiddata=zeros(size(outputdata.avgrawdata));
outputdata.spectradata=zeros(size(outputdata.avgrawdata));
for n=1:size(outputdata.avgrawdata,outputdata.Parameters.Index.channelIndex)
    outputdata.fftfiddata(:,n,:,:,:)=Phasecorrection(SpatialFFT(squeeze(outputdata.avgrawdata(:,n,:,:,:))));
    outputdata.spectradata(:,n,:,:,:)=fftshift(fft(outputdata.fftfiddata(:,n,:,:,:),[],1),1);
end
disp('Spatial FFT applied.')

% Phase each spectra
outputdata.AutoPhasedFID=Phase31PSpectra(outputdata.fftfiddata,outputdata.xaxis,outputdata.Parameters);
outputdata.AutophasedSPEC=fftshift(fft(outputdata.AutoPhasedFID,[],1),1).*outputdata.Parameters.FirstOrdPhaseFunct;

outputdata.PhasedFIDQ2=Phase31PSpectraQ2(outputdata.fftfiddata,outputdata.Parameters);
outputdata.AutoPhasedSPECQ2=fftshift(fft(outputdata.PhasedFIDQ2,[],1),1).*outputdata.Parameters.FirstOrdPhaseFunct;

disp('Automatic phasing applied.')

% Channel combination
[outputdata.RoemerEqualfid, outputdata.RoemerEqNoiseSensitivity_maps]=RoemerEqualNoise_withmap_input(outputdata.AutoPhasedFID,0,outputdata.NoiseCov,outputdata.Parameters.Index.channelIndex);        disp('Roemer equeal noise channel combination applied.')
[outputdata.RoemerEqualfidQ2, outputdata.RoemerEqNoiseSensitivity_mapsQ2]=RoemerEqualNoise_withmap_input(outputdata.PhasedFIDQ2,0,outputdata.NoiseCov,outputdata.Parameters.Index.channelIndex);        disp('Roemer equeal noise channel combination applied.')

[outputdata.WSVD.WSVDcombined, outputdata.WSVD.wsvdQuality, outputdata.WSVD.wsvdWeights]=WSVDcomb(outputdata.AutophasedSPEC,outputdata.NoiseCov,outputdata.Parameters.Index.channelIndex);        disp('WSVD channel combination applied.')
[outputdata.WSVDQ2.WSVDcombined, outputdata.WSVDQ2.wsvdQuality, outputdata.WSVDQ2.wsvdWeights]=WSVDcomb(outputdata.AutoPhasedSPECQ2,outputdata.NoiseCov,outputdata.Parameters.Index.channelIndex);        disp('WSVD channel combination applied.')

[outputdata.WSVDNoPhasing.WSVDcombined, outputdata.WSVDNoPhasing.wsvdQuality, outputdata.WSVDNoPhasing.wsvdWeights]=WSVDcomb(outputdata.spectradata,outputdata.NoiseCov,outputdata.Parameters.Index.channelIndex);        disp('WSVD channel combination applied.')

outputdata.WSVD.WSVDcombined=fftshift(fft(Phase31PSpectra(ifft(ifftshift(outputdata.WSVD.WSVDcombined./outputdata.Parameters.FirstOrdPhaseFunct,1),[],1),outputdata.xaxis,outputdata.Parameters),[],1),1).*outputdata.Parameters.FirstOrdPhaseFunct;
outputdata.WSVDQ2.WSVDcombined=fftshift(fft(Phase31PSpectraQ2(ifft(ifftshift(outputdata.WSVDQ2.WSVDcombined./outputdata.Parameters.FirstOrdPhaseFunct,1),[],1),outputdata.Parameters),[],1),1).*outputdata.Parameters.FirstOrdPhaseFunct;

disp('Channel combination applied.')

% Circular shift
%In RL direction
outputdata.RoemerEqualfid=circshift(outputdata.RoemerEqualfid,-1,2); disp('Circularshift in AP(-1) for Roemer.')
outputdata.WSVD.WSVDcombined=circshift(outputdata.WSVD.WSVDcombined,-1,2); disp('Circularshift in AP(-1) for WSVD.')
outputdata.WSVDQ2.WSVDcombined=circshift(outputdata.WSVDQ2.WSVDcombined,-1,2); disp('Circularshift in AP(-1) for WSVD.')

outputdata.WSVDNOphasing.WSVDcombined=circshift(outputdata.WSVDQ2.WSVDcombined,-1,2); disp('Circularshift in AP(-1) for WSVD.')


outputdata.Roemercombinedspectra=fftshift(fft(Phase31PSpectra(outputdata.RoemerEqualfid,outputdata.xaxis,outputdata.Parameters),[],1),1).*outputdata.Parameters.FirstOrdPhaseFunct; % To be used in SNR mask for linear prediction
% outputdata.WSVDcombinedspectra=fftshift(fft(outputdata.WSVD.WSVDcombined,[],1),1); % To be used in SNR mask for linear prediction

%% PCr based SNR

outputdata.noisewindow=find(outputdata.xaxis>10 & outputdata.xaxis<20);
outputdata.PCrwindow=find(outputdata.xaxis>0 & outputdata.xaxis<2.5);

outputdata.noisemap.Roemercombinedspectra=squeeze(std(real(outputdata.Roemercombinedspectra(outputdata.noisewindow,:,:,:))));
outputdata.noisemap.WSVDcombined=squeeze(std(real(outputdata.WSVD.WSVDcombined(outputdata.noisewindow,:,:,:))));
outputdata.noisemap.WSVDcombinedQ2=squeeze(std(real(outputdata.WSVDQ2.WSVDcombined(outputdata.noisewindow,:,:,:))));

outputdata.SNR.Roemercombinedspectra=squeeze(max(real(outputdata.Roemercombinedspectra(outputdata.PCrwindow,:,:,:)),[],1))./outputdata.noisemap.Roemercombinedspectra;
outputdata.SNR.WSVDcombined=squeeze(max(real(outputdata.WSVD.WSVDcombined(outputdata.PCrwindow,:,:,:)),[],1))./outputdata.noisemap.WSVDcombined;

%% Export as spar/sdat
% Not included yet

end