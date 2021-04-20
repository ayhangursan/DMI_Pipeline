function outputdata=testDMIfunctions(rawdata,refscan,refNCov,waterloc)

% A script to check deconvolution and denoising
disp(strcat('Processing CSI data:',rawdata));
%Read Philips .data file
currentfolder=cd;
outputdata.datapath=strcat(currentfolder,'\',rawdata);
[outputdata.data, outputdata.list]=csi_loadData(outputdata.datapath);

%% Parameters
outputdata.Parameters.gyromagneticratio=6.53569*10^6;
outputdata.Parameters.Tesla=7;
outputdata.Parameters.Freq=outputdata.Parameters.gyromagneticratio*outputdata.Parameters.Tesla;
outputdata.Parameters.BW=5000;
outputdata.Parameters.ppmwindow=((outputdata.Parameters.Freq/outputdata.Parameters.BW)^-1)*10^6;
outputdata.Parameters.NP=outputdata.list.F_resolution;
outputdata.Parameters.time= linspace(0,outputdata.Parameters.NP/outputdata.Parameters.BW,outputdata.Parameters.NP);
outputdata.Parameters.zerofill= outputdata.Parameters.NP*2;

outputdata.xaxis=linspace(-outputdata.Parameters.ppmwindow/2,outputdata.Parameters.ppmwindow/2,outputdata.Parameters.NP)+4.7;
outputdata.xaxiszerofill=linspace(-outputdata.Parameters.ppmwindow/2,outputdata.Parameters.ppmwindow/2,outputdata.Parameters.zerofill)+4.7;
outputdata.Parameters.Apodizationparameter=5; % Line broadening parameter(in Hertz)
outputdata.Parameters.apodfunc=exp(-outputdata.Parameters.Apodizationparameter*pi*outputdata.Parameters.time).'; %Lorentzian apodization. This value could be decreased as linewidths are a bit broad now.
% Index
outputdata.Parameters.Index.channelIndex = find(contains(outputdata.data.labels,'chan'));
outputdata.Parameters.Index.averageIndex = find(contains(outputdata.data.labels,'aver'));
outputdata.Parameters.Index.kxIndex = find(contains(outputdata.data.labels,'kx'));
outputdata.Parameters.Index.kyIndex = find(contains(outputdata.data.labels,'ky'));
outputdata.Parameters.Index.kzIndex = find(contains(outputdata.data.labels,'kz'));
outputdata.Parameters.dims=size(outputdata.data.raw);
outputdata.Parameters.CSIdims=[outputdata.Parameters.dims(outputdata.Parameters.Index.kxIndex) outputdata.Parameters.dims(outputdata.Parameters.Index.kyIndex) outputdata.Parameters.dims(outputdata.Parameters.Index.kzIndex)];
% Noise covariance matrix!!! additional step may required to get noise from
% the end of the spectrum of fid
% outputdata.NCov=cov(outputdata.data.noise); % Fetch noise data from noise scans
%% Average data
if outputdata.list.number_of_signal_averages==1
    outputdata.avgrawdata=flip(squeeze((outputdata.data.raw)),outputdata.Parameters.Index.kyIndex); disp('NSA=1')
else
    outputdata.avgrawdata=flip(squeeze(mean(outputdata.data.raw,outputdata.Parameters.Index.averageIndex)),outputdata.Parameters.Index.kyIndex);
end
dataorder=1:ndims(outputdata.avgrawdata);
% dataorder([outputdata.Parameters.Index.kyIndex outputdata.Parameters.Index.kxIndex])=dataorder([outputdata.Parameters.Index.kxIndex outputdata.Parameters.Index.kyIndex]);
outputdata.avgrawdata=permute(outputdata.avgrawdata,dataorder);
outputdata.fftfiddata=zeros(size(outputdata.avgrawdata));
outputdata.spectradata=zeros(size(outputdata.avgrawdata));

%% Spatial filter-Hamming for 3D
outputdata.acqpattern=acquisitionpatterncheck(outputdata.data.raw);
if isempty(outputdata.Parameters.Index.channelIndex)==0
    outputdata.acqpattern=squeeze(outputdata.acqpattern(1,:,:,:)/max(outputdata.acqpattern,[],'all'));
end
[outputdata.idealhammingwindow,outputdata.Correctionmask,outputdata.avgrawdata]=Hammingfilterpostcorrection(outputdata.avgrawdata,outputdata.acqpattern,outputdata.Parameters);
%% Apply half-voxel shift in FH
Gridshift=[0 0 -0.5];
outputdata.avgrawdata=Applyvoxelshift(outputdata.avgrawdata,Gridshift,outputdata.Parameters);
disp(strcat('Voxel shift applied. X:',num2str(Gridshift(1)),' Y:',num2str(Gridshift(2)),' Z:',num2str(Gridshift(3))))
%% FFT for each channel
if  isempty(outputdata.Parameters.Index.channelIndex)==1
    outputdata.fftfiddata(:,:,:,:)=Phasecorrection(ifftshift(fftshift(ifft(fftn(squeeze(outputdata.avgrawdata(:,:,:,:))),[],1)),1));
    outputdata.spectradata(:,:,:,:)=fftshift(fft(outputdata.fftfiddata(:,:,:,:),[],1),1);
else
    for n=1:size(outputdata.avgrawdata,outputdata.Parameters.Index.channelIndex)
        outputdata.fftfiddata(:,n,:,:,:)=Phasecorrection(ifftshift(fftshift(ifft(fftn(squeeze(outputdata.avgrawdata(:,n,:,:,:))),[],1)),1));
        outputdata.spectradata(:,n,:,:,:)=fftshift(fft(outputdata.fftfiddata(:,n,:,:,:),[],1),1);
    end
end
disp('Spatial FFT applied.')


%% Channel combination
if isempty(outputdata.Parameters.Index.channelIndex)==0
    
    disp('Input Noise Covariance matrix is used');
    [outputdata.RoemerEqualfid, outputdata.RoemerEqNoiseSensitivity_maps]=RoemerEqualNoise_withmap_input(outputdata.fftfiddata,refscan,refNCov,outputdata.Parameters.Index.channelIndex);        disp('Roemer equeal noise channel combination applied.')
%     [outputdata.WSVD.WSVDcombined, outputdata.WSVD.wsvdQuality, outputdata.WSVD.wsvdWeights]=WSVDcomb(outputdata.fftfiddata,refNCov,outputdata.Parameters.Index.channelIndex);        disp('WSVD channel combination applied.')
end
%% Circular shift
%In RL direction
outputdata.RoemerEqualfid=circshift(outputdata.RoemerEqualfid,-1,3); disp('Circularshift in RL for Roemer comb(fid).')
% outputdata.WSVD.WSVDcombined=circshift(outputdata.WSVD.WSVDcombined,-1,3); disp('Circularshift in RL for WSVD.')

%% Deconvolve combined data
disp('Deconvolution:')
outputdata.RoemerDeconfid=DeconvCSI(outputdata.RoemerEqualfid,outputdata.Parameters);disp('Deconvolution applied for Roemer equal noise combination');
% outputdata.WSVDDeconfid=DeconvCSI(outputdata.WSVD.WSVDcombined,outputdata.Parameters);disp('Deconvolution applied for WSVD combination');
%% SpectralFFT
disp('Spectral FFT')
outputdata.combinedspectradata=fftshift(fft(outputdata.RoemerEqualfid,[],1),1); % To be used in SNR mask for linear prediction
% outputdata.WSVDSpec=fftshift(fft(outputdata.WSVD.WSVDcombined,[],1),1);
outputdata.RoemerDeconSpec=fftshift(fft(outputdata.RoemerDeconfid,[],1),1);
% outputdata.WSVDDeconSpec=fftshift(fft(outputdata.WSVDDeconfid,[],1),1);

%% Calculate SNR for each dataset
disp('Calculating water signal based SNR. Noise is sampled between 20 and 30 ppm.')
noisewindow=find(outputdata.xaxis>20 & outputdata.xaxis<30);

noisemapRoemerSpec=zeros(outputdata.Parameters.CSIdims);
noisemapRoemerDeconvSpec=zeros(outputdata.Parameters.CSIdims);

% noisemapWSVDSpec=zeros(outputdata.Parameters.CSIdims);
% noisemapWSVDDeconvSpec=zeros(outputdata.Parameters.CSIdims);

% Estimate noise in noise window
noisemapRoemerSpec=squeeze(std(real(outputdata.combinedspectradata(noisewindow,:,:,:)))); % Using real part of the spectra to estimate noise in accordance with MRS consensus paper(doi:10.1002/nbm.4347)
noisemapRoemerDeconvSpec=squeeze(std(real(outputdata.RoemerDeconSpec(noisewindow,:,:,:)))); % Using real part of the spectra to estimate noise in accordance with MRS consensus paper(doi:10.1002/nbm.4347)

% noisemapWSVDSpec=squeeze(std(real(outputdata.WSVDSpec(noisewindow,:,:,:)))); % Using real part of the spectra to estimate noise in accordance with MRS consensus paper(doi:10.1002/nbm.4347)
% noisemapWSVDDeconvSpec=squeeze(std(real(outputdata.WSVDDeconSpec(noisewindow,:,:,:)))); % Using real part of the spectra to estimate noise in accordance with MRS consensus paper(doi:10.1002/nbm.4347)

% Water based SNR
waterwindow=find(outputdata.xaxis>4 & outputdata.xaxis<5.4);

outputdata.SNRRoemer=squeeze(max(real(outputdata.combinedspectradata(waterwindow,:,:,:)),[],1))./noisemapRoemerSpec;
outputdata.SNRRoemerDeconv=squeeze(max(real(outputdata.RoemerDeconSpec(waterwindow,:,:,:)),[],1))./noisemapRoemerDeconvSpec;

% outputdata.SNRWSVD=squeeze(max(real(outputdata.WSVDSpec(waterwindow,:,:,:)),[],1))./noisemapWSVDSpec;
% outputdata.SNRWSVDDeconv=squeeze(max(real(outputdata.WSVDDeconSpec(waterwindow,:,:,:)),[],1))./noisemapWSVDDeconvSpec;
%% Calcualte T2star times
disp('Calculating T2^* times')
[outputdata.RoemerT2starmap, outputdata.RoemerHzmap,  ~]=SpectralLinewidth(outputdata.combinedspectradata,outputdata.xaxis,outputdata.Parameters);
% [outputdata.WSVDT2starmap, ~,  ~]=SpectralLinewidth(outputdata.WSVDSpec,outputdata.xaxis,outputdata.Parameters);
[outputdata.RoemerDeconvT2starmap, outputdata.RoemerDeconHzmap,  ~]=SpectralLinewidth(outputdata.RoemerDeconSpec,outputdata.xaxis,outputdata.Parameters);
% [outputdata.WSVDDeconvT2starmap, ~,  ~]=SpectralLinewidth(outputdata.WSVDDeconSpec,outputdata.xaxis,outputdata.Parameters);

% Find water signal for reference
if waterloc=0
    [outputdata.spectra_aligned, outputdata.waterloc]=SpectralFrequencyAlignment(outputdata.combinedspectradata,outputdata.xaxis,outputdata.Parameters,Referenceloc);
    
else
    outputdata.waterloc=waterloc;
end

%% AMARES fit
disp('AMARES fitting')
% % For Roemer equal noise channel combination
% outputdata.AMARES_RoemerEqualfid=DMI_AMARES(outputdata.RoemerEqualfid,outputdata.Parameters,outputdata.xaxis,outputdata.waterloc);
% outputdata.AMARES_RoemerDeconfid=DMI_AMARES(outputdata.RoemerDeconfid,outputdata.Parameters,outputdata.xaxis);

% % For WSVD channel combination
% outputdata.AMARES_WSVD=DMI_AMARES(outputdata.WSVD.WSVDcombined,outputdata.Parameters,outputdata.xaxis,outputdata.waterloc);
% outputdata.AMARES_WSVDDeconfid=DMI_AMARES(outputdata.WSVDDeconfid,outputdata.Parameters,outputdata.xaxis);
end
