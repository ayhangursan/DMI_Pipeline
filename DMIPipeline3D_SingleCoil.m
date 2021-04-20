function outputdata=DMIPipeline3D_SingleCoil(rawdata,TE,Gridshift,Referenceloc)
%%
disp(strcat('Processing CSI data:',rawdata));
%Read Philips .data file
currentfolder=cd;
outputdata.datapath=strcat(currentfolder,'\',rawdata);
[outputdata.data, outputdata.list]=csi_loadData(outputdata.datapath);

%% Parameters
outputdata.Parameters.gyromagneticratio=6.53569*10^6;
outputdata.Parameters.Tesla=7;
outputdata.Parameters.Freq=outputdata.Parameters.gyromagneticratio*outputdata.Parameters.Tesla;
outputdata.Parameters.BW=5000;disp(strcat('Spectral BW:',num2str(outputdata.Parameters.BW)));
outputdata.Parameters.ppmwindow=((outputdata.Parameters.Freq/outputdata.Parameters.BW)^-1)*10^6;
outputdata.Parameters.NP=outputdata.list.F_resolution;
outputdata.Parameters.time= linspace(0,outputdata.Parameters.NP/outputdata.Parameters.BW,outputdata.Parameters.NP);
outputdata.Parameters.zerofill= outputdata.Parameters.NP*2;

outputdata.xaxis=linspace(-outputdata.Parameters.ppmwindow/2,outputdata.Parameters.ppmwindow/2,outputdata.Parameters.NP)+4.7;
outputdata.xaxiszerofill=linspace(-outputdata.Parameters.ppmwindow/2,outputdata.Parameters.ppmwindow/2,outputdata.Parameters.zerofill)+4.7;
outputdata.Parameters.TE=TE*10^-3; % Echo time as an input. Reverse it later !!!
disp(strcat('Acquisition echo time:', num2str(outputdata.Parameters.TE*10^3),'ms'))
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
% Noise covariance matrix!!! additional step may required to get noise from
% the end of the spectrum of fid
% outputdata.NCov=cov(outputdata.data.noise); % Fetch noise data from noise scans
%% Average data
disp('Averaging data.')
if outputdata.list.number_of_signal_averages==1
    outputdata.avgrawdata=flip(squeeze((outputdata.data.raw)),outputdata.Parameters.Index.kyIndex); disp('NSA=1')
else
    outputdata.avgrawdata=flip(squeeze(mean(outputdata.data.raw,outputdata.Parameters.Index.averageIndex)),outputdata.Parameters.Index.kyIndex);
end
dataorder=1:ndims(outputdata.avgrawdata);
outputdata.avgrawdata=permute(outputdata.avgrawdata,dataorder);
outputdata.fftfiddata=zeros(size(outputdata.avgrawdata));
outputdata.spectradata=zeros(size(outputdata.avgrawdata));

%% Spatial filter-Hamming for 3D
outputdata.acqpattern=acquisitionpatterncheck(outputdata.data.raw);
[outputdata.idealhammingwindow,outputdata.Correctionmask,outputdata.avgrawdata]=Hammingfilterpostcorrection(outputdata.avgrawdata,outputdata.acqpattern,outputdata.Parameters);

%% Apply half-voxel shift in FH
outputdata.avgrawdata=Applyvoxelshift(outputdata.avgrawdata,Gridshift,outputdata.Parameters);
disp(strcat('Voxel shift applied. X:',num2str(Gridshift(1)),' Y:',num2str(Gridshift(2)),' Z:',num2str(Gridshift(3))))
%% Spatial FFT for each channel
outputdata.fftfiddata=Phasecorrection(ifftshift(fftshift(ifft(fftn(squeeze(outputdata.avgrawdata)),[],1)),1));
outputdata.spectradata=fftshift(fft(outputdata.fftfiddata,[],1),1);
disp('Spatial FFT applied.')

%% Circular shift
%In RL direction

outputdata.fftfiddata=circshift(outputdata.fftfiddata,-1,3); disp('Circularshift in RL for fid.')
outputdata.spectradata=circshift(outputdata.spectradata,-1,3); disp('Circularshift in RL for spectra.')

%% First order phase correction with linear prediction SVD
tic
disp('Linear Prediction with SVD:');
if  isempty(outputdata.Parameters.Index.channelIndex)==1
    if outputdata.Parameters.NP > 512
        disp(strcat('Using first half of the FID=',num2str(outputdata.Parameters.NP/2),' time points.'));
        outputdata.LPfid=LinearPredSVD(outputdata.fftfiddata(1:outputdata.Parameters.NP/2,:,:,:),outputdata.Parameters);
        outputdata.LPfid((outputdata.Parameters.NP/2)+1:(outputdata.Parameters.NP),:,:,:)=outputdata.fftfiddata((outputdata.Parameters.NP/2)-outputdata.Parameters.missingpoints+1:outputdata.Parameters.NP-outputdata.Parameters.missingpoints,:,:,:);
    end
end
toc

%% Denoising, Spectral apodization and zerofilling

outputdata.LPfid=outputdata.Parameters.apodfunc.*outputdata.LPfid;
outputdata.LPfid(end:outputdata.Parameters.zerofill,:)=0;
disp(strcat('Spectral apodization(Lorentzian-',num2str(outputdata.Parameters.Apodizationparameter),' Hz) and Zero filling(x2) applied.'))

%% Denoising
outputdata.denoised=PCA_CSIdenoising(outputdata.LPfid,3,outputdata.Parameters); % Patch size is fixed to 3
disp('Denoising is applied to  in the final step.');

%% Final spectra
outputdata.LPspec=fftshift(fft(Phasecorrection(outputdata.LPfid),[],1),1);
%% SNR calculations
disp('Standard deviation of spectral noise (real) between 20 and 30 ppm downfield is used for noise estimation.')
outputdata.noisewindowzf=find(outputdata.xaxiszerofill>20 & outputdata.xaxiszerofill<30);
outputdata.waterwindowzf=find(outputdata.xaxiszerofill>4 & outputdata.xaxiszerofill<5.4);

outputdata.noisemap.LPspec=zeros(outputdata.Parameters.CSIdims);
outputdata.noisemap.LPspec=squeeze(std(real(outputdata.LPspec(outputdata.noisewindowzf,:,:,:))));
outputdata.SNR=squeeze(max(real(outputdata.LPspec(outputdata.waterwindowzf,:,:,:)),[],1))./outputdata.noisemap.LPspec;

% Check noise levels
% Create noise maps(zeros)
% outputdata.noisemap.FinalspectraDenoised=zeros(outputdata.Parameters.CSIdims);
% Estimate noise in noise window
% outputdata.noisemap.denoised=squeeze(std(real(outputdata.denoised(outputdata.noisewindowzf,:,:,:))));
% outputdata.SNR.Denoised=squeeze(max(real(outputdata.denoised(outputdata.waterwindowzf,:,:,:)),[],1))./outputdata.noisemap.denoised;

%% clear raw data to open space on the RAM
outputdata=rmfield(outputdata,'data');
outputdata=rmfield(outputdata,'avgrawdata');
%% Quantification
[outputdata.spectradata_aligned, outputdata.Referenceloc]=SpectralFrequencyAlignment(outputdata.spectradata,outputdata.xaxis,outputdata.Parameters,Referenceloc);
outputdata.AMARES_Results=BrainDMI_AMARES(outputdata.fftfiddata,outputdata.Parameters,outputdata.xaxis,outputdata.Referenceloc);
end