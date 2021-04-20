function outputdata=DMIPipeline_Singleloop(rawdata,TE)
%% Pipeline modified for Brain DMI data with single loop
Referenceloc=0;
SpatialFilter=0; %k-space weighting will be in Hamming
% rawdata='raw_020_3DCSI_baseline.data';
% TE=1.95;
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
outputdata.Parameters.TE=TE*10^-3; % Echo time as an input. Reverse it later !!!
disp(strcat('Acquisition echo time:', num2str(outputdata.Parameters.TE*10^3),'ms'))
outputdata.Parameters.missingpoints=ceil(outputdata.Parameters.BW*outputdata.Parameters.TE);
outputdata.Parameters.Apodizationparameter=1; % Line broadening parameter(in Hertz)
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

if SpatialFilter==0
    [outputdata.idealhammingwindow,outputdata.Correctionmask,outputdata.avgrawdata]=Hammingfilterpostcorrection(outputdata.avgrawdata,outputdata.acqpattern,outputdata.Parameters);
else
    [outputdata.idealGaussianwindow,outputdata.Correctionmask,outputdata.avgrawdata]=HammingtoGaussian(outputdata.avgrawdata,outputdata.acqpattern,outputdata.Parameters);
end

%% FFT for each channel
if  isempty(outputdata.Parameters.Index.channelIndex)==1
    outputdata.fftfiddata(:,:,:,:)=Phasecorrection(ifftshift(fftshift(ifft(fftn(squeeze(outputdata.avgrawdata(:,:,:,:))),[],1)),1));
    outputdata.spectradata(:,:,:,:)=fftshift(fft(outputdata.fftfiddata(:,:,:,:),[],1),1);
end
disp('Spatial FFT applied.')

%% Circular shift
%In RL direction
if isempty(outputdata.Parameters.Index.channelIndex)==0
    outputdata.RoemerEqualfid=circshift(outputdata.RoemerEqualfid,-1,3); disp('Circularshift in RL.')
    outputdata.combinedspectradata=circshift(outputdata.combinedspectradata,-1,3);
    outputdata.WSVD.WSVDcombined=circshift(outputdata.WSVD.WSVDcombined,-1,3);
end
%In FH direction
% outputdata.RoemerEqualfid=circshift(outputdata.RoemerEqualfid,1,4); disp('Circularshift in FH.')
% outputdata.combinedspectradata=circshift(outputdata.combinedspectradata,1,4);

%% SNR mask before linear prediction of missing points and apodization
if isempty(outputdata.Parameters.Index.channelIndex)==0
    
    outputdata.noisemap.combinedspectradata=zeros([outputdata.Parameters.CSIdims outputdata.Parameters.dims(end)]);
    outputdata.noisemap.combinedspectradata=squeeze(std(outputdata.combinedspectradata((7/8)*end:end,:,:,:,:)));
    outputdata.SNR.combinedspectradata=squeeze(max(real(outputdata.combinedspectradata),[],1))./outputdata.noisemap.combinedspectradata;
end

%% First order phase correction with linear prediction SVD
tic
disp('Linear Prediction with SVD:');
if  isempty(outputdata.Parameters.Index.channelIndex)==1
    outputdata.LPfid=LinearPredSVD(outputdata.fftfiddata,outputdata.Parameters);
    
else
    if outputdata.Parameters.NP > 512
        disp(strcat('Using first half of the FID=',num2str(outputdata.Parameters.NP/2),' time points.'));
        % Roemer equal noise channel combination
        outputdata.LPfidRoemer=LinearPredSVD(outputdata.RoemerEqualfid(1:outputdata.Parameters.NP/2,:,:,:),outputdata.Parameters);
        outputdata.LPfidRoemer((outputdata.Parameters.NP/2)+1:(outputdata.Parameters.NP),:,:,:)=outputdata.RoemerEqualfid((outputdata.Parameters.NP/2)-outputdata.Parameters.missingpoints+1:outputdata.Parameters.NP-outputdata.Parameters.missingpoints,:,:,:);
        disp('Linear Prediction of missing points with SVD applied for Roemer equal noise channel combined data.')
        
        % WSVD channel combination
        
        outputdata.LPfidWSVD=LinearPredSVD(outputdata.WSVD.WSVDcombined(1:outputdata.Parameters.NP/2,:,:,:),outputdata.Parameters);
        outputdata.LPfidWSVD((outputdata.Parameters.NP/2)+1:(outputdata.Parameters.NP),:,:,:)=outputdata.WSVD.WSVDcombined((outputdata.Parameters.NP/2)-outputdata.Parameters.missingpoints+1:outputdata.Parameters.NP-outputdata.Parameters.missingpoints,:,:,:);
        disp('Linear Prediction of missing points with SVD applied for WSVD channel combined data.')
        
    else
        outputdata.LPfidRoemer=LinearPredSVD(outputdata.RoemerEqualfid,outputdata.Parameters);disp('Using whole FID for linear prediction');
        
    end
    
    disp('Time domain data could be found under LPfidRoemer.')
    
end
toc

%% Denoising, Spectral apodization and zerofilling

% outputdata.denoisedafterLP=PCA_CSIdenoising(outputdata.LPfidRoemer,3,outputdata.Parameters); % Patch size is fixed to 3
if  isempty(outputdata.Parameters.Index.channelIndex)==1
    outputdata.LPfid=outputdata.Parameters.apodfunc.*outputdata.LPfid;
    outputdata.LPfid(end:outputdata.Parameters.zerofill,:)=0;
end

disp(strcat('Spectral apodization(Lorentzian-',num2str(outputdata.Parameters.Apodizationparameter),' Hz) and Zero filling(x2) applied.'))

%% Denoising
% if isempty(outputdata.Parameters.Index.channelIndex)==0
%     outputdata.denoised=PCA_CSIdenoising(outputdata.LPRoemerApod,3,outputdata.Parameters); % Patch size is fixed to 3
%     disp('Denoising is applied in the final step.')
% end
%% Final spectra

if  isempty(outputdata.Parameters.Index.channelIndex)==1
    outputdata.LPspec=fftshift(fft(Phasecorrection(outputdata.LPfid),[],1),1);
    [outputdata.LPspec, outputdata.Referenceloc]=SpectralFrequencyAlignment(outputdata.LPspec,outputdata.xaxiszerofill,outputdata.Parameters,Referenceloc);
end


%% SNR calculations
disp('Standard deviation of spectral noise (real) between 20 and 30 ppm downfield is used for noise estimation.')
outputdata.noisewindow=find(outputdata.xaxiszerofill>20 & outputdata.xaxiszerofill<30);
outputdata.waterwindow=find(outputdata.xaxiszerofill>4 & outputdata.xaxiszerofill<5.4);

if  isempty(outputdata.Parameters.Index.channelIndex)==1
    outputdata.noisemap.LPspec=zeros(outputdata.Parameters.CSIdims);
    outputdata.noisemap.LPspec=squeeze(std(real(outputdata.LPspec(outputdata.noisewindow,:,:,:))));
    outputdata.SNR=squeeze(max(real(outputdata.LPspec(outputdata.waterwindow,:,:,:)),[],1))./outputdata.noisemap.LPspec;
  
end

%% clear raw data to open space on the RAM
outputdata=rmfield(outputdata,'data');
% outputdata=rmfield(outputdata,'avgrawdata');
%% Quantification
% tic
% disp('HSVD based quantification:Roemer')
% [outputdata.Fittedtimesignal.Roemer, outputdata.Fittedspectra.Roemer, outputdata.quantifiedsignal.Roemer]=HSVD_quantification(Phasecorrection(outputdata.LPfidRoemer),outputdata.Parameters);
% disp('HSVD based quantification:WSVD')
% [outputdata.Fittedtimesignal.WSVD, outputdata.Fittedspectra.WSVD, outputdata.quantifiedsignal.WSVD]=HSVD_quantification(Phasecorrection(outputdata.LPfidWSVD),outputdata.Parameters);
%
% toc
