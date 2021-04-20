function outputdata=DMIPipeline3D_refscan(rawdata,refscan,refNCov,TE,Gridshift,Referenceloc)
SpatialFilter=0; %k-space weighting is done in Hamming filter fashion

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
outputdata.Parameters.FirstOrdPhaseFunct=exp(-1i * (2* pi * (outputdata.xaxis-4.7).'*(outputdata.Parameters.Freq/(10^6)) * outputdata.Parameters.TE));

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
%% Apply half-voxel shift in FH
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
    if refNCov==0
        outputdata.NCov=cov(outputdata.fftfiddata((7/8)*end:end,:,1));disp('Noise Covariance matrix generated'); % Calculate noise covariance matrix from first voxel. From August 2020 we acquire a noise scan before CSI scans to generate noise covariance matrix.
        [outputdata.RoemerEqualfid, outputdata.RoemerEqNoiseSensitivity_maps]=RoemerEqualNoise(outputdata.fftfiddata,outputdata.NCov,outputdata.Parameters.Index.channelIndex);
        outputdata.combinedspectradata=fftshift(fft(outputdata.RoemerEqualfid,[],1),1); % To be used in SNR mask for linear prediction
        
    else
        disp('Input Noise Covariance matrix is used');
        [outputdata.RoemerEqualfid, outputdata.RoemerEqNoiseSensitivity_maps]=RoemerEqualNoise_withmap_input(outputdata.fftfiddata,refscan,refNCov,outputdata.Parameters.Index.channelIndex);        disp('Roemer equeal noise channel combination applied.')
        [outputdata.WSVD.WSVDcombined, outputdata.WSVD.wsvdQuality, outputdata.WSVD.wsvdWeights]=WSVDcomb(outputdata.fftfiddata,refNCov,outputdata.Parameters.Index.channelIndex);        disp('WSVD channel combination applied.')
    end
end
%% Circular shift
%In RL direction
if isempty(outputdata.Parameters.Index.channelIndex)==0
    
    outputdata.RoemerEqualfid=circshift(outputdata.RoemerEqualfid,-1,3); disp('Circularshift in RL for Roemer comb(fid).')
    outputdata.WSVD.WSVDcombined=circshift(outputdata.WSVD.WSVDcombined,-1,3); disp('Circularshift in RL for WSVD.')
end
outputdata.combinedspectradata=fftshift(fft(outputdata.RoemerEqualfid,[],1),1); % To be used in SNR mask for linear prediction
outputdata.WSVDcombinedspectradata=fftshift(fft(outputdata.WSVD.WSVDcombined,[],1),1); % To be used in SNR mask for linear prediction

%% SNR mask before linear prediction of missing points and apodization
if isempty(outputdata.Parameters.Index.channelIndex)==0
    outputdata.noisewindow=find(outputdata.xaxis>20 & outputdata.xaxis<30);
    outputdata.waterwindow=find(outputdata.xaxis>4 & outputdata.xaxis<5.4);
    outputdata.noisemap.combinedspectradata=zeros([outputdata.Parameters.CSIdims]);
    outputdata.noisemap.combinedspectradata=squeeze(std(outputdata.combinedspectradata(outputdata.noisewindow,:,:,:)));
    outputdata.noisemap.WSVDcombinedspectradata=squeeze(std(outputdata.WSVDcombinedspectradata(outputdata.noisewindow,:,:,:)));
    
    outputdata.SNR.combinedspectradata=squeeze(max(real(outputdata.combinedspectradata(outputdata.waterwindow,:,:,:)),[],1))./outputdata.noisemap.combinedspectradata;
    outputdata.SNR.WSVDcombinedspectradata=squeeze(max(real(outputdata.WSVDcombinedspectradata(outputdata.waterwindow,:,:,:)),[],1))./outputdata.noisemap.WSVDcombinedspectradata;
    
    %     outputdata.SNRmask.combined=zeros(outputdata.Parameters.CSIdims);
    %     outputdata.SNRmask.combined(find(outputdata.SNR.combinedspectradata(:,:,:,1)>max(outputdata.SNR.combinedspectradata(:,:,:,1),[],'all')*0.05))=1;
    %     outputdata.SNRmask.combinedspectral=permute(repmat(outputdata.SNRmask.combined,[1 1 1 outputdata.Parameters.NP outputdata.Parameters.dyn]),[4 1 2 3 5]);
    %
    %     outputdata.RoemerEqualfidmasked=NaN(size(outputdata.SNRmask.combinedspectral));
    %     outputdata.RoemerEqualfidmasked(find(outputdata.SNRmask.combinedspectral))=outputdata.RoemerEqualfid(find(outputdata.SNRmask.combinedspectral));
    %
    %
    % end
    
    %% First order phase correction with linear prediction SVD
    tic
    disp('Linear Prediction with SVD:');
    if  isempty(outputdata.Parameters.Index.channelIndex)==1
        if outputdata.Parameters.NP > 512
            disp(strcat('Using first half of the FID=',num2str(outputdata.Parameters.NP/2),' time points.'));
            outputdata.LPfid=LinearPredSVD(outputdata.fftfiddata(1:outputdata.Parameters.NP/2,:,:,:),outputdata.Parameters);
            outputdata.LPfid((outputdata.Parameters.NP/2)+1:(outputdata.Parameters.NP),:,:,:)=outputdata.fftfiddata((outputdata.Parameters.NP/2)-outputdata.Parameters.missingpoints+1:outputdata.Parameters.NP-outputdata.Parameters.missingpoints,:,:,:);
        end
        
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
            %
        else
            outputdata.LPfidRoemer=LinearPredSVD(outputdata.RoemerEqualfid,outputdata.Parameters);disp('Using whole FID for linear prediction');
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
    else
        
        outputdata.LPRoemerApod=outputdata.Parameters.apodfunc.*outputdata.LPfidRoemer;
        outputdata.LPRoemerApod(end:outputdata.Parameters.zerofill,:)=0;
        
        outputdata.LPWSVDApod=outputdata.Parameters.apodfunc.*outputdata.LPfidWSVD;
        outputdata.LPWSVDApod(end:outputdata.Parameters.zerofill,:)=0;
    end
    disp(strcat('Spectral apodization(Lorentzian-',num2str(outputdata.Parameters.Apodizationparameter),' Hz) and Zero filling(x2) applied.'))
    
    %% Denoising
    if  isempty(outputdata.Parameters.Index.channelIndex)==1
        outputdata.denoised=PCA_CSIdenoising(outputdata.LPfid,3,outputdata.Parameters); % Patch size is fixed to 3
        
    else
        outputdata.denoisedRoemer=PCA_CSIdenoising_V2(outputdata.combinedspectradata,5,outputdata.Parameters).*outputdata.Parameters.FirstOrdPhaseFunct; % Patch size is fixed to 3
        outputdata.denoisedWSVD=PCA_CSIdenoising_V2(outputdata.WSVDcombinedspectradata,5,outputdata.Parameters).*outputdata.Parameters.FirstOrdPhaseFunct; % Patch size is fixed to 3
        
    end
    disp('Denoising is applied to  in the final step.')
    
    %% Final spectra
    
    if  isempty(outputdata.Parameters.Index.channelIndex)==1
        outputdata.LPspec=fftshift(fft(Phasecorrection(outputdata.LPfid),[],1),1);
        %     [outputdata.LPspec, outputdata.Referenceloc]=SpectralFrequencyAlignment(outputdata.LPspec,outputdata.xaxiszerofill,outputdata.Parameters,Referenceloc);
        
    else
        outputdata.FinalspectraRoemer=fftshift(fft(Phasecorrection(outputdata.LPRoemerApod),[],1),1);
        outputdata.FinalspectraDenoised=fftshift(fft(Phasecorrection(outputdata.denoisedRoemer),[],1),1);
        %     [outputdata.FinalspectraRoemer_aligned, outputdata.Referenceloc]=SpectralFrequencyAlignment(outputdata.FinalspectraRoemer,outputdata.xaxiszerofill,outputdata.Parameters,Referenceloc);
        %     [outputdata.FinalspectraDenoised_aligned, outputdata.Referenceloc]=SpectralFrequencyAlignment(outputdata.FinalspectraDenoised,outputdata.xaxiszerofill,outputdata.Parameters,Referenceloc);
        
        outputdata.FinalspectraWSVD=fftshift(fft(Phasecorrection(outputdata.LPWSVDApod),[],1),1);
    end
    
    
    %% SNR calculations
    disp('Standard deviation of spectral noise (real) between 20 and 30 ppm downfield is used for noise estimation.')
    outputdata.noisewindowzf=find(outputdata.xaxiszerofill>20 & outputdata.xaxiszerofill<30);
    outputdata.waterwindowzf=find(outputdata.xaxiszerofill>4 & outputdata.xaxiszerofill<5.4);
    
    if  isempty(outputdata.Parameters.Index.channelIndex)==1
        outputdata.noisemap.LPspec=zeros(outputdata.Parameters.CSIdims);
        outputdata.noisemap.LPspec=squeeze(std(real(outputdata.LPspec(outputdata.noisewindow,:,:,:))));
        outputdata.SNR=squeeze(max(real(outputdata.LPspec(outputdata.waterwindow,:,:,:)),[],1))./outputdata.noisemap.LPspec;
        
    else
        
        % Check noise levels
        % Create noise maps(zeros)
        outputdata.noisemap.FinalspectraRoemer=zeros(outputdata.Parameters.CSIdims);
        outputdata.noisemap.FinalspectraWSVD=zeros(outputdata.Parameters.CSIdims);
        outputdata.noisemap.FinalspectraDenoised=zeros(outputdata.Parameters.CSIdims);
        % Estimate noise in noise window
        outputdata.noisemap.FinalspectraRoemer=squeeze(std(real(outputdata.FinalspectraRoemer(outputdata.noisewindowzf,:,:,:)))); % Using real part of the spectra to estimate noise in accordance with MRS consensus paper(doi:10.1002/nbm.4347)
        outputdata.noisemap.FinalspectraWSVD=squeeze(std(real(outputdata.FinalspectraWSVD(outputdata.noisewindowzf,:,:,:)))); % Using real part of the spectra to estimate noise in accordance with MRS consensus paper(doi:10.1002/nbm.4347)
        outputdata.noisemap.FinalspectraDenoised=squeeze(std(real(outputdata.FinalspectraDenoised(outputdata.noisewindow,:,:,:))));
        %% Water based SNR
        
        outputdata.SNR.Roemer=squeeze(max(real(outputdata.FinalspectraRoemer(outputdata.waterwindowzf,:,:,:)),[],1))./outputdata.noisemap.FinalspectraRoemer;
        outputdata.SNR.WSVD=squeeze(max(real(outputdata.FinalspectraWSVD(outputdata.waterwindowzf,:,:,:)),[],1))./outputdata.noisemap.FinalspectraWSVD;
        outputdata.SNR.Denoised=squeeze(max(real(outputdata.FinalspectraDenoised(outputdata.waterwindow,:,:,:)),[],1))./outputdata.noisemap.FinalspectraDenoised;
        
        % Generating SNR mask
        outputdata.SNRmask.Roemer=find(outputdata.SNR.Roemer(:,:,:,1)>max(outputdata.SNR.Roemer(:,:,:,1),[],'all')*0.1);
        outputdata.SNRmask.Denoised=find(outputdata.SNR.Denoised(:,:,:,1)>max(outputdata.SNR.Roemer(:,:,:,1),[],'all')*0.1);
        
    end
    [outputdata.combinedspectradata_aligned, outputdata.Referenceloc]=SpectralFrequencyAlignment(outputdata.combinedspectradata,outputdata.xaxis,outputdata.Parameters,Referenceloc);
    
    %% clear raw data to open space on the RAM
    outputdata=rmfield(outputdata,'data');
%     outputdata=rmfield(outputdata,'avgrawdata');
    %% Quantification
    % if  isempty(outputdata.Parameters.Index.channelIndex)==1
    %     outputdata.AMARES=DMI_AMARES(outputdata.fftfiddata,outputdata.Parameters,outputdata.xaxis);
    %
    % else
    % % tic
    % % disp('HSVD based quantification:Roemer')
    % % [outputdata.Fittedtimesignal.Roemer, outputdata.Fittedspectra.Roemer, outputdata.quantifiedsignal.Roemer]=HSVD_quantification(Phasecorrection(outputdata.LPfidRoemer),outputdata.Parameters);
    % % disp('HSVD based quantification:WSVD')
    % % [outputdata.Fittedtimesignal.WSVD, outputdata.Fittedspectra.WSVD, outputdata.quantifiedsignal.WSVD]=HSVD_quantification(Phasecorrection(outputdata.LPfidWSVD),outputdata.Parameters);
    % %
    % % toc
    
    %
%     outputdata.AMARES_Results=DMI_AMARES(outputdata.RoemerEqualfid,outputdata.Parameters,outputdata.xaxis,outputdata.Referenceloc);
    % end
end