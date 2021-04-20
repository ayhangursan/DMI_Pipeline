function outputdata=DMIPipeline3D(rawdata)
%% Read Philips .data file
currentfolder=cd;
outputdata.datapath=strcat(currentfolder,'\',rawdata);
[outputdata.data, outputdata.list]=csi_loadData(outputdata.datapath);

%% Set parameters
outputdata.Parameters.gyromagneticratio=6.53569*10^6;
outputdata.Parameters.Tesla=7;
outputdata.Parameters.Freq=outputdata.Parameters.gyromagneticratio*outputdata.Parameters.Tesla;
outputdata.Parameters.BW=5000;
outputdata.Parameters.ppmwindow=((outputdata.Parameters.Freq/outputdata.Parameters.BW)^-1)*10^6;
outputdata.Parameters.NP=outputdata.list.F_resolution;
outputdata.Parameters.time= linspace(0,outputdata.Parameters.NP/outputdata.Parameters.BW,outputdata.Parameters.NP);
outputdata.Parameters.zerofill= outputdata.Parameters.NP*2;

outputdata.xaxis=linspace(-outputdata.Parameters.ppmwindow/2,outputdata.Parameters.ppmwindow/2,outputdata.Parameters.NP);
outputdata.xaxiszerofill=linspace(-outputdata.Parameters.ppmwindow/2,outputdata.Parameters.ppmwindow/2,outputdata.Parameters.zerofill);
outputdata.Parameters.TE=1.95*10^-3;
disp(strcat('Acquisition echo time:', num2str(outputdata.Parameters.TE*10^3),'ms'))
outputdata.Parameters.missingpoints=floor(outputdata.Parameters.BW*outputdata.Parameters.TE);

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
% outputdata.NCov=cov(outputdata.data.noise);disp('Noise covariance matrix from noise scan.') % Fetch noise data from noise scans

%% Average data
outputdata.avgrawdata=flip(squeeze(mean(outputdata.data.raw,outputdata.Parameters.Index.averageIndex)),outputdata.Parameters.Index.kyIndex); % If state needed for 1 NSA case
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
[outputdata.idealhammingwindow,outputdata.Correctionmask,outputdata.avgrawdata]=Hammingfilterpostcorrection(outputdata.avgrawdata,outputdata.acqpattern,outputdata.Parameters);disp('Hamming filter correction applied.')

% outputdata.avgrawdata=Hammingfilter3D(outputdata.avgrawdata,outputdata.Parameters);disp('Hamming filter applied.')

%% FFT for each channel
if  isempty(outputdata.Parameters.Index.channelIndex)==1
    outputdata.fftfiddata(:,:,:,:)=Phasecorrection(ifftshift(fftshift(ifft(fftn(squeeze(outputdata.avgrawdata(:,:,:,:))),[],1)),1));
    outputdata.spectradata(:,:,:,:)=flip(fftshift(fft(outputdata.fftfiddata(:,:,:,:),[],1),1),1);
else
    for n=1:size(outputdata.avgrawdata,outputdata.Parameters.Index.channelIndex)
        outputdata.fftfiddata(:,n,:,:,:)=Phasecorrection(ifftshift(fftshift(ifft(fftn(squeeze(outputdata.avgrawdata(:,n,:,:,:))),[],1)),1));
        outputdata.spectradata(:,n,:,:,:)=flip(fftshift(fft(outputdata.fftfiddata(:,n,:,:,:),[],1),1),1);
    end
end
disp('Spatial FFT applied.')
%% Channel combination
if isempty(outputdata.Parameters.Index.channelIndex)==0
    outputdata.NCov=cov(outputdata.fftfiddata((3/4)*end:end,:,1));disp(strcat('Noise covariance matrix from first fid-Last quarter:',num2str(outputdata.Parameters.NP/4))) % Calculate noise covariance matrix from first voxel
    [outputdata.RoemerEqualfid, outputdata.RoemerEqNoiseSensitivity_maps]=RoemerEqualNoise(outputdata.fftfiddata,outputdata.NCov,outputdata.Parameters.Index.channelIndex);
    %     [outputdata.WSVD.WSVDcombined, outputdata.WSVD.wsvdQuality, outputdata.WSVD.wsvdWeights]=WSVDcomb(outputdata.fftfiddata,outputdata.NCov,outputdata.Parameters.Index.channelIndex);
    outputdata.combinedspectradata=flip(fftshift(fft(outputdata.RoemerEqualfid,[],1),1),1); % To be used in SNR mask for linear prediction
    
    disp('Roemer equal noise combinations applied.')
end

%% SNR mask before linear prediction of missing points and apodization
if  isempty(outputdata.Parameters.Index.channelIndex)==0
    
    outputdata.noisemap.combinedspectradata=zeros([outputdata.Parameters.CSIdims outputdata.Parameters.dims(end)]);
    outputdata.noisemap.combinedspectradata=squeeze(std(outputdata.combinedspectradata((15/16)*end:end,:,:,:)));
    outputdata.SNR.combinedspectradata=squeeze(max(real(outputdata.combinedspectradata),[],1))./outputdata.noisemap.combinedspectradata;
    
    outputdata.SNRmask.combined=zeros(outputdata.Parameters.CSIdims);
    outputdata.SNRmask.combined(find(outputdata.SNR.combinedspectradata>max(outputdata.SNR.combinedspectradata,[],'all')*0.05))=1;
    
    outputdata.RoemerEqualfidmasked=NaN(size(outputdata.SNRmask.combined));
    outputdata.RoemerEqualfidmasked(find(outputdata.SNRmask.combined))=outputdata.RoemerEqualfid(find(outputdata.SNRmask.combined));
end

%% First order phase correction with linear prediction SVD
tic
disp('Linear Prediction with SVD:');
if  isempty(outputdata.Parameters.Index.channelIndex)==1
    outputdata.LPfid=LinearPredSVD(outputdata.fftfiddata,outputdata.Parameters);
    
else
    outputdata.LPfidRoemer=LinearPredSVD(outputdata.RoemerEqualfid,outputdata.Parameters);
    %     outputdata.LPfidWSVD=LinearPredSVD(outputdata.WSVD.WSVDcombined,outputdata.Parameters);
    
end
toc
disp('Linear Prediction of missing points with SVD applied.')

%% Spectral apodization and zerofilling
outputdata.Parameters.apodfunc=exp(-20*pi*outputdata.Parameters.time).'; %20 Hz Lorentzian apodization
if  isempty(outputdata.Parameters.Index.channelIndex)==1
    outputdata.LPfid=outputdata.Parameters.apodfunc.*outputdata.LPfid;
    outputdata.LPfid(end:outputdata.Parameters.zerofill,:)=0;
else
    outputdata.LPfidRoemer_apod=outputdata.Parameters.apodfunc.*outputdata.LPfidRoemer;
    outputdata.LPfidRoemer_apod(end:outputdata.Parameters.zerofill,:)=0;
    
    %     outputdata.LPfidWSVD=outputdata.Parameters.apodfunc.*outputdata.LPfidWSVD;
    %     outputdata.LPfidWSVD(end:outputdata.Parameters.zerofill,:)=0;
end
disp('Spectral apodization(20 Hz) and Zero filling(x2) applied.')

%% Final spectra

if  isempty(outputdata.Parameters.Index.channelIndex)==1
    outputdata.LPspec=flip(fftshift(fft(Phasecorrection(outputdata.LPfid),[],1),1),1);
else
    outputdata.FinalspectraRoemer=flip(fftshift(fft(Phasecorrection(outputdata.LPfidRoemer_apod),[],1),1),1);
    %     outputdata.FinalspectraWSVD=flip(fftshift(fft(Phasecorrection(outputdata.LPfidWSVD),[],1),1),1);
end
%% SNR calculation
if  isempty(outputdata.Parameters.Index.channelIndex)==1
    
    outputdata.noisemap.LPspec=zeros(outputdata.Parameters.CSIdims);
    outputdata.noisemap.LPspec=squeeze(std(outputdata.LPspec((15/16)*end:end,:,:,:)));
    outputdata.SNR=squeeze(max(real(outputdata.LPspec),[],1))./outputdata.noisemap.LPspec;
    
    snrcoln=4;snrrown=1;snrslicen=4;
    colorbarlimits=[min([outputdata.SNR],[],'all') max([outputdata.SNR],[],'all')];
    
    figure;
    for m=1:snrslicen
        subplot(snrrown,snrcoln,m);
        slice=floor(outputdata.Parameters.CSIdims(3)/2)-floor(snrslicen/2)+m;
        imagesc(outputdata.SNR(:,:,slice));caxis(colorbarlimits);colorbar;daspect([outputdata.Parameters.CSIdims(1) outputdata.Parameters.CSIdims(2) 1]);
        title(strcat('SNR Slice:',num2str(slice)))
    end
else
    
    outputdata.noisemap.FinalspectraRoemer=zeros(outputdata.Parameters.CSIdims);
    %     outputdata.noisemap.FinalspectraWSVD=zeros(outputdata.Parameters.CSIdims);
    
    outputdata.noisemap.FinalspectraRoemer=squeeze(std(outputdata.FinalspectraRoemer((15/16)*end:end,:,:,:)));
    %     outputdata.noisemap.FinalspectraWSVD=squeeze(std(outputdata.FinalspectraWSVD((15/16)*end:end,:,:,:)));
    
    outputdata.SNR.Roemer=squeeze(max(real(outputdata.FinalspectraRoemer),[],1))./outputdata.noisemap.FinalspectraRoemer;
    %     outputdata.SNR.WSVD=squeeze(max(real(outputdata.FinalspectraWSVD),[],1))./outputdata.noisemap.FinalspectraWSVD;
    snrcoln=4;snrrown=2;snrslicen=4;
    %     colorbarlimits=[min([outputdata.SNR.Roemer outputdata.SNR.WSVD],[],'all') max([outputdata.SNR.Roemer outputdata.SNR.WSVD],[],'all')];
    
    %     figure;
    %     for m=1:snrslicen
    %         subplot(snrrown,snrcoln,m);
    %         slice=floor(outputdata.Parameters.CSIdims(3)/2)-floor(snrslicen/2)+m;
    %         imagesc(outputdata.SNR.Roemer(:,:,slice));caxis(colorbarlimits);colorbar;daspect([outputdata.Parameters.CSIdims(1) outputdata.Parameters.CSIdims(2) 1]);
    %         title(strcat('SNR Roemer Slice:',num2str(slice)))
    %         subplot(snrrown,snrcoln,m+snrcoln);
    %         imagesc(outputdata.SNR.WSVD(:,:,slice));caxis(colorbarlimits);colorbar;daspect([outputdata.Parameters.CSIdims(1) outputdata.Parameters.CSIdims(2) 1]);
    %         title(strcat('SNR WSVD Slice:',num2str(slice)))
    %     end
    outputdata.SNRmask.FinalSpectra=find(outputdata.SNR.Roemer(:,:,:,1)>max(outputdata.SNR.Roemer(:,:,:,1),[],'all')*0.1);
    disp('SNR images and mask created from final spectra.')
    
end



%% Quantification
% disp('Quantification step is starting.')
% tic
% if  isempty(outputdata.Parameters.Index.channelIndex)==1
%     [outputdata.metabmap, outputdata.metabfits, outputdata.FWHMmap]=Lorentzianfitgen(outputdata.LPspec,outputdata.xaxiszerofill, outputdata.Parameters.BW);
% else
%     [outputdata.metabmap.Roemer, outputdata.metabfits.Roemer, outputdata.FWHMmap.Roemer]=Lorentzianfitgen(outputdata.FinalspectraRoemer,outputdata.xaxiszerofill, outputdata.Parameters.BW);
%     [outputdata.metabmap.WSVD, outputdata.metabfits.WSVD, outputdata.FWHMmap.WSVD]=Lorentzianfitgen(outputdata.FinalspectraWSVD,outputdata.xaxiszerofill, outputdata.Parameters.BW);
% end
% toc
% disp('Quantification step is finished.')
end