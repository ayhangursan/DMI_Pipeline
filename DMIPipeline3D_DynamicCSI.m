function outputdata=DMIPipeline3D_DynamicCSI(rawdata,noisescan,TE)
disp(strcat('Processing CSI data:',rawdata));
disp(strcat('Noise scan:',noisescan));
%Read seperate noise scan
outputdata.NoiseScan=DMINoiseCovariance(noisescan);
outputdata.NoiseCov=outputdata.NoiseScan.noisecovariance;
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
outputdata.Parameters.dyn=outputdata.list.number_of_extra_attribute_1_values;

outputdata.xaxis=linspace(-outputdata.Parameters.ppmwindow/2,outputdata.Parameters.ppmwindow/2,outputdata.Parameters.NP)+4.7;
outputdata.xaxiszerofill=linspace(-outputdata.Parameters.ppmwindow/2,outputdata.Parameters.ppmwindow/2,outputdata.Parameters.zerofill)+4.7;
outputdata.Parameters.TE=TE*10^-3;
disp(strcat('Acquisition echo time:', num2str(outputdata.Parameters.TE*10^3),' ms'))
outputdata.Parameters.missingpoints=ceil(outputdata.Parameters.BW*outputdata.Parameters.TE);
outputdata.Parameters.Apodizationparameter=5; % Line broadening parameter(in Hertz)
% Index
outputdata.Parameters.apodfunc=exp(-pi*outputdata.Parameters.Apodizationparameter.*outputdata.Parameters.time).'; %Lorentzian apodization. This value could be decreased as linewidths are a bit broad now.% Index
outputdata.Parameters.Index.dynIndex = find(contains(outputdata.data.labels,'extr1'));
outputdata.Parameters.Index.channelIndex = find(contains(outputdata.data.labels,'chan'));
outputdata.Parameters.Index.averageIndex = find(contains(outputdata.data.labels,'aver'));
outputdata.Parameters.Index.kxIndex = find(contains(outputdata.data.labels,'kx'));
outputdata.Parameters.Index.kyIndex = find(contains(outputdata.data.labels,'ky'));
outputdata.Parameters.Index.kzIndex = find(contains(outputdata.data.labels,'kz'));
outputdata.Parameters.dims=size(outputdata.data.raw);
outputdata.Parameters.CSIdims=[outputdata.Parameters.dims(outputdata.Parameters.Index.kxIndex) outputdata.Parameters.dims(outputdata.Parameters.Index.kyIndex) outputdata.Parameters.dims(outputdata.Parameters.Index.kzIndex)];
outputdata.Parameters.FirstOrdPhaseFunct=exp(-1i * (2* pi * (outputdata.xaxis-4.7).'*(outputdata.Parameters.Freq/(10^6)) * outputdata.Parameters.TE));


%% Average data
outputdata.avgrawdata=flip(squeeze(mean(outputdata.data.raw,outputdata.Parameters.Index.averageIndex)),outputdata.Parameters.Index.kyIndex); % If state needed for 1 NSA case
dataorder=1:ndims(outputdata.avgrawdata);
% dataorder([outputdata.Parameters.Index.kyIndex outputdata.Parameters.Index.kxIndex])=dataorder([outputdata.Parameters.Index.kxIndex outputdata.Parameters.Index.kyIndex]);
outputdata.avgrawdata=permute(outputdata.avgrawdata,dataorder);
outputdata.fftfiddata=zeros(size(outputdata.avgrawdata));
outputdata.rawspectradata=zeros(size(outputdata.avgrawdata));
%% Spatial filter-Hamming for 3D
if isempty(outputdata.Parameters.Index.channelIndex)==0
    outputdata.acqpattern=acquisitionpatterncheckSeries(outputdata.data.raw);
end
outputdata.normalizedacqpattern=outputdata.acqpattern/outputdata.Parameters.dims(outputdata.Parameters.Index.averageIndex);
[outputdata.idealhammingwindow,outputdata.Correctionmask,outputdata.avgrawdata]=HammingfilterpostcorrectionSeries(outputdata.avgrawdata,outputdata.normalizedacqpattern,outputdata.Parameters);disp('Hamming filter correction applied.')

%% FFT for each channel
for dyn=1:outputdata.Parameters.dims(outputdata.Parameters.Index.dynIndex)
    if  isempty(outputdata.Parameters.Index.channelIndex)==1
        outputdata.fftfiddata(:,:,:,:,dyn)=Phasecorrection(ifftshift(fftshift(ifft(fftn(squeeze(outputdata.avgrawdata(:,:,:,:,dyn))),[],1)),1));
        outputdata.rawspectradata(:,:,:,:,dyn)=flip(fftshift(fft(outputdata.fftfiddata(:,:,:,:,dyn),[],1),1),1);
    else
        for n=1:size(outputdata.avgrawdata,outputdata.Parameters.Index.channelIndex)
            outputdata.fftfiddata(:,n,:,:,:,dyn)=Phasecorrection(ifftshift(fftshift(ifft(fftn(squeeze(outputdata.avgrawdata(:,n,:,:,:,dyn))),[],1)),1));
            outputdata.rawspectradata(:,n,:,:,:,dyn)=flip(fftshift(fft(outputdata.fftfiddata(:,n,:,:,:,dyn),[],1),1),1);
        end
    end
end
disp('Spatial FFT applied.')
%% Channel combination
if isempty(outputdata.Parameters.Index.channelIndex)==0
    outputdata.RoemerEqualfid=zeros([outputdata.Parameters.NP outputdata.Parameters.CSIdims outputdata.Parameters.dims(end)]);
    [outputdata.RoemerEqualfid(:,:,:,:,1), outputdata.RoemerEqNoiseSensitivity_maps]=RoemerEqualNoise(squeeze(outputdata.fftfiddata(:,:,:,:,:,1)),outputdata.NoiseCov,outputdata.Parameters.Index.channelIndex);
    for m=2:outputdata.Parameters.dims(end)
        disp(strcat('Scan dyn:',num2str(m)))
        [outputdata.RoemerEqualfid(:,:,:,:,m), ~]=RoemerEqualNoise_withmap_input(squeeze(outputdata.fftfiddata(:,:,:,:,:,m)),outputdata.RoemerEqNoiseSensitivity_maps,outputdata.NoiseCov,outputdata.Parameters.Index.channelIndex);
        
    end
    
    outputdata.combinedspectradata=fftshift(fft(outputdata.RoemerEqualfid,[],1),1).*outputdata.Parameters.FirstOrdPhaseFunct; % First order phase correction is incorporated in this line.
    disp('Roemer equeal noise channel combination applied.')
end

% %% First order phase correction with linear prediction SVD
% tic
% disp('Linear Prediction with SVD:');
%
%     outputdata.LPfidRoemer=LinearPredSVD(outputdata.RoemerEqualfid,outputdata.Parameters);
%
% toc
% disp('Linear Prediction of missing points with SVD applied.')
% %% Better visualization
% Spectral apodization and zerofilling
% outputdata.Parameters.apodfunc=exp(-5*outputdata.Parameters.time).'; %20 Hz
% outputdata.LPRoemerApod=outputdata.Parameters.apodfunc.*outputdata.LPfidRoemer;
% outputdata.LPRoemerApod(end:outputdata.Parameters.zerofill,:)=0;
% disp('Spectral apodization(20 Hz) and Zero filling(x2) applied.')
% %% Final spectra
%
% outputdata.FinalspectraRoemer=flip(fftshift(fft(Phasecorrection(outputdata.LPfidRoemer),[],1),1),1);
%
%% Water based SNR
% Estimate noise in noise window
outputdata.noisewindowzf=find(outputdata.xaxis>20 & outputdata.xaxis<30);
outputdata.waterwindowzf=find(outputdata.xaxis>4 & outputdata.xaxis<5.4);
outputdata.noisemap.CombinedSpec=squeeze(std(real(outputdata.combinedspectradata(outputdata.noisewindowzf,:,:,:,:)))); % Using real part of the spectra to estimate noise in accordance with MRS consensus paper(doi:10.1002/nbm.4347)

outputdata.SNR.Roemer=squeeze(max(real(outputdata.combinedspectradata(outputdata.waterwindowzf,:,:,:,:)),[],1))./outputdata.noisemap.CombinedSpec;

%% Linewidths

outputdata.T2starmap=zeros(size(outputdata.SNR.Roemer));
outputdata.Hzmap=zeros(size(outputdata.SNR.Roemer));

for dyn=1:outputdata.Parameters.dims(outputdata.Parameters.Index.dynIndex)
[outputdata.T2starmap(:,:,:,dyn), outputdata.Hzmap(:,:,:,dyn), ~]=SpectralLinewidth(outputdata.combinedspectradata,outputdata.xaxis,outputdata.Parameters);
end
end

% %%
% AP=6;
% RL=3;
% FH=7;
% figure('units','normalized','outerposition',[0 0 1 1],'WindowState', 'maximized');
% fig1=waterfall(outputdata.xaxis,[1:11]*2,real(squeeze(outputdata.combinedspectradata(:,AP,RL,FH,:)).*exp(-1i * (2* pi * (outputdata.xaxis-4.7).'*45.7 * 0.00195))).');
% axes1 = fig1.Parent;
% set(gca,'XDir','reverse','ZtickLabel',[],'LineWidth',2,'Colormap',[0 0 0]);
% ax=gca;
% ax.ZAxis.Color='b';
% ylabel('Time (min)','FontSize',28,'Rotation',40);
% xlabel('^2H frequency (ppm)','HorizontalAlignment','right','FontSize',28);
% fig1.LineWidth=1.8;
% xlim([0 8]); % ppm axis limits
% view([24.1875 43.3418]);
% set(axes1,'FontSize',24,'FontWeight','bold','XDir','reverse','ZTick',[],'YGrid','off','XGrid','off','ZColor',[1 1 1],'Color',[1 1 1]);
% %%
% AP=6;
% RL=3;
% FH=6;
% figure('WindowState', 'maximized')
% title('HDO peak height','FontSize',28)
% HDOintensity=max(real(squeeze(outputdata.combinedspectradata(:,AP,RL,FH,:)).*exp(-1i * (2* pi * (outputdata.xaxis-4.7).'*45.7 * 0.00195))),[],1)./max(real(squeeze(outputdata.combinedspectradata(:,AP,RL,FH,1)).*exp(-1i * (2* pi * (outputdata.xaxis-4.7).'*45.7 * 0.00195))),[],1);
% plot([1:11]*2,HDOintensity,'-o','LineWidth',1.9)
% ylim([0 1.1])
% xlim([0 23])
% xlabel('Time (min)','FontSize',28);
% ylabel('Peak intensity (a.u)','FontSize',28);
% ax=gca;
% ax.XAxis.FontSize=24;
% ax.YAxis.FontSize=24;
