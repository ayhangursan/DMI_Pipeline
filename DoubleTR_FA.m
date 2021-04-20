function [FA_map, TRshort, TRlong, Parameters]=DoubleTR_FA(TRs,TRsdata,TRl,TRldata)
%% Read Philips .data file
currentfolder=cd;
% load data for short TR
TRshort.datapath=strcat(currentfolder,'\',TRsdata);
[TRshort.data, TRshort.list]=csi_loadData(TRshort.datapath);
% load data for short TR
TRlong.datapath=strcat(currentfolder,'\',TRldata);
[TRlong.data, TRlong.list]=csi_loadData(TRlong.datapath);

%% Parameters
% TR1>TR2
Parameters.gyromagneticratio=6.53569*10^6;
Parameters.Tesla=7;
Parameters.Freq=Parameters.gyromagneticratio*Parameters.Tesla;
Parameters.BW=5000;
Parameters.ppmwindow=((Parameters.Freq/Parameters.BW)^-1)*10^6;
Parameters.NP=TRshort.list.F_resolution;
Parameters.time= linspace(0,Parameters.NP/Parameters.BW,Parameters.NP);
Parameters.zerofill= Parameters.NP*2;

Parameters.xaxis=linspace(-Parameters.ppmwindow/2,Parameters.ppmwindow/2,Parameters.NP);
Parameters.xaxiszerofill=linspace(-Parameters.ppmwindow/2,Parameters.ppmwindow/2,Parameters.zerofill);
Parameters.TE=1.98*10^-3;
Parameters.missingpoints=floor(Parameters.BW*Parameters.TE);
% Index
Parameters.Index.channelIndex = find(contains(TRshort.data.labels,'chan'));
Parameters.Index.averageIndex = find(contains(TRshort.data.labels,'aver'));
Parameters.Index.kxIndex = find(contains(TRshort.data.labels,'kx'));
Parameters.Index.kyIndex = find(contains(TRshort.data.labels,'ky'));
Parameters.Index.kzIndex = find(contains(TRshort.data.labels,'kz'));
Parameters.dims=size(TRshort.data.raw);
Parameters.CSIdims=[Parameters.dims(Parameters.Index.kxIndex) Parameters.dims(Parameters.Index.kyIndex) Parameters.dims(Parameters.Index.kzIndex)];


Parameters.T1=350; %Water  T1 value (ms)
%% Average data
TRshort.avgrawdata=squeeze(mean(TRshort.data.raw,Parameters.Index.averageIndex)); % If state needed for 1 NSA case
TRlong.avgrawdata=squeeze(mean(TRlong.data.raw,Parameters.Index.averageIndex)); % If state needed for 1 NSA case

%% Spatial filter-Hamming for 3D
% outputdata.avgrawdata=Hammingfilter3D(outputdata.avgrawdata,outputdata.Parameters);disp('Hamming filter applied.')

%% Spatial FFT for each channel
dataorder=1:ndims(TRshort.avgrawdata);
% FFT for short TR
TRshort.avgrawdata=permute(TRshort.avgrawdata,dataorder);
TRshort.fftfiddata=zeros(size(TRshort.avgrawdata));
TRshort.spectradata=zeros(size(TRshort.avgrawdata));

% FFT for long TR
TRlong.avgrawdata=permute(TRlong.avgrawdata,dataorder);
TRlong.fftfiddata=zeros(size(TRlong.avgrawdata));
TRlong.spectradata=zeros(size(TRlong.avgrawdata));

%% Spectral fft
% FFT for short TR
for n=1:size(TRshort.avgrawdata,Parameters.Index.channelIndex)
    TRshort.fftfiddata(:,n,:,:,:)=Phasecorrection(ifftshift(fftshift(ifft(fftn(squeeze(TRshort.avgrawdata(:,n,:,:,:))),[],1)),1));
    TRshort.spectradata(:,n,:,:,:)=flip(fftshift(fft(TRshort.fftfiddata(:,n,:,:,:),[],1),1),1); %To check data in intermediate steps
end
% FFT for long TR
for n=1:size(TRlong.avgrawdata,Parameters.Index.channelIndex)
    TRlong.fftfiddata(:,n,:,:,:)=Phasecorrection(ifftshift(fftshift(ifft(fftn(squeeze(TRlong.avgrawdata(:,n,:,:,:))),[],1)),1));
    TRlong.spectradata(:,n,:,:,:)=flip(fftshift(fft(TRlong.fftfiddata(:,n,:,:,:),[],1),1),1); %To check data in intermediate steps
end

%% Channel combination
% Channel combination for short TR
if isempty(Parameters.Index.channelIndex)==0
    TRshort.NCov=cov(TRshort.fftfiddata((15/16)*end:end,:,1)); % Calculate noise covariance matrix from first voxel
    [TRshort.RoemerEqualfid, TRshort.RoemerEqNoiseSensitivity_maps]=RoemerEqualNoise(TRshort.fftfiddata,TRshort.NCov,Parameters.Index.channelIndex);
end
% Channel combination for long TR
if isempty(Parameters.Index.channelIndex)==0
    TRlong.NCov=cov(TRlong.fftfiddata((15/16)*end:end,:,1)); % Calculate noise covariance matrix from first voxel
    [TRlong.RoemerEqualfid, TRlong.RoemerEqNoiseSensitivity_maps]=RoemerEqualNoise_withmap_input(TRlong.fftfiddata,TRshort.RoemerEqNoiseSensitivity_maps,TRlong.NCov,Parameters.Index.channelIndex);
end

%% First order phase correction with linear prediction SVD
disp('Linear Prediction with SVD')
tic
% Linear prediction-SVD for short TR
   TRshort.RoemerEqualfid=LinearPredSVD(TRshort.RoemerEqualfid,Parameters);
% Linear prediction-SVD for long TR
   TRlong.RoemerEqualfid=LinearPredSVD(TRlong.RoemerEqualfid,Parameters);  
toc

%% Spectral apodization and zerofilling
Parameters.apodfunc=exp(-20*Parameters.time).'; %20 Hz
% Apod. and zerofill for short TR

TRshort.RoemerEqualfid=Parameters.apodfunc.*TRshort.RoemerEqualfid;
TRshort.RoemerEqualfid(end:Parameters.zerofill,:)=0;

% Apod. and zerofill for long TR

TRlong.RoemerEqualfid=Parameters.apodfunc.*TRlong.RoemerEqualfid;
TRlong.RoemerEqualfid(end:Parameters.zerofill,:)=0;

%% Final spectra
% Spectral FFT for short TR
TRshort.FinalspectraRoemer=flip(fftshift(fft(Phasecorrection(TRshort.RoemerEqualfid),[],1),1),1);
% Spectral FFT for long TR
TRlong.FinalspectraRoemer=flip(fftshift(fft(Phasecorrection(TRlong.RoemerEqualfid),[],1),1),1);

%% Watermap by peak height
TRshort.watermap=squeeze(max(real(TRshort.FinalspectraRoemer(Parameters.zerofill/2-50:Parameters.zerofill/2+50,:,:,:))));

TRlong.watermap=squeeze(max(real(TRlong.FinalspectraRoemer(Parameters.zerofill/2-50:Parameters.zerofill/2+50,:,:,:))));


%% Flip angle calculation
%TR2 is short TR scan.
%TR1 is long TR scan.
r=TRshort.watermap./TRlong.watermap;
e1=exp(-TRl/Parameters.T1);
e2=exp(-TRs/Parameters.T1);
A=((1-e2-r+r*e1)./(e1-r*e2+r*e1*e2-e1*e2));
NutationAngle=acosd(A);
FA_map = reshape(NutationAngle, [size(TRshort.watermap)]);

end