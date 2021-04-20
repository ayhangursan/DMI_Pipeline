
% Script to calculate coefficient of variation(cov) in various datasets.
% CoV was calculated on abundance water signal.

% V9164 - high spatial resolution with isotropic voxels
% V9158 - high in-plane spatial resolution with anisotropic voxels
% V9077 - Low spatial resolution
% V9073 - Coefficient of variation calculation on natural
% V9066 - low spatial resolution + EPSI Scans
% Date:09-04-2021

%% V9073
cd('C:\Users\agursan.DS\surfdrive\Deuterium_CSI\Exams\2020-06-15 2H in-vivo V9073\listdata');
V9073.Noise=GlobalFID_NCovGeneration('raw_020.data');
Gridshift=[0 0 0];
V9073.Scan1=DMIPipeline3D_refscan('raw_032.data',0,V9073.Noise.noisecovariance,1.95,Gridshift,0);
V9073.Scan2=DMIPipeline3D_refscan('raw_033.data',V9073.Scan1.RoemerEqNoiseSensitivity_maps,V9073.Noise.noisecovariance,1.95,Gridshift,0);
V9073.Scan3=DMIPipeline3D_refscan('raw_034.data',V9073.Scan1.RoemerEqNoiseSensitivity_maps,V9073.Noise.noisecovariance,1.95,Gridshift,0);
V9073.Scan4=DMIPipeline3D_refscan('raw_035.data',V9073.Scan1.RoemerEqNoiseSensitivity_maps,V9073.Noise.noisecovariance,1.95,Gridshift,0);
V9073.Scan5=DMIPipeline3D_refscan('raw_036.data',V9073.Scan1.RoemerEqNoiseSensitivity_maps,V9073.Noise.noisecovariance,1.95,Gridshift,0);
V9073.Scan6=DMIPipeline3D_refscan('raw_037.data',V9073.Scan1.RoemerEqNoiseSensitivity_maps,V9073.Noise.noisecovariance,1.95,Gridshift,0);

% CoV calculations for 5 scans
% define water window
waterwindowppm=[4 5.4];
waterwindowindex=find(V9073.Scan1.xaxiszerofill> waterwindowppm(1) & V9073.Scan1.xaxiszerofill< waterwindowppm(2));
%% CoV by peak heights
V9073.PeakHeights=zeros([5 V9073.Scan1.Parameters.CSIdims]);
V9073.PeakHeights(1,:,:,:)=squeeze(max(real(V9073.Scan1.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9073.PeakHeights(2,:,:,:)=squeeze(max(real(V9073.Scan2.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9073.PeakHeights(3,:,:,:)=squeeze(max(real(V9073.Scan3.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9073.PeakHeights(4,:,:,:)=squeeze(max(real(V9073.Scan4.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9073.PeakHeights(5,:,:,:)=squeeze(max(real(V9073.Scan5.FinalspectraRoemer(waterwindowindex,:,:,:))));
% V9073.PeakHeights(6,:,:,:)=squeeze(max(real(V9073.Scan6.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9073.CoV=100*squeeze(std(V9073.PeakHeights,[],1)./mean(V9073.PeakHeights,1));

%% CoV by peak heights - Raw data
PhaFunct=exp(-1i * (2* pi * (V9073.Scan1.xaxiszerofill-4.7).'*(V9073.Scan1.Parameters.Freq/(10^6)) * V9073.Scan1.Parameters.TE));
for mm=1:5
    eval(['TimepointFID=padarray(V9073.Scan',num2str(mm),'.RoemerEqualfid.*V9073.Scan1.Parameters.apodfunc,[1024 0 0 0],0,''post'');'])
    TimepointSPEC=real(fftshift(fft(TimepointFID,[],1),1).*PhaFunct);
    eval(['V9073.RawPeakHeights(',num2str(mm),',:,:,:)=squeeze(max(TimepointSPEC(waterwindowindex,:,:,:)));'])
    clear TimepointSPEC TimepointFID;
end
V9073.RawCoV=100*squeeze(std(V9073.RawPeakHeights,[],1)./mean(V9073.RawPeakHeights,1));
%% CoV by AMARES fits
for mm=1:5
    eval(['TimepointAMARES=cell2mat(V9073.Scan',num2str(mm),'.AMARES_Results);'])
    b=zeros(numel(TimepointAMARES),2);
    for kk=1:numel(TimepointAMARES)
        b(kk,:)=TimepointAMARES(kk).amplitude;
    end
    wateramp=reshape(b(:,2),[size(TimepointAMARES)]);
    eval(['V9073.AMARESAmplitude(',num2str(mm),',:,:,:)=wateramp;'])
    clear TimepointAMARES wateramp;
end
V9073.AMARESCoV=100*squeeze(std(V9073.AMARESAmplitude,[],1)./mean(V9073.AMARESAmplitude,1));
%% V9164
% Voxel size 25 x 25 x 25 mm
cd('C:\Users\agursan.DS\surfdrive\Deuterium_CSI\Exams\2020-10-02 2H in-vivo V9164 25mm iso voxels\listdata');
Halfvoxelshift=[0 0 -0.5];
V9164.Noise=DMINoiseCovariance('raw_202.data');
V9164.Scan1=DMIPipeline3D_refscan('raw_030.data',0,V9164.Noise.noisecovariance,1.95,Halfvoxelshift,0);
V9164.Scan2=DMIPipeline3D_refscan('raw_031.data',V9164.Scan1.RoemerEqNoiseSensitivity_maps,V9164.Noise.noisecovariance,1.95,Halfvoxelshift,0);
V9164.Scan3=DMIPipeline3D_refscan('raw_032.data',V9164.Scan1.RoemerEqNoiseSensitivity_maps,V9164.Noise.noisecovariance,1.95,Halfvoxelshift,0);
V9164.Scan4=DMIPipeline3D_refscan('raw_033.data',V9164.Scan1.RoemerEqNoiseSensitivity_maps,V9164.Noise.noisecovariance,1.95,Halfvoxelshift,0);
V9164.Scan5=DMIPipeline3D_refscan('raw_034.data',V9164.Scan1.RoemerEqNoiseSensitivity_maps,V9164.Noise.noisecovariance,1.95,Halfvoxelshift,0);
%% CoV by peak heights
V9164.PeakHeights=zeros([5 V9164.Scan1.Parameters.CSIdims]);
V9164.PeakHeights(1,:,:,:)=squeeze(max(real(V9164.Scan1.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9164.PeakHeights(2,:,:,:)=squeeze(max(real(V9164.Scan2.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9164.PeakHeights(3,:,:,:)=squeeze(max(real(V9164.Scan3.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9164.PeakHeights(4,:,:,:)=squeeze(max(real(V9164.Scan4.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9164.PeakHeights(5,:,:,:)=squeeze(max(real(V9164.Scan5.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9164.CoV=100*squeeze(std(V9164.PeakHeights,[],1)./mean(V9164.PeakHeights,1));
%% CoV by peak heights - Raw data
PhaFunct=exp(-1i * (2* pi * (V9164.Scan1.xaxiszerofill-4.7).'*(V9164.Scan1.Parameters.Freq/(10^6)) * V9164.Scan1.Parameters.TE));
for mm=1:5
    eval(['TimepointFID=padarray(V9164.Scan',num2str(mm),'.RoemerEqualfid.*V9164.Scan1.Parameters.apodfunc,[1024 0 0 0],0,''post'');'])
    TimepointSPEC=real(fftshift(fft(TimepointFID,[],1),1).*PhaFunct);
    eval(['V9164.RawPeakHeights(',num2str(mm),',:,:,:)=squeeze(max(TimepointSPEC(waterwindowindex,:,:,:)));'])
    clear TimepointSPEC TimepointFID;
end
V9164.RawCoV=100*squeeze(std(V9164.RawPeakHeights,[],1)./mean(V9164.RawPeakHeights,1));
%% CoV by AMARES fits
for mm=1:5
    eval(['TimepointAMARES=cell2mat(V9164.Scan',num2str(mm),'.AMARES_Results);'])
    b=zeros(numel(TimepointAMARES),2);
    for kk=1:numel(TimepointAMARES)
        b(kk,:)=TimepointAMARES(kk).amplitude;
    end
    wateramp=reshape(b(:,2),[size(TimepointAMARES)]);
    eval(['V9164.AMARESAmplitude(',num2str(mm),',:,:,:)=wateramp;'])
    clear TimepointAMARES wateramp;
end
V9164.AMARESCoV=100*squeeze(std(V9164.AMARESAmplitude,[],1)./mean(V9164.AMARESAmplitude,1));
%% V9158
% Voxel size 22 x 22 x 30
cd('C:\Users\agursan.DS\surfdrive\Deuterium_CSI\Exams\2020-09-25 2H in-vivo V9158 IncreasedRes in plane\listdata')
Halfvoxelshift=[0 0 -0.5];
V9158.Noise=DMINoiseCovariance('raw_200.data');
V9158.Scan1=DMIPipeline3D_refscan('raw_030.data',0,V9158.Noise.noisecovariance,1.95,Halfvoxelshift,0);
V9158.Scan2=DMIPipeline3D_refscan('raw_031.data',V9158.Scan1.RoemerEqNoiseSensitivity_maps,V9158.Noise.noisecovariance,1.95,Halfvoxelshift,0);
V9158.Scan3=DMIPipeline3D_refscan('raw_032.data',V9158.Scan1.RoemerEqNoiseSensitivity_maps,V9158.Noise.noisecovariance,1.95,Halfvoxelshift,0);
V9158.Scan4=DMIPipeline3D_refscan('raw_033.data',V9158.Scan1.RoemerEqNoiseSensitivity_maps,V9158.Noise.noisecovariance,1.95,Halfvoxelshift,0);
V9158.Scan5=DMIPipeline3D_refscan('raw_034.data',V9158.Scan1.RoemerEqNoiseSensitivity_maps,V9158.Noise.noisecovariance,1.95,Halfvoxelshift,0);
%% V9158 CoV by peak heights
V9158.PeakHeights=zeros([5 V9158.Scan1.Parameters.CSIdims]);
V9158.PeakHeights(1,:,:,:)=squeeze(max(real(V9158.Scan1.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9158.PeakHeights(2,:,:,:)=squeeze(max(real(V9158.Scan2.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9158.PeakHeights(3,:,:,:)=squeeze(max(real(V9158.Scan3.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9158.PeakHeights(4,:,:,:)=squeeze(max(real(V9158.Scan4.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9158.PeakHeights(5,:,:,:)=squeeze(max(real(V9158.Scan5.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9158.CoV=100*squeeze(std(V9158.PeakHeights,[],1)./mean(V9158.PeakHeights,1));
%% CoV by peak heights - Raw data
PhaFunct=exp(-1i * (2* pi * (V9158.Scan1.xaxiszerofill-4.7).'*(V9158.Scan1.Parameters.Freq/(10^6)) * V9158.Scan1.Parameters.TE));
for mm=1:5
    eval(['TimepointFID=padarray(V9158.Scan',num2str(mm),'.RoemerEqualfid.*V9158.Scan1.Parameters.apodfunc,[1024 0 0 0],0,''post'');'])
    TimepointSPEC=real(fftshift(fft(TimepointFID,[],1),1).*PhaFunct);
    eval(['V9158.RawPeakHeights(',num2str(mm),',:,:,:)=squeeze(max(TimepointSPEC(waterwindowindex,:,:,:)));'])
    clear TimepointSPEC TimepointFID;
end
V9158.RawCoV=100*squeeze(std(V9158.RawPeakHeights,[],1)./mean(V9158.RawPeakHeights,1));
%% CoV by AMARES fits
for mm=1:5
    eval(['TimepointAMARES=cell2mat(V9158.Scan',num2str(mm),'.AMARES_Results);'])
    b=zeros(numel(TimepointAMARES),2);
    for kk=1:numel(TimepointAMARES)
        b(kk,:)=TimepointAMARES(kk).amplitude;
    end
    wateramp=reshape(b(:,2),[size(TimepointAMARES)]);
    eval(['V9158.AMARESAmplitude(',num2str(mm),',:,:,:)=wateramp;'])
    clear TimepointAMARES wateramp;
end
V9158.AMARESCoV=100*squeeze(std(V9158.AMARESAmplitude,[],1)./mean(V9158.AMARESAmplitude,1));
%% V9077
cd('C:\Users\agursan.DS\surfdrive\Deuterium_CSI\Exams\2020-06-22 2H in-vivo V9077 3DCSI with ref phantom\listdata');
V9077.Noise=GlobalFID_NCovGeneration('raw_020.data');
V9077.Scan1=DMIPipeline3D_refscan('raw_041.data',0,V9077.Noise.noisecovariance,1.95,Halfvoxelshift,0);
V9077.Scan2=DMIPipeline3D_refscan('raw_042.data',V9077.Scan1.RoemerEqNoiseSensitivity_maps,V9077.Noise.noisecovariance,1.95,Halfvoxelshift,0);
V9077.Scan3=DMIPipeline3D_refscan('raw_043.data',V9077.Scan1.RoemerEqNoiseSensitivity_maps,V9077.Noise.noisecovariance,1.95,Halfvoxelshift,0);
V9077.Scan4=DMIPipeline3D_refscan('raw_044.data',V9077.Scan1.RoemerEqNoiseSensitivity_maps,V9077.Noise.noisecovariance,1.95,Halfvoxelshift,0);
V9077.Scan5=DMIPipeline3D_refscan('raw_045.data',V9077.Scan1.RoemerEqNoiseSensitivity_maps,V9077.Noise.noisecovariance,1.95,Halfvoxelshift,0);
%% V9077 CoV by peak heights

V9077.PeakHeights=zeros([5 V9077.Scan1.Parameters.CSIdims]);
V9077.PeakHeights(1,:,:,:)=squeeze(max(real(V9077.Scan1.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9077.PeakHeights(2,:,:,:)=squeeze(max(real(V9077.Scan2.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9077.PeakHeights(3,:,:,:)=squeeze(max(real(V9077.Scan3.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9077.PeakHeights(4,:,:,:)=squeeze(max(real(V9077.Scan4.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9077.PeakHeights(5,:,:,:)=squeeze(max(real(V9077.Scan5.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9077.CoV=100*squeeze(std(V9077.PeakHeights,[],1)./mean(V9077.PeakHeights,1));
%% CoV by peak heights - Raw data
PhaFunct=exp(-1i * (2* pi * (V9077.Scan1.xaxiszerofill-4.7).'*(V9077.Scan1.Parameters.Freq/(10^6)) * V9077.Scan1.Parameters.TE));
for mm=1:5
    eval(['TimepointFID=padarray(V9077.Scan',num2str(mm),'.RoemerEqualfid.*V9077.Scan1.Parameters.apodfunc,[1024 0 0 0],0,''post'');'])
    TimepointSPEC=real(fftshift(fft(TimepointFID,[],1),1).*PhaFunct);
    eval(['V9077.RawPeakHeights(',num2str(mm),',:,:,:)=squeeze(max(TimepointSPEC(waterwindowindex,:,:,:)));'])
    clear TimepointSPEC TimepointFID;
end
V9077.RawCoV=100*squeeze(std(V9077.RawPeakHeights,[],1)./mean(V9077.RawPeakHeights,1));
%% CoV by AMARES fits
for mm=1:5
    eval(['TimepointAMARES=cell2mat(V9077.Scan',num2str(mm),'.AMARES_Results);'])
    b=zeros(numel(TimepointAMARES),2);
    for kk=1:numel(TimepointAMARES)
        b(kk,:)=TimepointAMARES(kk).amplitude;
    end
    wateramp=reshape(b(:,2),[size(TimepointAMARES)]);
    eval(['V9077.AMARESAmplitude(',num2str(mm),',:,:,:)=wateramp;'])
    clear TimepointAMARES wateramp;
end
V9077.AMARESCoV=100*squeeze(std(V9077.AMARESAmplitude,[],1)./mean(V9077.AMARESAmplitude,1));
%% V9066
cd('C:\Users\agursan.DS\surfdrive\Deuterium_CSI\Exams\2020-06-09 2H in-vivo V9066\listdata')
Halfvoxelshift=[0 0 -0.5];
V9066.Noise=GlobalFID_NCovGeneration('raw_020.data');
V9066.Scan1=DMIPipeline3D_refscan('raw_031.data',0,V9066.Noise.noisecovariance,1.95,Halfvoxelshift,0);
V9066.Scan2=DMIPipeline3D_refscan('raw_032.data',V9066.Scan1.RoemerEqNoiseSensitivity_maps,V9066.Noise.noisecovariance,1.95,Halfvoxelshift,0);
V9066.Scan3=DMIPipeline3D_refscan('raw_033.data',V9066.Scan1.RoemerEqNoiseSensitivity_maps,V9066.Noise.noisecovariance,1.95,Halfvoxelshift,0);
V9066.Scan4=DMIPipeline3D_refscan('raw_034.data',V9066.Scan1.RoemerEqNoiseSensitivity_maps,V9066.Noise.noisecovariance,1.95,Halfvoxelshift,0);
V9066.Scan5=DMIPipeline3D_refscan('raw_035.data',V9066.Scan1.RoemerEqNoiseSensitivity_maps,V9066.Noise.noisecovariance,1.95,Halfvoxelshift,0);
%% V9066 CoV by peak heights
V9066.PeakHeights=zeros([5 V9066.Scan1.Parameters.CSIdims]);
V9066.PeakHeights(1,:,:,:)=squeeze(max(real(V9066.Scan1.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9066.PeakHeights(2,:,:,:)=squeeze(max(real(V9066.Scan2.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9066.PeakHeights(3,:,:,:)=squeeze(max(real(V9066.Scan3.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9066.PeakHeights(4,:,:,:)=squeeze(max(real(V9066.Scan4.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9066.PeakHeights(5,:,:,:)=squeeze(max(real(V9066.Scan5.FinalspectraRoemer(waterwindowindex,:,:,:))));
V9066.CoV=100*squeeze(std(V9066.PeakHeights,[],1)./mean(V9066.PeakHeights,1));
%% CoV by peak heights - Raw data
PhaFunct=exp(-1i * (2* pi * (V9066.Scan1.xaxiszerofill-4.7).'*(V9066.Scan1.Parameters.Freq/(10^6)) * V9066.Scan1.Parameters.TE));
for mm=1:5
    eval(['TimepointFID=padarray(V9066.Scan',num2str(mm),'.RoemerEqualfid.*V9066.Scan1.Parameters.apodfunc,[1024 0 0 0],0,''post'');'])
    TimepointSPEC=real(fftshift(fft(TimepointFID,[],1),1).*PhaFunct);
    eval(['V9066.RawPeakHeights(',num2str(mm),',:,:,:)=squeeze(max(TimepointSPEC(waterwindowindex,:,:,:)));'])
    clear TimepointSPEC TimepointFID;
end
V9066.RawCoV=100*squeeze(std(V9066.RawPeakHeights,[],1)./mean(V9066.RawPeakHeights,1));
%% CoV by AMARES fits
for mm=1:5
    eval(['TimepointAMARES=cell2mat(V9066.Scan',num2str(mm),'.AMARES_Results);'])
    b=zeros(numel(TimepointAMARES),2);
    for kk=1:numel(TimepointAMARES)
        b(kk,:)=TimepointAMARES(kk).amplitude;
    end
    wateramp=reshape(b(:,2),[size(TimepointAMARES)]);
    eval(['V9066.AMARESAmplitude(',num2str(mm),',:,:,:)=wateramp;'])
    clear TimepointAMARES wateramp;
end
V9066.AMARESCoV=100*squeeze(std(V9066.AMARESAmplitude,[],1)./mean(V9066.AMARESAmplitude,1));


%% V9046 - slice mismatch EXCLUDE
% cd('C:\Users\agursan.DS\surfdrive\Deuterium_CSI\Exams\2020-05-14 2H in-vivo V9046 low res series\listdata')
% V9046.Noise=GlobalFID_NCovGeneration('raw_020.data');
% V9046.Scan1=DMIPipeline3D_refscan('raw_031.data',0,V9046.Noise.noisecovariance,1.95,Halfvoxelshift,0);
% V9046.Scan2=DMIPipeline3D_refscan('raw_032.data',V9046.Scan1.RoemerEqNoiseSensitivity_maps,V9046.Noise.noisecovariance,1.95,Halfvoxelshift,0);
% V9046.Scan3=DMIPipeline3D_refscan('raw_033.data',V9046.Scan1.RoemerEqNoiseSensitivity_maps,V9046.Noise.noisecovariance,1.95,Halfvoxelshift,0);
% V9046.Scan4=DMIPipeline3D_refscan('raw_034.data',V9046.Scan1.RoemerEqNoiseSensitivity_maps,V9046.Noise.noisecovariance,1.95,Halfvoxelshift,0);
% V9046.Scan5=DMIPipeline3D_refscan('raw_035.data',V9046.Scan1.RoemerEqNoiseSensitivity_maps,V9046.Noise.noisecovariance,1.95,Halfvoxelshift,0);
% %% V9046 CoV by peak heights
% V9046.PeakHeights=zeros([5 V9046.Scan1.Parameters.CSIdims]);
% V9046.PeakHeights(1,:,:,:)=squeeze(max(real(V9046.Scan1.FinalspectraRoemer(waterwindowindex,:,:,:))));
% V9046.PeakHeights(2,:,:,:)=squeeze(max(real(V9046.Scan2.FinalspectraRoemer(waterwindowindex,:,:,:))));
% V9046.PeakHeights(3,:,:,:)=squeeze(max(real(V9046.Scan3.FinalspectraRoemer(waterwindowindex,:,:,:))));
% V9046.PeakHeights(4,:,:,:)=squeeze(max(real(V9046.Scan4.FinalspectraRoemer(waterwindowindex,:,:,:))));
% V9046.PeakHeights(5,:,:,:)=squeeze(max(real(V9046.Scan5.FinalspectraRoemer(waterwindowindex,:,:,:))));
% V9046.CoV=100*squeeze(std(V9046.PeakHeights,[],1)./mean(V9046.PeakHeights,1));
%% Script for ROI drawing on image and overlay
% For V9158 ROI is saved
% For V9164 ROI is saved

% cd('C:\Users\agursan.DS\surfdrive\Deuterium_CSI\Exams\2020-06-15 2H in-vivo V9073\V9073\DICOM');dicominput='IM_0012';FOV=[240 300]; % Dicom V9073
% cd('C:\Users\agursan.DS\surfdrive\Deuterium_CSI\Exams\2020-10-02 2H in-vivo V9164 25mm iso voxels\DICOM');dicominput='IM_0018';FOV=[250 300]; % Dicom V9164
% cd('C:\Users\agursan.DS\surfdrive\Deuterium_CSI\Exams\2020-09-25 2H in-vivo V9158 IncreasedRes in plane\Multix\DICOM');dicominput='IM_0010';FOV=[242 308]; % Dicom V9158
% cd('C:\Users\agursan.DS\surfdrive\Deuterium_CSI\Exams\2020-06-22 2H in-vivo V9077 3DCSI with ref phantom\dicom\DICOM');dicominput='IM_0012';FOV=[240 300]; % Dicom V9077
cd('C:\Users\agursan.DS\surfdrive\Deuterium_CSI\Exams\2020-06-09 2H in-vivo V9066\dicom\DICOM');dicominput='IM_0012';FOV=[240 300]; % Dicom V9077
% cd('C:\Users\agursan.DS\surfdrive\Deuterium_CSI\Exams\2020-05-14 2H in-vivo V9046 low res series\dicom\DICOM');dicominput='IM_0016';FOV=[240 300]; % Dicom V9077


Dataset=V9066;
Parameters=Dataset.Scan1.Parameters;
inputCSI=Parameters.CSIdims;
% CoV=Dataset.RawCoV;


% End Debug
Images.DicomImage=squeeze(double(dicomread(dicominput)));
Images.DicomTagInfo=dicominfo(dicominput);
NoSlices=size(Images.DicomImage,3);
Midslice=ceil(NoSlices/2);
OrderedSlices=[Midslice+1:NoSlices 1:Midslice];
Images.DicomImage=Images.DicomImage(:,:,OrderedSlices);% Reorder linearly
% Crop with respect to FoV
imageratio=FOV(1)/FOV(2);
ImageSize=[size(Images.DicomImage,1) size(Images.DicomImage,2)];
Filledpart=ImageSize(1)*(1-imageratio);
Images.DicomImage=squeeze(Images.DicomImage(ceil(Filledpart/2):(ImageSize(1)-floor(Filledpart/2)),:,:));
Images.DrawnROImask=NaN(inputCSI);
%%
Axialslice=6;
figure('WindowState','maximized')
imagesc(Images.DicomImage(:,:,Axialslice))
daspect([1 1 1])
set(gca, 'Layer','top')
axis on;
[rows, columns, numberOfColorChannels] = size(Images.DicomImage(:,:,1));
hold on;
for row = 1 : rows ./ inputCSI(1) : rows
    line([1, columns], [row, row], 'Color', 'r','Linewidth',1.5);
end
for col = 1 : columns ./ inputCSI(2) : columns
    line([col, col], [1, rows], 'Color', 'r','Linewidth',1.5);
end
colormap gray
Images.DrawnROI=drawfreehand;
Images.DrawnROImask(:,:,Axialslice)=imresize(createMask(Images.DrawnROI),[inputCSI(1) inputCSI(2)]);
%%
% CoV=Dataset.RawCoV;
% CoV=Dataset.CoV;
CoV=Dataset.AMARESCoV;
figure('WindowState','maximized')
clear TotalCoV
TotalCoV(1)=NaN;
tilenum=1;
plotnum=1;
for mm=4:6
        ax(plotnum)=subplot(4,16,[tilenum:tilenum+3 tilenum+16:tilenum+19]);
%     ax(plotnum)=subplot(2,16,[tilenum:tilenum+3]);
    
    imagesc(Images.DicomImage(:,:,mm).^0.8)
    daspect([1 1 1])
    set(gca, 'Layer','top')
    axis on;
    hold on;
    for row = 1 : rows ./ inputCSI(1) : rows
        line([1, columns], [row, row], 'Color', 'r');
    end
    for col = 1 : columns ./ inputCSI(2) : columns
        line([col, col], [1, rows], 'Color', 'r');
    end
    colormap(ax(plotnum),gray)
    title(['Axial Image slice: ',num2str(mm)],'FontSize',20)
    set(gca,'YTickLabel',[],'XTickLabel',[])
    
        ax(plotnum+4)=subplot(4,16,[tilenum+32:tilenum+35 tilenum+48:tilenum+51]);
%     ax(plotnum+4)=subplot(2,16,[tilenum+16:tilenum+16+3]);
    
    imagesc(CoV(:,:,mm).*Images.DrawnROImask(:,:,mm))
    daspect([1 1 1])
    set(gca, 'Layer','top')
    axis on;
    hold on;
    for row = 1:inputCSI(1)
        line([0, inputCSI(2)+1], [row-0.5, row-0.5], 'Color', 'r');
    end
    for col = 1 : inputCSI(2)
        line([col+0.5, col+0.5], [0, row+0.5], 'Color', 'r');
    end
    colormap(ax(plotnum+4),parula)
    %     colorbar('FontSize',16)
    caxis([0 15])
    sliceCoV=CoV(:,:,mm);
    sliceCoV=sliceCoV(logical(Images.DrawnROImask(:,:,mm)));
    slicemeanCoV=mean(sliceCoV,'all','omitnan');
    slicestdCoV=std(sliceCoV,[],'omitnan');
    titlecalculation=['Mean: ',num2str(slicemeanCoV,3),' (\pm',num2str(slicestdCoV,3),')'];
    title(["Masked CoV",titlecalculation],'FontSize',20)
    %     title(["You can do it","with a string array too"])
    
    TotalCoV=[TotalCoV sliceCoV.'];
    tilenum=tilenum+4;
    plotnum=plotnum+1;
end
TotalmeanCoV=mean(TotalCoV,'all','omitnan');
TotalstdCoV=std(TotalCoV,[],'omitnan');
sgtitle(['V9164-CoV in all mask mean(std): ',num2str(TotalmeanCoV,3),'(',num2str(TotalstdCoV,3),')'],'FontSize',24)

size(TotalCoV)