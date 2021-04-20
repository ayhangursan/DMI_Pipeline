function [Splittings]=HSVD_QuadrupolarSplittigCalculation(FID,Parameters)

%% time domain signal is used as input

%DEBUG
% clear FID;
% FID=squeeze(V9259.LegCSI.RoemerEqualfid(:,3:4,3:4,6)); % DEBUG FID
%Blue boxes
% Scan 1
% FID(:,1,1)=squeeze(V9259.LegCSI.RoemerEqualfid(:,4,7,6));
% FID(:,1,2)=squeeze(V9259.LegCSI.RoemerEqualfid(:,5,7,6));
% FID(:,2,1)=squeeze(V9259.LegCSI.RoemerEqualfid(:,6,6,6));
% FID(:,2,2)=squeeze(V9259.LegCSI.RoemerEqualfid(:,7,5,6));
%Red boxes
% Scan 1
% FID(:,1,1)=squeeze(V9259.LegCSI.RoemerEqualfid(:,5,9,6));
% FID(:,1,2)=squeeze(V9259.LegCSI.RoemerEqualfid(:,6,9,6));
% FID(:,2,1)=squeeze(V9259.LegCSI.RoemerEqualfid(:,7,8,6));
% FID(:,2,2)=squeeze(V9259.LegCSI.RoemerEqualfid(:,8,8,6));

numberofvoxels=numel(FID)./size(FID,1);
NumPoints=size(FID,1);
tacq = NumPoints/(Parameters.BW);               % Acquisition time (in s)
dt = tacq/NumPoints;                       % Dwell-time (in s)
acqtime = 0:dt:(NumPoints-1)*dt;              % Time base of FID
acqtime = reshape(acqtime,NumPoints,1);
xaxis=linspace(-Parameters.ppmwindow/2,Parameters.ppmwindow/2,Parameters.NP)+4.7;
missingpoints=Parameters.missingpoints;
Lmax = round(0.4*NumPoints);               % Dimension 1 for LxM SVD matrix
Mmax = NumPoints+1-Lmax;                   % Dimension 2 for LxM SVD matrix

H = zeros(Lmax,Mmax);
Splittings=zeros([size(FID,2) size(FID,3) size(FID,4)]);

for m=1:numberofvoxels
    if m==floor(numberofvoxels/4)
        disp('25% of the voxels are quantified')
    elseif m==floor(numberofvoxels/2)
        disp('50% of the voxels are quantified')
    elseif m==floor(3*numberofvoxels/4)
        disp('75% of the voxels are quantified')
    elseif m==numberofvoxels
        disp('100% of the voxels are quantified')
    end
    %
%     datasvd=Phasecorrection(FID(:,m));
datasvd=FID(:,m); % No 0th order phase correction
    for L = 1:1:Lmax
        M = 1:1:Mmax;
        H(L,M) = datasvd(L+M-1);
    end
    
    %**********************************************
    % Perform SVD on Hankel matrix
    %**********************************************
    [U,~,~] = svd(H);
    nHSVD=32;
    
    Uup = zeros(Lmax-1,nHSVD);
    Udown = zeros(Lmax-1,nHSVD);
    %**********************************************
    % Calculate truncated SVD matrix
    %**********************************************
    Utr = U(:,1:1:nHSVD); % Matrix truncated to 32 elements
    
    for kk1 = 2:Lmax
        for kk2 = 1:nHSVD
            Uup(kk1-1,kk2) = Utr(kk1,kk2);
            Udown(kk1-1,kk2) = Utr(kk1-1,kk2);
        end
    end
    
    q = log(eig(pinv(Udown)*Uup));
    
    % Frequency (in Hz)
    frq = imag(q)/(2*pi*dt);
    % Time constant (in s)
    decay = real(q)/dt;
    % Calculate basis functions
    basis = zeros(NumPoints,nHSVD);
    for kk1 = 1:1:nHSVD
        basis(:,kk1) = exp((decay(kk1)+2*pi*1i*frq(kk1))*acqtime);
    end
    
    % Amplitude estimates
    ampcomplex = pinv(basis)*datasvd;
    clear basis;
    amp = abs(ampcomplex);
    phs = atan2(imag(ampcomplex),real(ampcomplex));
    
    %***********************************************************************
    % Extract signals in the frequency range [250, 250]
    % Eliminate signals with T2 < 0.001 s
    %***********************************************************************
    coor = find((frq > 250*-1) & (frq < 250*1) & (abs(decay) < 1/(pi*0.001)));
    
    ampred = amp(coor);
    frqred = frq(coor);
    decayred = decay(coor);
    phsred = phs(coor);
    
    [~ , ampredorder]=sort(ampred);
    disp('[Freq(Hz)  Freq(ppm)   Amplitude   Phase 1st-Phase(rad)]')
    disp([(frqred(ampredorder)) ((frqred(ampredorder))/45.7)+4.7 ampred(ampredorder) phsred(ampredorder) 2*pi*Parameters.TE*(frqred(ampredorder))])
        disp([(frqred(ampredorder)) ((frqred(ampredorder))/45.7)+4.7 ampred(ampredorder) phsred(ampredorder)-2*pi*Parameters.TE*(frqred(ampredorder))])

    % Total fit
    nHSVD = length(ampred);
    Fittedsignal=zeros(size(datasvd));
    for c1 = 1:1:nHSVD
        Fittedsignal = Fittedsignal + ampred(c1).*exp(2*pi*1i*frqred(c1)*acqtime).*exp(acqtime*decayred(c1)).*exp(-1i*phsred(c1));
    end
    
    %     % Plot raw data and fit
    figure('WindowState','maximized')
    sgtitle(strcat('Voxel number:',num2str(m)))
    subplot(2,2,1)
    plot(xaxis,abs(fftshift(fft(datasvd))),'LineWidth',2)
    hold on
    plot(xaxis,real(fftshift(fft(datasvd))),'LineWidth',2)
    plot(xaxis,real(fftshift(fft(Fittedsignal))),'LineWidth',2)
    hold off
    xlim([0 10])
    set(gca,'XDir','reverse')
    legend('Raw-Mag','Raw-Real','HSVD','FontSize',20)
    title('Default Phase')
    
    subplot(2,2,2)
    plot(xaxis,abs(fftshift(fft(datasvd))),'LineWidth',2)
    hold on
    plot(xaxis,real(fftshift(fft(datasvd))),'LineWidth',2)
    
    phase=atan2(imag(Fittedsignal(2)),real(Fittedsignal(2)));
    plot(xaxis,real(fftshift(fft(Fittedsignal.*(exp(-1i*phase))))),'LineWidth',2)
    hold off
    xlim([0 10])
    set(gca,'XDir','reverse')
    legend('Raw-Mag','Raw-Real','HSVD','FontSize',20)
    title('Extra Phase correction applied')
    
    
    if numel(frqred)<2
        Splittings(m)=NaN;
    else
        Splittings(m)=abs(frqred(ampredorder(end))-frqred(ampredorder(end-1)));
    end
    
    nHSVD = length(ampred);
        
        dt = 1/(Parameters.BW);
        time = -missingpoints*dt:dt:(NumPoints-1)*dt;
        
        datafit = zeros(1,NumPoints+missingpoints);
        for c1 = 1:1:nHSVD
            datafit = datafit + ampred(c1).*exp(2*pi*1i*frqred(c1)*time).*exp(time*decayred(c1)).*exp(1i*phsred(c1));
        end
        
        % Combine the extrapolated FID points with the original FID
        datalp(1:missingpoints) = datafit(1:missingpoints);
        datalp(missingpoints+1:NumPoints) = datasvd(1:NumPoints-(missingpoints));
        phase=atan2(imag(datalp(2)),real(datalp(2)));

        subplot(2,2,3)
        plot(xaxis,abs(fftshift(fft(datasvd))),'LineWidth',2)
        hold on
        plot(xaxis,real(fftshift(fft(datasvd))),'LineWidth',2)
        plot(xaxis,real(fftshift(fft(datalp))),'LineWidth',2)
        plot(xaxis,real(fftshift(fft(datalp.*(exp(-1i*phase))))),'LineWidth',2)
        hold off
        xlim([0 10])
        set(gca,'XDir','reverse')
        legend('Raw-Mag','Raw-Real','HSVD-predicted','HSVD-predict+0Ph corr','FontSize',18,'Location','eastoutside')
        title('Extra Phase correction applied')
        
        
end

end