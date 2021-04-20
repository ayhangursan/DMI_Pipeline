function [Fittedtimesignal, Fittedspectra, quantifiedsignal]=HSVD_quantification(FID,waterloc,Parameters,isbaseline)

%% time domain signal is used as input
if nargin < 4
    isbaseline=0;
    disp('HSVD based quantification of metabolites: water, Glucose, Glx and Lipid.')
    
elseif isbaseline==1
    disp('Input is a baseline scan.')
    disp('HSVD based quantification of metabolites: water and Lipid.')
end

    
numberofvoxels=prod(Parameters.CSIdims);
NumPoints=size(FID,1);
tacq = NumPoints/(Parameters.BW);               % Acquisition time (in s)
dt = tacq/NumPoints;                       % Dwell-time (in s)
time = 0:dt:(NumPoints-1)*dt;              % Time base of FID
time = reshape(time,NumPoints,1);

Lmax = round(0.4*NumPoints);               % Dimension 1 for LxM SVD matrix
Mmax = NumPoints+1-Lmax;                   % Dimension 2 for LxM SVD matrix

H = zeros(Lmax,Mmax);

Fittedtimesignal.water=zeros([NumPoints Parameters.CSIdims]);
Fittedtimesignal.Glucose=zeros([NumPoints Parameters.CSIdims]);
Fittedtimesignal.Glx=zeros([NumPoints Parameters.CSIdims]);
Fittedtimesignal.Lipid=zeros([NumPoints Parameters.CSIdims]);
Fittedtimesignal.datafit = zeros([NumPoints Parameters.CSIdims]);

if waterloc==0
    waterloc=4.7*ones(Parameters.CSIdims);
end

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
    
    datasvd=Phasecorrection(FID(:,m));
    
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
        basis(:,kk1) = exp((decay(kk1)+2*pi*1i*frq(kk1))*time);
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
    coor = find((frq > 300*-1) & (frq < 300*1) & (abs(decay) < 1/(pi*0.001)));
    
    ampred = amp(coor);
    frqred = frq(coor);
    decayred = decay(coor);
    phsred = phs(coor);
    % Total fit
    nHSVD = length(ampred);
    Fittedsignal=zeros(size(datasvd));
    for c1 = 1:1:nHSVD
        Fittedsignal = Fittedsignal + ampred(c1).*exp(2*pi*1i*frqred(c1)*time).*exp(time*decayred(c1)).*exp(1i*phsred(c1));
    end
    Fittedtimesignal.datafit(:,m)=Fittedsignal;
    
    % Find fit for each metabolite
    % Water
    waterfreq=((waterloc(m)-4.7)+[-0.3 0.3])*45.7;
    watercoor = find((frqred > waterfreq(1)) & (frqred < waterfreq(2)));
    waterfit=zeros(size(datasvd));
    if isempty(watercoor)==0
        
        for mm = 1:numel(watercoor)
            c1=watercoor(mm);
            waterfit(:) = waterfit(:)+ampred(c1).*exp(2*pi*1i*frqred(c1)*time).*exp(time*decayred(c1)).*exp(1i*phsred(c1));
        end
        
        Fittedtimesignal.water(:,m)=waterfit;
    end
    
    % Glucose
    Glucosefreq=((waterloc(m)-4.7)+[-0.3 0.3]-0.9)*45.7;
    Glucosecoor = find((frqred > Glucosefreq(1)) & (frqred < Glucosefreq(2)));
    Glucosefit=zeros(size(datasvd));
    if isempty(Glucosecoor)==0
        
        for mm = 1:numel(Glucosecoor)
            c1=Glucosecoor(mm);
            Glucosefit(:) = Glucosefit(:)+ampred(c1).*exp(2*pi*1i*frqred(c1)*time).*exp(time*decayred(c1)).*exp(1i*phsred(c1));
        end
        Fittedtimesignal.Glucose(:,m)=Glucosefit(:);
    end
    
    % Glx
    Glxfreq=((waterloc(m)-4.7)+[-0.3 0.3]-2.3)*45.7;
    Glxcoor = find((frqred > Glxfreq(1)) & (frqred < Glxfreq(2)));
    Glxfit=zeros(size(datasvd));
    if isempty(Glxcoor)==0
        
        for mm = 1:numel(Glxcoor)
            c1=Glxcoor(mm);
            Glxfit(:) = Glxfit(:)+ampred(c1).*exp(2*pi*1i*frqred(c1)*time).*exp(time*decayred(c1)).*exp(1i*phsred(c1));
        end
        Fittedtimesignal.Glx(:,m)=Glxfit(:);
    end
    
    
    % Lipid
    Lipidfreq=((waterloc(m)-4.7)+[-0.3 0.3]-3.3)*45.7;
    Lipidcoor = find((frqred > Lipidfreq(1)) & (frqred < Lipidfreq(2)));
    Lipidfit=zeros(size(datasvd));
    if isempty(Lipidcoor)==0
        
        for mm = 1:numel(Lipidcoor)
            c1=Lipidcoor(mm);
            Lipidfit(:) = Lipidfit(:)+ampred(c1).*exp(2*pi*1i*frqred(c1)*time).*exp(time*decayred(c1)).*exp(1i*phsred(c1));
        end
        Fittedtimesignal.Lipid(:,m)=Lipidfit(:);
    end
    
    
    
end

Fittedspectra.water=real(fftshift(fft(Phasecorrection(Fittedtimesignal.water),[],1),1));
Fittedspectra.Glucose=real(fftshift(fft(Phasecorrection(Fittedtimesignal.Glucose),[],1),1));
Fittedspectra.Glx=real(fftshift(fft(Phasecorrection(Fittedtimesignal.Glx),[],1),1));
Fittedspectra.Lipid=real(fftshift(fft(Phasecorrection(Fittedtimesignal.Lipid),[],1),1));
Fittedspectra.datafit = real(fftshift(fft(Phasecorrection(Fittedtimesignal.datafit),[],1),1));


quantifiedsignal.water=reshape(sum(real(Fittedspectra.water),1),[Parameters.CSIdims]);
quantifiedsignal.Glucose=reshape(sum(real(Fittedspectra.Glucose),1),[Parameters.CSIdims]);
quantifiedsignal.Glx=reshape(sum(real(Fittedspectra.Glx),1),[Parameters.CSIdims]);
quantifiedsignal.Lipid=reshape(sum(real(Fittedspectra.Lipid),1),[Parameters.CSIdims]);

end