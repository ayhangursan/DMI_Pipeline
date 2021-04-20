function [FID3D]=LinearPredSVD(fiddata,Parameters)
% Linear Prediction of missing point
% Implementation from DMIWizard package
% Original code belongs to Robin de Graaf(Yale University)
dimension=size(fiddata);
FID3D=zeros(size(fiddata));
CSIgrid=dimension(2:end);
numberofvoxels=prod(CSIgrid);
missingpoints=Parameters.missingpoints+1; % Missing points +1 for the first point of the fid
disp(strcat('Predicted points:',num2str(missingpoints)));
for m=1:numberofvoxels
    if isnan(fiddata(2:end,m))==0
        %% SVD part
        nHSVD = 32; %HSVD size-number of components
        data=horzcat(fiddata(2:end,m).',0).'; % First point of FID is omitted
        th = 0;
        count = dimension(1); %Number of spectral points
        while (th < 1)
            if (abs(data(count)) == abs(data(count-1)))
                count = count - 1;
            else
                th = 1;
            end
        end
        
        npzero = dimension(1) - count;
        datasvd=zeros(size(data));
        datasvd = data(1:dimension(1)-npzero);
        
        npnonzero = length(datasvd);
        
        tacq = npnonzero/(Parameters.BW);               % Acquisition time (in s)
        dt = tacq/npnonzero;                       % Dwell-time (in s)
        time = 0:dt:(npnonzero-1)*dt;              % Time base of FID
        time = reshape(time,npnonzero,1);
        
        Lmax = round(0.4*npnonzero);               % Dimension 1 for LxM SVD matrix
        Mmax = npnonzero+1-Lmax;                   % Dimension 2 for LxM SVD matrix
        
        %*****************
        % Allocate memory
        %*****************
        H = zeros(Lmax,Mmax);
        
        %**********************************************
        % Create Hankel matrix from original FID data
        %**********************************************
        for L = 1:1:Lmax;
            M = 1:1:Mmax;
            H(L,M) = datasvd(L+M-1);
        end;
        
        %**********************************************
        % Perform SVD on Hankel matrix
        %**********************************************
        [U,~,~] = svd(H);
        
        Uup = zeros(Lmax-1,nHSVD);
        Udown = zeros(Lmax-1,nHSVD);
        
        %**********************************************
        % Calculate truncated SVD matrix
        %**********************************************
        Utr = U(:,1:1:nHSVD);
        
        for kk1 = 2:Lmax;
            for kk2 = 1:nHSVD;
                Uup(kk1-1,kk2) = Utr(kk1,kk2);
                Udown(kk1-1,kk2) = Utr(kk1-1,kk2);
            end;
        end;
        
        q = log(eig(pinv(Udown)*Uup));
        % Frequency (in Hz)
        frq = imag(q)/(2*pi*dt);
        % Time constant (in s)
        decay = real(q)/dt;
        
        time = 0:dt:(npnonzero-1)*dt;              % Time base of FID
        time = reshape(time,npnonzero,1);
        
        % Calculate basis functions
        basis = zeros(npnonzero,nHSVD);
        for kk1 = 1:1:nHSVD;
            basis(:,kk1) = exp((decay(kk1)+2*pi*1i*frq(kk1))*time);
        end;
        
        % Amplitude estimates
        ampcomplex = pinv(basis)*datasvd;
        clear basis;
        amp = abs(ampcomplex);
        phs = atan2(imag(ampcomplex),real(ampcomplex));
        
        %***********************************************************************
        % Extract signals in the frequency range [freqmin, freqmax]
        % Eliminate signals with T2 < 0.001 s
        %***********************************************************************
        coor = find((frq > 800*-1) & (frq < 800*1) & (abs(decay) < 1/(pi*0.001)));
        
        ampred = amp(coor);
        frqred = frq(coor);
        decayred = decay(coor);
        phsred = phs(coor);
        
        %% Linear Prediction of missing points
        nHSVD = length(ampred);
        
        dt = 1/(Parameters.BW);
        time = -missingpoints*dt:dt:(dimension(1)-1)*dt;
        
        datafit = zeros(1,dimension(1)+missingpoints);
        for c1 = 1:1:nHSVD;
            datafit = datafit + ampred(c1).*exp(2*pi*1i*frqred(c1)*time).*exp(time*decayred(c1)).*exp(1i*phsred(c1));
        end;
        
        % Combine the extrapolated FID points with the original FID
        datalp(1:missingpoints) = datafit(1:missingpoints);
        datalp(missingpoints+1:dimension(1)) = data(1:dimension(1)-(missingpoints));
        FID3D(:,m)=datalp;
        clearvars -except FID3D m numberofvoxels fiddata Parameters nHSVD missingpoints dimension;
    end
end
end