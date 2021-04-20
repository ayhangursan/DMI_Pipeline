function [WSVDcombined, wsvdQuality, wsvdWeights]=WSVDcomb(fiddata,NCov,channelIndex)
%% WSVD implementation
% WSVD meant to used in spectral domain but it worked for me in DMI data.
% input multichannel dataset=fiddata(spec,channel,kx,ky,kz)
% NCov=noise covariance matrix-> Make sure you have enough amount of points
% to generate noise covariance matrix
% Channel index=usually 2. Indicate which dimension is for channels.
% Noise whitening
[noiseVec, noiseVal] = eig(NCov);
scaleMatrixFft = noiseVec*diag(sqrt(0.5)./sqrt(diag(noiseVal)));
invScaleMatrixFft = inv(scaleMatrixFft);
%% Channel combination
dims=size(fiddata);
nvox=prod(dims(3:end));
ncoils=dims(channelIndex);
dims(channelIndex)=[];
WSVDcombined=zeros(squeeze(dims));
wsvdQuality=zeros(dims(2:end));
wsvdWeights=zeros([ncoils dims(2:end)]);
Dataformed=reshape(fiddata,[dims(1) ncoils nvox]);
for k=1:nvox
    scaledSpectra=squeeze(Dataformed(:,:,k))*scaleMatrixFft;
    %Optimal reconstruction with SVD
    [u,s,v]=svd(scaledSpectra,'econ');
    wsvdQuality(k) = ((s(1,1)/norm(diag(s)))*sqrt(ncoils)-1)/(sqrt(ncoils)-1);
    
    wsvdCoilAmplitudes=v(:,1)'*invScaleMatrixFft;
    
    svdRescale = norm(wsvdCoilAmplitudes)*normalise(wsvdCoilAmplitudes(1));
    wsvdCoilAmplitudes=wsvdCoilAmplitudes / svdRescale;
    
    
    
    WSVDcombined(:,k) = u(:,1)*s(1,1)*svdRescale;
    wsvdWeights(:,k) = 0.5*pinv(NCov) * wsvdCoilAmplitudes' * conj(svdRescale) * svdRescale;
    
end
end
