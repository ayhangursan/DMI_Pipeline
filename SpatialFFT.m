function LocalizedSignal=SpatialFFT(KSpaceSignal)
% Function to use for MR Spectroscopy data.
%  Input data order (Num Points,k1,k2,k3)
% Be careful on extra domains! Such as channels or dynamics
    LocalizedSignal=ifftshift(fftshift(ifft(fftn(squeeze(KSpaceSignal)),[],1)),1);
