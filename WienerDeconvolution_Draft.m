% function WienerDeconvCSI(inputkspace,Parameters)
% inputkspace=fftshift(fft(ifftshift(ifftn(V9139.BaselineCSI1.RoemerEqualfid)),[],1),1);
input=V9226.Scan1.RoemerEqualfid;
SNR=V9226.Scan1.SNR.Roemer;
time=V9139.BaselineCSI1.Parameters.time;
xaxis=V9226.Scan1.xaxis;
Parameters=V9226.Scan1.Parameters;
% inputkspace=fftshift(ifftshift(fft(ifftn(squeeze(V9139.BaselineCSI1.RoemerEqualfid)),[],1)),1); % Spatial to kspace
CSIdims=size(SNR);
%%
% FIDsignal=repmat((exp(-pi*time*20).*exp(time* 1i*2*pi*0.2*(45.7).')).',[1 CSIdims(1) CSIdims(2) CSIdims(3)]);
% input=zeros(size(FIDsignal));
% input=permute(repmat(SNR,[1 1 1 1024]),[4 1 2 3]).*FIDsignal;
inputkspace=fftshift(ifftshift(fft(ifftn(squeeze(input)),[],1)),1); % Spatial to kspace

%%
figure
subplot(2,2,1)
plot(abs(input(:,5,5,5)))
title('FID-Voxel from center of the image')
subplot(2,2,2)
plot(abs(inputkspace(:,5,5,5)))
title('FID-center of the kspace')

subplot(2,2,3)
imagesc(squeeze(abs(sum(input(:,:,:,5),1))));
title('SNR map')

subplot(2,2,4)
imagesc(squeeze(abs(sum(inputkspace(:,:,:,5),1))));
title('Initial kspace')


ZPadkspace=padarray(inputkspace,[0 CSIdims(1) CSIdims(2) CSIdims(3)],0,'both');
ZPadImage=Phasecorrection(ifft(fftn(ZPadkspace),[],1));

% Show image to check
figure;
subplot(1,2,1)
imagesc(squeeze((abs(ZPadImage(1,:,:,floor(CSIdims(3)*1.5))))));
title('ZPad Image')
daspect([1 1 1])
subplot(1,2,2)
plot(abs(ZPadImage(:,floor(CSIdims(1)*1.5),floor(CSIdims(2)*1.5),floor(CSIdims(3)*1.5))));
title('Fid from center ZPad Image')

MagZPadImage=abs(ZPadImage);
PhaZPadImage=atan2(imag(ZPadImage),real(ZPadImage));

% Calculate ideal hamming filter

for k=1:prod(CSIdims)
    
    [x,y,z]= ind2sub(CSIdims,k);
    filterfuncx=0.54+0.46*cos(2*pi*x/CSIdims(1));
    filterfuncy=0.54+0.46*cos(2*pi*y/CSIdims(2));
    filterfuncz=0.54+0.46*cos(2*pi*z/CSIdims(3));
    hammingwindow(k)=filterfuncx*filterfuncy*filterfuncz;
end

hammingwindow=fftshift(reshape(hammingwindow,CSIdims));

idealhammingwindow=circshift(hammingwindow,1,1);
idealhammingwindow=circshift(idealhammingwindow,1,2);
idealhammingwindow=circshift(idealhammingwindow,1,3);
%% Zeropad hamming filter
ZPadidealhammingwindow=padarray(idealhammingwindow,[CSIdims(1) CSIdims(2) CSIdims(3)],0,'both');
ZPadPSF=fftshift(fftn(fftshift(ZPadidealhammingwindow)));
% Show hamming window to check
figure;
subplot(1,3,1)
imagesc(squeeze((abs(ZPadidealhammingwindow(:,:,floor(CSIdims(3)*1.5))))));
title('Zeropadded hamming filter')
daspect([1 1 1])
colorbar;
subplot(1,3,2)
imagesc(real(ZPadPSF(:,:,floor(CSIdims(3)*1.5))));
title('After FFT')
daspect([1 1 1])
colorbar;

subplot(1,3,3)
plot(real(ZPadPSF(floor(CSIdims(1)*1.5),:,floor(CSIdims(3)*1.5))))
title('1D of the PSF')
%% Deconvolution
DeconImg=zeros(size(MagZPadImage));
for m=1:size(MagZPadImage,1)
% DeconImg(m,:,:,:) =
% deconvwnr(squeeze(MagZPadImage(m,:,:,:)),abs(ZPadPSF),0.007); % Wiener
% concolution falied. Proceeded with blind deconvolution

[DeconImg(m,:,:,:),psfr] = deconvblind(squeeze(MagZPadImage(m,:,:,:)),abs(ZPadPSF),10);
end

% Show deconvolved image for check
figure;
subplot(2,2,1)
imagesc(squeeze(MagZPadImage(1,:,:,floor(CSIdims(3)*1.5))))
title('MagZPadImage')
daspect([1 1 1]);colorbar

subplot(2,2,2)
imagesc(squeeze(abs(ZPadPSF(:,:,floor(CSIdims(3)*1.5)))))
title('ZPadPSF')
daspect([1 1 1]);colorbar

subplot(2,2,3)
imagesc(squeeze(sum(abs(DeconImg(:,:,:,floor(CSIdims(3)*1.5))),1)))
title('DeconImg')
daspect([1 1 1]);colorbar

% Re-add phase
ComplexDecImg=DeconImg.*exp(1i*PhaZPadImage);

subplot(2,2,4)
imagesc(squeeze(sum(abs(ComplexDecImg(:,:,:,floor(CSIdims(3)*1.5))),1)))
title('DeconImg with phase(=same')
daspect([1 1 1]);colorbar
%% Truncate
DeconZPadKspace=fftshift(ifftshift(fft(ifftn(squeeze(ComplexDecImg)),[],1)),1);
DeconKspace=fftshift(DeconZPadKspace(:,CSIdims(1)+1:end-CSIdims(1),CSIdims(2)+1:end-CSIdims(2),CSIdims(3)+1:end-CSIdims(3)));
% Final image
FinalImag=Phasecorrection(ifftshift(ifft(fftshift(fftn(DeconKspace)),[],1)));

figure
subplot(1,3,1)
imagesc(squeeze(sum(abs(DeconZPadKspace(:,:,:,13)),1)));colorbar
title('DeconZPadKspace')
daspect([1 1 1])

subplot(2,3,2)
imagesc(squeeze(sum(abs(inputkspace(:,:,:,5)),1)));colorbar
title('Initital Kspace')
daspect([1 1 1])

subplot(2,3,5)
imagesc(squeeze(sum(abs(DeconKspace(:,:,:,5)),1)));colorbar
title('DeconKspace')
daspect([1 1 1])

subplot(2,3,3)
imagesc(squeeze(sum(abs(input(:,:,:,5)),1)));colorbar
title('Initial Imag')
daspect([1 1 1]);colorbar

subplot(2,3,6)
imagesc(squeeze(sum(abs(FinalImag(:,:,:,5)),1)));colorbar
title('Final Imag')
daspect([1 1 1]);colorbar
%%
row=4;
col=5;
slice=5;
figure
subplot(1,2,1)
% plot(xaxis,real(fftshift(fft(input(:,row,col,slice)))))
plot(abs(((input(:,row,col,slice)))))

% xlim([0 10])
% set(gca,'XDir','reverse')
title('input')

subplot(1,2,2)
plot(abs(((FinalImag(:,row,col,slice)))))
% xlim([0 10])
% set(gca,'XDir','reverse')
title('Final Imag')

%%
figure
subplot(1,2,1)
imagesc(squeeze(sum(abs(input(1:10,:,:,5)),1)));colorbar
title('Input')
daspect([1 1 1])
subplot(1,2,2)
imagesc(squeeze(sum(abs(FinalImag(1:10,:,:,5)),1)));colorbar
title('Final Imag')
daspect([1 1 1])