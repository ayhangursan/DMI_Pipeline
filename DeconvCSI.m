function FinalImag=DeconvCSI(spatialfid,Parameters)
input=spatialfid;
CSIdims=Parameters.CSIdims;
inputkspace=fftshift(ifftshift(fft(ifftn(squeeze(input)),[],1)),1); % Spatial to kspace

ZPadkspace=padarray(inputkspace,[0 CSIdims(1) CSIdims(2) CSIdims(3)],0,'both');
ZPadImage=Phasecorrection(ifft(fftn(ZPadkspace),[],1));

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
% Zeropad hamming filter
ZPadidealhammingwindow=padarray(idealhammingwindow,[CSIdims(1) CSIdims(2) CSIdims(3)],0,'both');
ZPadPSF=fftshift(fftn(fftshift(ZPadidealhammingwindow)));

% Deconvolution
DeconImg=zeros(size(MagZPadImage));
for m=1:size(MagZPadImage,1)
DeconImg(m,:,:,:) =deconvwnr(squeeze(MagZPadImage(m,:,:,:)),abs(ZPadPSF),0.01); % Wiener
% concolution falied. Proceeded with blind deconvolution
% Regularization parameter taken from Martijn Froeling doi:10.1002/mrm.28654
% [DeconImg(m,:,:,:),psfr] = deconvblind(squeeze(MagZPadImage(m,:,:,:)),abs(ZPadPSF),10);
end

% Re-add phase
ComplexDecImg=DeconImg.*exp(1i*PhaZPadImage);

% Truncate k-space
DeconZPadKspace=fftshift(ifftshift(fft(ifftn(squeeze(ComplexDecImg)),[],1)),1);
DeconKspace=DeconZPadKspace(:,CSIdims(1)+1:end-CSIdims(1),CSIdims(2)+1:end-CSIdims(2),CSIdims(3)+1:end-CSIdims(3));

% Final image
FinalImag=Phasecorrection(ifft(fftn(DeconKspace),[],1));

%% Test the function
% figure
% subplot(2,3,1)
% imagesc(squeeze(sum(abs(inputkspace(1:15,:,:,ceil(CSIdims(3)/2))))))
% title('input kspace')
% subplot(2,3,2)
% 
% imagesc(squeeze(sum(abs(spatialfid(1:15,:,:,ceil(CSIdims(3)/2))))))
% title('input spatial fid')
% 
% subplot(2,3,3)
% imagesc(squeeze(sum(abs(DeconZPadKspace(1:15,:,:,ceil(3*CSIdims(3)/2))))))
% title('zeropad kspace')
% 
% subplot(2,3,4)
% imagesc(squeeze(sum(abs(DeconKspace(1:15,:,:,ceil(CSIdims(3)/2))))))
% title('kspace')
% 
% subplot(2,3,5)
% imagesc(squeeze(sum(real(FinalImag(:,:,:,ceil(CSIdims(3)/2))))))
% title('Final img')

end