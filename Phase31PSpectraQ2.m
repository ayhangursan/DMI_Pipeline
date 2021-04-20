function PhasedFIDQ2=Phase31PSpectraQ2(FID,Parameters)
NumLines=numel(FID)./size(FID,1);
SPECTRA=reshape(fftshift(fft(FID,[],1),1),[size(FID,1) NumLines]);
FirstOrdPhaseFunct=Parameters.FirstOrdPhaseFunct;
SPECTRA=SPECTRA.*FirstOrdPhaseFunct;
PhasedSPECTRA=zeros(size(SPECTRA));
for kk=1:NumLines
    spectrum=SPECTRA(:,kk);
    phasevector=angle(spectrum);
    [~, ind]=max(abs(spectrum));
    PhasedSPECTRA(:,kk)=complex(abs(spectrum).*cos(phasevector-phasevector(ind)) ,abs(spectrum).*sin(phasevector-phasevector(ind)) );
end
PhasedSPECTRA=PhasedSPECTRA./FirstOrdPhaseFunct;
PhasedSPECTRA=reshape(PhasedSPECTRA,size(FID));
PhasedFIDQ2=ifft(ifftshift(PhasedSPECTRA,1),[],1);
