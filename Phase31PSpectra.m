function PhasedFID=Phase31PSpectra(FID,ppmaxis,Parameters)
%%Function to phase 31P spectra
dimensions=size(FID);
CSIgrid=dimensions(2:end);
NumLines=numel(FID)./size(FID,1);
SPECTRA=fftshift(fft(reshape(FID,[size(FID,1) NumLines]),[],1),1);
FirstOrdPhaseFunct=Parameters.FirstOrdPhaseFunct;
SPECTRA=SPECTRA.*FirstOrdPhaseFunct;
Phasevec=zeros(1,NumLines);
signalwindow=[find(ppmaxis>-0.5 & ppmaxis<1) find(ppmaxis>-1.7 & ppmaxis<-1.9) find(ppmaxis>-6. & ppmaxis<-8) find(ppmaxis>-16.5 & ppmaxis<-14.6) find(ppmaxis>4.8 & ppmaxis<6.4)];%Pcr ATP peaks Pi

for m=1:NumLines %#ok<BDSCI>
    spectrum=SPECTRA(:,m);
    
    
    f = @(ph0)sum(real(spectrum(signalwindow).*(exp(1i*ph0))));
    g = @(ph0)-f(ph0);
    Phaselim0 = 0;
    [optimphase,~] = fminsearch(g,Phaselim0);
    Phasevec(m)=optimphase;
end

PhasedFID=FID.*shiftdim(repmat(exp(1i*reshape(Phasevec,CSIgrid)),[ones(1,size(CSIgrid,2)) size(FID,1)]),size(CSIgrid,2));
%% Debug
% plotLW=1.5;
% figure('WindowState','maximized')
% plot(ppmaxis,abs(spectrum),'k','LineWidth',plotLW)
% hold on
% % plot(ppmaxis,real(fftshift(fft(FID(:,4,7),[],1),1)),'LineWidth',plotLW)
% plot(ppmaxis,real(spectrum),'b','LineWidth',plotLW)
% plot(ppmaxis,real(spectrum.*(exp(1i*xmin))),'r','LineWidth',plotLW)
% % plot(ppmaxis,real(spectrum(real(spectrum.*(exp(1i*xmin)))>0)),'r','LineWidth',plotLW)
%
% hold off
% pbaspect([16 9 1])
% legend('Raw-mag','1^s^tPhased-real','Autophased-real','FontSize',20,'Location','southoutside')
% title('SPECTRA','FontSize',20)
% set(gca,'XDir','reverse')
% ax = gca;
% ax.XAxis.FontSize=28;
% ax.YTick=[];
% xlabel('^3^1P frequency (ppm)','FontSize',24)
