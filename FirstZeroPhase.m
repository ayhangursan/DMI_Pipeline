function FirstZeroPhase(FID,ppmaxis,Parameters)
%Function to phase spectra
dimensions=size(FID);
CSIgrid=dimensions(2:end);
ImagingFreq=45.7;
TE=Parameters.TE
waterwindow=[3.8 5.6];
waterind=ppmaxis>waterwindow(1) & ppmaxis<waterwindow(2);
SPECTRA=fftshift(fft(FID,[],1),1);

FirstOrdPhaseFunct=exp(-1i * (2* pi * (ppmaxis-4.7).'*ImagingFreq * TE));

SPECTRA=SPECTRA.*FirstOrdPhaseFunct;
AP=2;RL=2;
spectrum=SPECTRA(:,AP,RL);

f = @(ph0)sum(real(spectrum(waterind,:).*(exp(-1i*ph0))));
g = @(ph0)-f(ph0);
Phaselim0 = [0];
[xmin,gmin] = fminsearch(g,Phaselim0)


%%[SPECTRA PhasedFID]=
plotLW=1.5;
figure('WindowState','maximized')
plot(ppmaxis,abs(spectrum),'LineWidth',plotLW)
hold on
% plot(ppmaxis,real(fftshift(fft(FID(:,4,7),[],1),1)),'LineWidth',plotLW)
plot(ppmaxis,real(spectrum),'LineWidth',plotLW)
plot(ppmaxis,real(spectrum.*(exp(1i*xmin))),'k','LineWidth',plotLW)
% plot(ppmaxis,real(fftshift(fft(FID(:,4,7).*(exp(1i*xmin)),[],1),1)),'LineWidth',plotLW)

hold off
xlim([-5 15])
pbaspect([16 9 1])
legend('Raw-mag','1^s^tPhased-real','Autophased-real','FontSize',20,'Location','southoutside')
title('SPECTRA','FontSize',20)
set(gca,'XDir','reverse')
ax = gca;
ax.XAxis.FontSize=28;
ax.YTick=[];
xlabel('^2H frequency (ppm)','FontSize',24)
% 'Raw-real'
% 'Raw-0^t^h rephased',