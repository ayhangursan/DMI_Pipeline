function ShowAmaresfit(inputspec,DMIResults,row,col,slice,xaxis,Parameters,Referenceloc)

% inputfid=inputstruct.RoemerEqualfid;

    showrows=1;
fitResults=DMIResults{row,col,slice};
time=Parameters.time;
Freqshift=(Referenceloc(sub2ind(Parameters.CSIdims, row, col, slice))-4.7);
for metabnum=1:numel(fitResults.amplitude)
    lineshapes(:,metabnum)=exp( -abs(time)*fitResults.linewidth(metabnum).' * pi).*exp(-pi^2*time.^2*(fitResults.sigma(metabnum).').^2);
    Resultfid(:,metabnum)=Phasecorrection((fitResults.amplitude(metabnum).*lineshapes(:,metabnum).*exp(2*pi*1i*(fitResults.chemShift(metabnum)+Freqshift)*45.7*time).'.*exp(1i*fitResults.phase(metabnum))));
end


figure('WindowState','maximized')

subplot(showrows,3,1)
plot(xaxis,real(inputspec(:,row,col,slice)))
title('Input data','FontSize',22)
xlim([0 10])
set(gca,'XDir','reverse')
pbaspect([16 9 1])

subplot(showrows,3,2)
plot(xaxis,real(fftshift(fft(Resultfid))))
title('Fit','FontSize',22)
xlim([0 10])
set(gca,'XDir','reverse')
pbaspect([16 9 1])

subplot(showrows,3,3)
plot(xaxis,real(inputspec(:,row,col,slice)))
hold on;
plot(xaxis,real(fftshift(fft(Resultfid))))
hold off;
title('Spectra and Fit','FontSize',22)
xlim([0 10])
pbaspect([16 9 1])
set(gca,'XDir','reverse')

% if isempty(inputstruct.WSVD.WSVDcombined)==0
%     
%     inputfid=inputstruct.WSVD.WSVDcombined;
%     
%     subplot(showrows,3,4)
%     plot(xaxis,real(fftshift(fft(inputfid(:,row,col,slice)))))
%     title('Input data-WSVD','FontSize',22)
%     xlim([0 10])
%     set(gca,'XDir','reverse')
%     pbaspect([16 9 1])
%     
%     subplot(showrows,3,5)
%     plot(xaxis,real(fftshift(fft(Resultfid))))
%     title('Fit-WSVD','FontSize',22)
%     xlim([0 10])
%     set(gca,'XDir','reverse')
%     pbaspect([16 9 1])
%     
%     subplot(showrows,3,6)
%     plot(xaxis,real(fftshift(fft(inputfid(:,row,col,slice)))))
%     hold on;
%     plot(xaxis,real(fftshift(fft(Resultfid))))
%     hold off;
%     title('Spectra and Fit-WSVD','FontSize',22)
%     xlim([0 10])
%     pbaspect([16 9 1])
%     set(gca,'XDir','reverse')
% end

end