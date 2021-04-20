function ShowDMIAmaresfit(inputspec,DMIResults,row,col,slice,xaxis,Parameters,Referenceloc)

% inputfid=inputstruct.RoemerEqualfid;

fitResults=DMIResults{row,col,slice};
time=Parameters.time;
Freqshift=(Referenceloc(sub2ind(Parameters.CSIdims, row, col, slice))-4.7);
for metabnum=1:numel(fitResults.amplitude)
    lineshapes(:,metabnum)=exp( -abs(time)*fitResults.linewidth(metabnum).' * pi).*exp(-pi^2*time.^2*(fitResults.sigma(metabnum).').^2);
    Resultfid(:,metabnum)=Phasecorrection((fitResults.amplitude(metabnum).*lineshapes(:,metabnum).*exp(2*pi*1i*(fitResults.chemShift(metabnum))*45.7*time).'.*exp(1i*fitResults.phase(metabnum))));
end


rawspectra=inputspec(:,row,col,slice).*exp(-1i * (2* pi * (xaxis-4.7).'*45.7 * Parameters.TE));
% rawspectra=inputspec(:,row,col,slice); % NO 1st order phase correction

plotshift=max(real(inputspec(:,row,col,slice)));

figure('WindowState','maximized')

% plot(xaxis,real(inputspec(:,row,col,slice)),'LineWidth',1.25)
plot(xaxis-Freqshift,real(rawspectra),'LineWidth',1.25)

hold on;

plot(xaxis-Freqshift,real(fftshift(fft(sum(Resultfid,2)))),'LineWidth',1.25)
plot(xaxis-Freqshift,real(rawspectra)-real(fftshift(fft(sum(Resultfid,2)))),'LineWidth',1.25)

for metabnum=1:size(Resultfid,2)
    
    plot(xaxis,real(fftshift(fft(Resultfid(:,metabnum))))-plotshift*(metabnum),'LineWidth',1.25)
    
end
hold off
% pbaspect([16 9 1])
pbaspect([1 1 1])
title(strcat('1storder phase correction-Voxel row: ',num2str(row),' column: ',num2str(col),' slice: ',num2str(slice)),'FontSize',22)

xlim([-2 10])
xlabel('2H frequency (ppm)','FontSize',24)
set(gca,'XDir','reverse')
legend1=legend('Raw data','Fitted spectrum','Residual','Lipid','Glucose','Water','FontSize',24);
set(legend1,'Location','eastoutside','FontSize',24);
end