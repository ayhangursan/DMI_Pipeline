function spectraoverlay(spectra,slice,xaxis,globalmax)

dims=size(spectra);
row=dims(2);
col=dims(3);
figx=figure('units','normalized','outerposition',[0 0 1 1],'WindowState', 'maximized');
set(gcf,'PaperPositionMode','auto', 'Color', 'k');
set(gca, 'Color', 'k');

sgtitle(['Transverse Slice: ',num2str(slice)],'FontSize',24,'Color','red');
globalmin=min((real(spectra)),[],'all');
if nargin<4
    globalmax=max((real(spectra)),[],'all')
end
n=1;
for k=1:row
    for j=1:col
        subplot(row,col,n);
        plot(xaxis,squeeze((real(spectra(:,k,j,slice)))),'-g','LineWidth',3), axis([-1.3 10.7 globalmin globalmax]),box off,axis off;
        pbaspect([16 9 1]);
        n = n+1;
        set(gca,'XDir','reverse');
    end
end
% print(gcf, '-dtiff', '-r600', 'insert name here later');
