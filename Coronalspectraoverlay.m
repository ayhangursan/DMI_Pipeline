function Coronalspectraoverlay(spectra,coronalslice,xaxis)
dims=size(spectra);
spectra=flip(spectra,4);
row=dims(4);
col=dims(3);
figx=figure('units','normalized','outerposition',[0 0 1 1],'WindowState', 'maximized');
set(gcf,'PaperPositionMode','auto', 'Color', 'k');
set(gca, 'Color', 'k');
sgtitle(['Coronal Slice: ',num2str(coronalslice)],'FontSize',24,'Color','red');
globalmin=min((real(spectra)),[],'all');
globalmax=max((real(spectra)),[],'all');

n=1;
for k=1:row
    for j=1:col
        subplot(row,col,n);
        plot(xaxis,squeeze((real(spectra(:,coronalslice,j,k)))),'-g','LineWidth',3), axis([-0.3 9.7 globalmin globalmax]),box off,axis off;
        pbaspect([16 9 1]);
        n = n+1;
        set(gca,'XDir','reverse');
        
    end
end
%  print(gcf, '-dtiff', '-r600', strcat('Coronal slice:',num2str(coronalslice)));%Dont forget to comment it later

% for k=1:row
%     for j=1:col
%         subplot(row,col,n);
%         plot(xaxis,squeeze((real(spectra(:,coronalslice,j,k)))),'-g','LineWidth',3), axis([min(xaxis) max(xaxis) globalmin globalmax]),box off,axis off;
%         pbaspect([16 9 1]);
%         n = n+1;
%         set(gca,'XDir','reverse');
%         
%     end
% end
% disp('Displaying whole spectrum')