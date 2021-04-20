function gcaset=timeresolvedspectra(StructDataset,timeaxis,row,column,slice)
% timeaxis should be in  minutes
scansH=fieldnames(StructDataset);
xaxis=eval(strcat('StructDataset.',string(scansH(1)),'.xaxiszerofill;')); % With zerofill
% xaxis=eval(strcat('StructDataset.',string(scansH(1)),'.xaxis;'));

for n=1:numel(scansH)
    %     spectratime(:,n)=eval(strcat('StructDataset.',string(scansH(n)),'.FinalspectraRoemer_aligned(:,',num2str(row),',',num2str(column),',',num2str(slice),');'));disp('Using aligned spectra');
    %     spectratime(:,n)=eval(strcat('StructDataset.',string(scansH(n)),'.FinalspectraRoemer(:,',num2str(row),',',num2str(column),',',num2str(slice),');'));

    %         disp('Multi coil acquisition dataset.')
    %         fidtime=eval(strcat('StructDataset.',string(scansH(n)),'.RoemerEqualfid(:,',num2str(row),',',num2str(column),',',num2str(slice),').*exp(-5*pi*StructDataset.',string(scansH(n)),'.Parameters.time).'';'));
    %         fidtime(end:2048)=0;
    %         spectraNoApod(:,n)=fftshift(fft(fidtime)).*exp(-1i * (2* pi * (xaxis-4.7).'*45.7 * 0.00195));
    
    disp('Single coil dataset.')
    fidtime=eval(strcat('StructDataset.',string(scansH(n)),'.fftfiddata(:,',num2str(row),',',num2str(column),',',num2str(slice),').*exp(-5*pi*StructDataset.',string(scansH(n)),'.Parameters.time).'';'));
    fidtime(end:2048)=0;
    spectraNoApod(:,n)=fftshift(fft(fidtime)).*exp(-1i * (2* pi * (xaxis-4.7).'*45.7 * 0.0016));
    
    %     spectratime(:,n)=eval(strcat('StructDataset.',string(scansH(n)),'.combinedspectradata(:,',num2str(row),',',num2str(column),',',num2str(slice),');'));
    
end
% [h,m,s] = hms(hours(timeaxis(:,1))+minutes(timeaxis(:,2))+seconds(timeaxis(:,3)));
% timelabels=strcat(num2str(h(:),'%02.f'),':',num2str(m(:),'%02.f'),':',num2str(s(:),'%02.f'));
ylabeljump=3;
ylabelstart=3;
timelabels=strcat(num2str(timeaxis(ylabelstart:ylabeljump:end),'%02.f'));


figure('units','normalized','outerposition',[0 0 1 1],'WindowState', 'maximized');
fig1=waterfall(xaxis,[1:numel(timeaxis)],real(spectraNoApod).');
axes1 = fig1.Parent;
set(gca,'XDir','reverse','ZtickLabel',[],'LineWidth',2,'Colormap',[0 0 0]);%,'Colormap',[hot(1);hot(1);hot(1)]
ax=gca;
ax.ZAxis.Color='b';
% fig1=waterfall(xaxis,timeaxis,real(spectratime).');
yticklabels(timelabels);
% set(gca,'Yticklabel','Fontsize',24)
ylabel('Time (min)','FontSize',28,'Rotation',40);
xlabel('^2H frequency (ppm)','HorizontalAlignment','right','FontSize',28);

yticks([ylabelstart:ylabeljump:numel(scansH)]);

xlim([0 8]); % ppm axis limits

% define min and max for spectra to avoid offset
zlimitwindow=find(xaxis>0 & xaxis<10);
zlow=min(real(spectraNoApod(zlimitwindow,:)),[],'all');
zhigh=max(real(spectraNoApod(zlimitwindow,:)),[],'all');
zlim([zlow zhigh]);

view([24.1875 43.3418]);
fig1.LineWidth=1.8;
title(strcat('Voxel-Row:',num2str(row),' Column:',num2str(column),' Slice:',num2str(slice)),'FontSize',24);
% set(axes1,'FontSize',20,'FontWeight','bold','YTick',1:5:numel(scansH),'YTickLabel',{timelabels(1:5:numel(scansH),:)});
% set(axes1,'FontSize',20,'FontWeight','bold','XDir','reverse','YTick',[1 8 14 20 24],'YTickLabel',{'0','30','60','90','120'})
set(axes1,'FontSize',24,'FontWeight','bold','XDir','reverse','ZTick',[],'YGrid','off','XGrid','off','ZColor',[1 1 1],'Color',[1 1 1])

gcaset=gca;
end