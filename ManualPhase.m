function ManualPhase(FID,ppm_axis)

% Gui to play with phase of the spectra with two sliders(0th and 1st)
% Calculatons are done assuming DMI spectra acquired at 5k Hz bandwidth

FigH = figure('position',[50 500 900 900]);
sgtitle('Manual phasing','FontSize',24)
axes('XLim', [0 4*pi], 'units','pixels', ...
    'position',[100 300 700 550], 'NextPlot', 'add');
x     = ppm_axis;
y     = fftshift(fft(FID,[],1),1);
plot(x,abs(y),'LineWidth',1.8)
hold on
LineH = plot(x,real(y),'LineWidth',1.8);
hold off
xlim([0 10]); % Frequency limits
ylim([-max(abs(y),[],'all') 1.5*max(abs(y),[],'all')]);
set(gca,'XDir','reverse')
ax = gca;
ax.XAxis.FontSize=28;
ax.YTick=[];
xlabel('^2H frequency(ppm)','FontSize',24 )
legend('Mag Spectra','Real Spectra','FontSize',22)

annotation(FigH,'textbox',[0.15 0.11 0.1 0.1],'String',{'0^t^h Phase'},'FitBoxToText','on','FontSize',24);
annotation(FigH,'textbox', [0.45 0.11 0.1 0.1],'String',{'1^s^t Phase (by TE)'},'FitBoxToText','on','FontSize',24);

TextH1 = uicontrol('style','text',...
    'position',[100 80 250 50],'FontSize',24);
TextH2 = uicontrol('style','text',...
    'position',[450 80 200 50],'FontSize',24);

SliderH1 = uicontrol('style','slider','position',[50 50 300 20],...
    'min', -pi, 'max', pi);
SliderH2 = uicontrol('style','slider','position',[400 50 300 20],...
    'min', 0.00000, 'max', 0.00800);

addlistener(SliderH1, 'Value', 'PostSet', @callbackfn1);
addlistener(SliderH2, 'Value', 'PostSet', @callbackfn2);


movegui(FigH, 'center')
    function callbackfn1(source, eventdata)
        num0          = get(eventdata.AffectedObject, 'Value');
        LineH.YData  = real(fftshift(fft(FID.*(exp(-1i*num0)),[],1),1).*exp(-1i * (2* pi * (x-4.7).'*45.7 * SliderH2.Value)));
        TextH1.String = strcat(num2str(num0,3),' radians ');
    end

    function callbackfn2(source, eventdata)
        num1          = get(eventdata.AffectedObject, 'Value');
        LineH.YData  = real(fftshift(fft(FID.*(exp(-1i*SliderH1.Value)),[],1),1).*exp(-1i * (2* pi * (x-4.7).'*45.7 * num1)));
        TextH2.String = strcat(num2str(num1*1000,4),' ms');
    end
end
%%
% % Test data
% 
% clear SPECTRA
% clear FID
% 
% Slice=6;
% 
% % TA Scan 1
% FID=squeeze(V9259.LegCSI.RoemerEqualfid(:,3:4,3:4,Slice));
% % TA Scan 2
% FID(:,1,1)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,9,3,Slice));
% FID(:,1,2)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,8,3,Slice));
% FID(:,2,1)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,9,4,Slice));
% FID(:,2,2)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,8,4,Slice));

%%SOL Scan 1
% FID(:,1,1)=squeeze(V9259.LegCSI.RoemerEqualfid(:,4,7,Slice));
% FID(:,1,2)=squeeze(V9259.LegCSI.RoemerEqualfid(:,5,7,Slice));
% FID(:,2,1)=squeeze(V9259.LegCSI.RoemerEqualfid(:,6,6,Slice));
% FID(:,2,2)=squeeze(V9259.LegCSI.RoemerEqualfid(:,7,5,Slice));
%
% % SOL Scan 2
% FID(:,1,1)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,5,5,Slice));
% FID(:,1,2)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,6,6,Slice));
% FID(:,2,1)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,6,7,Slice));
% FID(:,2,2)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,7,7,Slice));

% GM Scan 1
% FID(:,1,1)=squeeze(V9259.LegCSI.RoemerEqualfid(:,5,9,Slice));
% FID(:,1,2)=squeeze(V9259.LegCSI.RoemerEqualfid(:,6,9,Slice));
% FID(:,2,1)=squeeze(V9259.LegCSI.RoemerEqualfid(:,7,8,Slice));
% % FID(:,2,2)=squeeze(V9259.LegCSI.RoemerEqualfid(:,8,8,Slice));
% 
% % GM Scan 2
% % FID(:,1,1)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,5,6,Slice));
% % FID(:,1,2)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,5,7,Slice));
% % FID(:,2,1)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,5,8,Slice));
% % FID(:,2,2)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,6,9,Slice));
% 
% % ApodizationParam=0;
% % FID=FID.*repmat(exp(-ApodizationParam*pi*V9259.LegCSI.Parameters.time).',[1 size(FID,2) size(FID,3)]); % Apodization
% % FID(end+1:2*size(FID,1),:,:,:)=0;
% 
% SPECTRA=fftshift(fft(FID,[],1),1);
% ppm_axis_shift=0.2;
% if size(SPECTRA,1)==1024
%     ppm_axis=V9259.LegCSI.xaxis+ppm_axis_shift; % Shift ppm axis to center water splits on 4.7 ppm
%     
% elseif size(SPECTRA,1)==2048
%     ppm_axis=V9259.LegCSI.xaxiszerofill+ppm_axis_shift; % Shift ppm axis to center water splits on 4.7 ppm
%     
% end
% close all;
% % Function
% % Parameters=V9259.LegCSI.Parameters;
% Parameters=V9259.LegCSI2.Parameters;
% 
% HSVD_QuadrupolarSplittigCalculation(FID(:,1,1),Parameters)%.*exp(-1i*0.25)
% % ManualPhase(FID(:,2,1),ppm_axis)
%