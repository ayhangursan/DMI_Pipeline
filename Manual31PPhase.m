function Manual31PPhase(FID,ppm_axis)
ImagFreq=120.76;
% Gui to play with phase of the spectra with two sliders(0th and 1st)
FigH = figure('position',[50 500 850 850]);
sgtitle('Manual phasing','FontSize',24)
axes('XLim', [0 4*pi], 'units','pixels', ...
    'position',[100 350 600 450], 'NextPlot', 'add');
x     = ppm_axis;
y     = fftshift(fft(FID,[],1),1);
plot(x,abs(y),'LineWidth',1.8)
hold on
LineH = plot(x,real(y),'LineWidth',1.8);
hold off
xlim([-25 25]);
ylim([-max(abs(y),[],'all') 1.5*max(abs(y),[],'all')]);
set(gca,'XDir','reverse')
ax = gca;
ax.XAxis.FontSize=28;
ax.YTick=[];
xlabel('^2H frequency(ppm)','FontSize',24 )
legend('Mag Spectra','Real Spectra','FontSize',24)

annotation(FigH,'textbox',[0.15 0.125 0.1 0.1],'String',{'0^t^h Phase'},'FitBoxToText','on','FontSize',24);
annotation(FigH,'textbox', [0.45 0.125 0.1 0.1],'String',{'1^s^t Phase (by TE)'},'FitBoxToText','on','FontSize',24);

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
        LineH.YData  = real(fftshift(fft(FID.*(exp(-1i*num0)),[],1),1).*exp(-1i * (2* pi * (x).'*ImagFreq * SliderH2.Value)));
        TextH1.String = strcat(num2str(num0,3),' radians ');
    end

    function callbackfn2(source, eventdata)
        num1          = get(eventdata.AffectedObject, 'Value');
        LineH.YData  = real(fftshift(fft(FID.*(exp(-1i*SliderH1.Value)),[],1),1).*exp(-1i * (2* pi * (x).'*ImagFreq * num1)));
        TextH2.String = strcat(num2str(num1*1000,4),' ms');
    end
end
