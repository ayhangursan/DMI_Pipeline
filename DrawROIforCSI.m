function [CSIMask]=DrawROIforCSI(dicominput,FOV,spectrainput)
% Draw ROI on anatomical images to generate a mask for CSI grid.
% As an Input T1w images are used to draw an ROI
%% Debug
cd('C:\Users\agursan.DS\surfdrive\Deuterium_CSI\Exams\2020-06-15 2H in-vivo V9073\V9073\DICOM')
dicominput='IM_0012';
% spectrainput=V9073Series.RoemerEqualfid;
inputCSI=size(spectrainput);
inputCSI=inputCSI(2:end);
FOV=[240 300];
% End Debug
CSIMask.DicomImage=squeeze(double(dicomread(dicominput)));
CSIMask.DicomTagInfo=dicominfo(dicominput);
NoSlices=size(CSIMask.DicomImage,3);
Midslice=ceil(NoSlices/2);

% In dicom files from multix system dicom files are stqarting from the mid
% slice. So stack of slices should be reordered to start from the Feet and move
% to the Head.
OrderedSlices=[Midslice+1:NoSlices 1:Midslice];
CSIMask.DicomImage=CSIMask.DicomImage(:,:,OrderedSlices);% Reorder linearly
% Crop with respect to FoV
imageratio=FOV(1)/FOV(2);
ImageSize=[size(CSIMask.DicomImage,1) size(CSIMask.DicomImage,2)];
Filledpart=ImageSize(1)*(1-imageratio);
CSIMask.DicomImage=squeeze(CSIMask.DicomImage(ceil(Filledpart/2):(ImageSize(1)-floor(Filledpart/2)),:,:));
%
figure('WindowState','maximized')
CSIMask.ax(1)=subplot(1,1,1);
imagesc(CSIMask.DicomImage(:,:,1))
daspect([1 1 1])
set(gca, 'Layer','top')
axis on;
[rows, columns, numberOfColorChannels] = size(CSIMask.DicomImage(:,:,1));
hold on;
for row = 1 : rows ./ inputCSI(1) : rows
    line([1, columns], [row, row], 'Color', 'r');
end
for col = 1 : columns ./ inputCSI(2) : columns
    line([col, col], [1, rows], 'Color', 'r');
end
colormap(CSIMask.ax(1),gray)


% Slider
SliderH1 = uicontrol('style','slider','position',[250 50 300 20],...
    'min',1, 'max', size(CSIMask.DicomImage,3),'Tag','slider1', 'Value',1,'SliderStep',[1/(size(CSIMask.DicomImage,3)-1) 0.2]);

TextH1 = uicontrol('style','text',...
    'position',[320 80 200 50],'FontSize',24,'String',strcat('Slice: ',num2str(SliderH1.Value)));

addlistener(SliderH1, 'Value', 'PostSet', @callbackfn1);

%% Start drawing ROI
ROICheck = uicontrol('Style','check','position',[800 50 300 20],'Value',0);
TextH2 = uicontrol('style','text','position',[700 80 200 50],'FontSize',20,'String',strcat('Draw ROI'));

addlistener(ROICheck, 'Value', 'PostSet', @callbackfn2);

%% Stop drawing ROI
ROIFinishCheck = uicontrol('Style','check','position',[1150 50 300 20],'Value',0);
TextH3 = uicontrol('style','text','position',[1050 80 300 50],'FontSize',20,'String',strcat('ROI completed'));

addlistener(ROIFinishCheck, 'Value', 'PostSet', @callbackfn3);

%% Save and Close GUI
CloseGUI = uicontrol('Style','pushbutton','position',[1450 50 100 20],'Value',0);
TextH4 = uicontrol('style','text','position',[1450 80 100 50],'FontSize',20,'String',strcat('Close'));

addlistener(CloseGUI, 'Value', 'PostSet', @callbackfn4);

    function callbackfn1(source, eventdata)
        slice          = round(get(eventdata.AffectedObject, 'Value'));
        
        ax(1)=subplot(1,1,1);
        imagesc(CSIMask.DicomImage(:,:,slice))
        daspect([1 1 1])
        set(gca, 'Layer','top')
        axis on;
        hold on;
        for row = 1 : rows ./ inputCSI(1) : rows
            line([1, columns], [row, row], 'Color', 'r');
        end
        for col = 1 : columns ./ inputCSI(2) : columns
            line([col, col], [1, rows], 'Color', 'r');
        end
        TextH1.String = strcat('Slice: ',num2str(slice));
                if ROICheck.Value==1
                    children = get(gca, 'children');
                    delete(children(1));
                    CSIMask.ax(1).Children.CData = CSIMask.DicomImage(:,:,slice);
                    TextH1.String = strcat('Slice: ',num2str(slice));
        

                elseif ROICheck.Value==0
                    CSIMask.ax(1).Children.CData = CSIMask.DicomImage(:,:,slice);
                    TextH1.String = strcat('Slice: ',num2str(slice));
                end
        
        
    end

    function callbackfn2(source, eventdata)
        DrawROI          = boolean(get(eventdata.AffectedObject, 'Value'));
        if DrawROI==1
            CSIMask.DrawnROI=drawfreehand(CSIMask.ax(1),'Deletable',1);
            CSIMask.DrawnROImask(:,:,SliderH1.Value)=createMask(CSIMask.DrawnROI);
            numel(find(createMask(CSIMask.DrawnROI)))

        else
            return
        end
    end

    function callbackfn3(source, eventdata)
        StopDrawROI          = boolean(get(eventdata.AffectedObject, 'Value'));
        if StopDrawROI==1
            ROICheck.Value=0;
            CSIMask.test=5;
            disp(strcat('Number of pixels included:',num2str(numel(find(CSIMask.DrawnROImask)))))
            std(CSIMask.B0ImagesScaled(CSIMask.DrawnROImask),'omitnan')
            %                                     annotation
            children = get(CSIMask.ax(1), 'children');
            if exist('children','Var')==1
            delete(children(1));
            end
            CSIMask.ax(1).Children.CData = CSIMask.DicomImage(:,:,SliderH1.Value);
%             CSIMask.ax(2).Children.CData = CSIMask.B0ImagesScaled(:,:,SliderH1.Value);

        end
    end

    function callbackfn4(source, eventdata)
        CloseGUI          = boolean(get(eventdata.AffectedObject, 'Value'));
%         save ('DrawnMask.mat','B0Results.DrawnROImask')
%         whos('-file','DrawnMask.mat')
        close all;
    end
%% Debug

end
% cd('C:\Users\agursan.DS\surfdrive\Deuterium_CSI\Exams\2020-12-04 2H in-vivo V9259 Leg\Multix\DICOM\')
