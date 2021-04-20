%% Draw masks Bended
function [B0Results]=DrawROIonB0(dicominput)
B0Results.DicomImage=double(dicomread(dicominput));
B0Results.DicomTagInfo=dicominfo(dicominput);
NoSlices=size(B0Results.DicomImage,4)/2;
Midslice=ceil(NoSlices/2);
OrderedSlices=[Midslice+1:NoSlices 1:Midslice];
% DrawnROImasksave=1;

RI = double(B0Results.DicomTagInfo.PerFrameFunctionalGroupsSequence.(['Item_',num2str(2)]).Private_2005_140f.Item_1.RescaleIntercept);
% Scale slope - RS
RS = double(B0Results.DicomTagInfo.PerFrameFunctionalGroupsSequence.(['Item_',num2str(2)]).Private_2005_140f.Item_1.RescaleSlope);

%% Flip and rotate if needed
% B0Results.DicomImage=flip(B0Results.DicomImage,1);
% B0Results.DicomImage=flip(B0Results.DicomImage,2);
% DON'T FORGET TO REVERSE
%% Seperate anatomical and B0 map images
B0Results.AnatomicalImages=squeeze(B0Results.DicomImage(:,:,1,1:2:end));
B0Results.AnatomicalImages=B0Results.AnatomicalImages(:,:,OrderedSlices);% Reorder linearly

B0Images=squeeze(B0Results.DicomImage(:,:,1,2:2:end));
B0Images=B0Images(:,:,OrderedSlices);% Reorder linearly
B0Results.B0ImagesScaled=(B0Images.*RS+RI);


figROI=figure('WindowState','maximized');
sgtitle('Draw ROI','FontSize',24)

B0Results.ax(1)=subplot(1,2,1);
imagesc(B0Results.AnatomicalImages(:,:,1))
daspect([1 1 1]);
title('Axial anatomical','FontSize',24)
colormap(B0Results.ax(1),gray)

B0Results.ax(2)=subplot(1,2,2);
imagesc(B0Results.B0ImagesScaled(:,:,1))
daspect([1 1 1]);
title('B0 map','FontSize',24)
colormap(B0Results.ax(2),jet)
colorbar
caxis([RI abs(RI)]);

%% Slider
SliderH1 = uicontrol('style','slider','position',[250 50 300 20],...
    'min',1, 'max', size(B0Results.AnatomicalImages,3),'Tag','slider1', 'Value',1,'SliderStep',[1/(size(B0Results.AnatomicalImages,3)-1) 0.1]);

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


% B0Results.DrawnROI=drawfreehand(B0Results.ax(1),'Deletable',1);
% B0Results.DrawnROImask(:,:,SliderH1.Value)=createMask(B0Results.DrawnROI);

    function callbackfn1(source, eventdata)
        slice          = round(get(eventdata.AffectedObject, 'Value'));
        
        if ROICheck.Value==1
            children = get(gca, 'children');
            delete(children(1));
            B0Results.ax(1).Children.CData = B0Results.AnatomicalImages(:,:,slice);
            B0Results.ax(2).Children.CData = B0Results.B0ImagesScaled(:,:,slice);
            TextH1.String = strcat('Slice: ',num2str(slice));
            
            B0Results.DrawnROI=drawfreehand(B0Results.ax(1),'Deletable',1);
            B0Results.DrawnROImask(:,:,SliderH1.Value)=createMask(B0Results.DrawnROI);
            
        elseif ROICheck.Value==0
            B0Results.ax(1).Children.CData = B0Results.AnatomicalImages(:,:,slice);
            B0Results.ax(2).Children.CData = B0Results.B0ImagesScaled(:,:,slice);
            TextH1.String = strcat('Slice: ',num2str(slice));
        end
        
        
    end

    function callbackfn2(source, eventdata)
        DrawROI          = boolean(get(eventdata.AffectedObject, 'Value'));
        if DrawROI==1
            B0Results.DrawnROI=drawfreehand(B0Results.ax(1),'Deletable',1);
            B0Results.DrawnROImask(:,:,SliderH1.Value)=createMask(B0Results.DrawnROI);
            numel(find(createMask(B0Results.DrawnROI)))
            
        else
            return
        end
    end

    function callbackfn3(source, eventdata)
        StopDrawROI          = boolean(get(eventdata.AffectedObject, 'Value'));
        if StopDrawROI==1
            ROICheck.Value=0;
            B0Results.test=5;
            
            %             B0Results.DrawnROImasksave=DrawnROImask;
            %                         global ROIMask
            %                         ROIMask=DrawnROImask;
            % %             Calculate Std(mask)
            disp(strcat('Number of pixels included:',num2str(numel(find(B0Results.DrawnROImask)))))
            std(B0Results.B0ImagesScaled(B0Results.DrawnROImask),'omitnan')
            %                                     annotation
            children = get(B0Results.ax(1), 'children');
            if exist('children','Var')==1
            delete(children(1));
            end
            B0Results.ax(1).Children.CData = B0Results.AnatomicalImages(:,:,SliderH1.Value);
            B0Results.ax(2).Children.CData = B0Results.B0ImagesScaled(:,:,SliderH1.Value);
            
        end
    end

    function callbackfn4(source, eventdata)
        CloseGUI          = boolean(get(eventdata.AffectedObject, 'Value'));
%         save ('DrawnMask.mat','B0Results.DrawnROImask')
%         whos('-file','DrawnMask.mat')
        close all;
    end
end
%% Debug
% cd('C:\Users\agursan.DS\surfdrive\Deuterium_CSI\Exams\2020-12-04 2H in-vivo V9259 Leg\Multix\DICOM\')
