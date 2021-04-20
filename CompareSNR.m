function CompareSNR(inputdataset)
MaxSNR=max([inputdataset.SNRRoemer inputdataset.SNRWSVD inputdataset.SNRRoemerDeconv inputdataset.SNRWSVDDeconv],[],'all')*0.9;
figure('WindowState','maximized');
sgtitle('SNR Comparsion','FontSize',24)

for m=ceil(inputdataset.Parameters.CSIdims(3)/2)
    subplot(2,2,1)
    imagesc(inputdataset.SNRRoemer(:,:,m))
    colorbar
    caxis([0 MaxSNR])
    daspect([1 1 1])
    title('Roemer eq noise','FontSize',24)
    
    subplot(2,2,2)
    imagesc(inputdataset.SNRWSVD(:,:,m))
    colorbar
    caxis([0 MaxSNR])
    daspect([1 1 1])
    title('WSVD','FontSize',24)
    
    subplot(2,2,3)
    imagesc(inputdataset.SNRRoemerDeconv(:,:,m))
    colorbar
    caxis([0 MaxSNR])
    daspect([1 1 1])
    title('Roemer eq noise-Deconv','FontSize',24)
    
    subplot(2,2,4)
    imagesc(inputdataset.SNRWSVDDeconv(:,:,m))
    colorbar
    caxis([0 MaxSNR])
    daspect([1 1 1])
    title('WSVD-Deconv','FontSize',24)
end