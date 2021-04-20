function CompareT2star(inputdataset)
% MaxSNR=max([inputdataset.SNRRoemer inputdataset.SNRWSVD inputdataset.SNRRoemerDeconv inputdataset.SNRWSVDDeconv],[],'all')*0.9;
figure('WindowState','maximized');
sgtitle('T2 star(ms)')

for m=ceil(inputdataset.Parameters.CSIdims(3)/2)
    subplot(2,2,1)
    imagesc(inputdataset.RoemerT2starmap(:,:,m))
    colorbar
    caxis([0 20])
    daspect([1 1 1])
    title('Roemer eq noise')
    
    subplot(2,2,2)
    imagesc(inputdataset.WSVDT2starmap(:,:,m))
    colorbar
    caxis([0 20])
    daspect([1 1 1])
    title('WSVD')
    
    subplot(2,2,3)
    imagesc(inputdataset.RoemerDeconvT2starmap(:,:,m))
    colorbar
    caxis([0 20])
    daspect([1 1 1])
    title('Roemer eq noise-Deconv')
    
    subplot(2,2,4)
    imagesc(inputdataset.WSVDDeconvT2starmap(:,:,m))
    colorbar
    caxis([0 20])
    daspect([1 1 1])
    title('WSVD-Deconv')
end