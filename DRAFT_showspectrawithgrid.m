
clear SPECTRA
clear FID
% SPECTRA=squeeze(V9139.Glucose22.combinedspectradata(:,:,:,5));
% SPECTRA=squeeze(V9245.Glucose8.combinedspectradata(:,:,:,7));
% SPECTRA=flip(permute(squeeze(V9259.LegCSI2.combinedspectradata(:,5,5,6)),[1 3 2]),3);

% All grid
% Scan 1
% FID=squeeze(V9259.LegCSI.RoemerEqualfid(:,:,:,6));
% Scan 2
% FID=flip(permute(squeeze(V9259.LegCSI2.RoemerEqualfid(:,3:end,1:11,6)),[1 3 2]),3);

% %%Green boxes
% Scan 1

% % % SPECTRA=squeeze(V9259.LegCSI.combinedspectradata(:,3:4,3:4,6));
% FID=squeeze(V9259.LegCSI.RoemerEqualfid(:,3:4,3:4,6));

% Scan 2
% FID(:,1,1)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,9,3,6));
% FID(:,1,2)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,8,3,6));
% FID(:,2,1)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,9,4,6));
% FID(:,2,2)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,8,4,6));

% %%Blue boxes
% Scan 1
 FID(:,1,1)=squeeze(V9259.LegCSI.RoemerEqualfid(:,4,7,6));
 FID(:,1,2)=squeeze(V9259.LegCSI.RoemerEqualfid(:,5,7,6));
 FID(:,2,1)=squeeze(V9259.LegCSI.RoemerEqualfid(:,6,6,6));
 FID(:,2,2)=squeeze(V9259.LegCSI.RoemerEqualfid(:,7,5,6));

% Scan 2
% FID(:,1,1)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,5,5,6));
% FID(:,1,2)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,6,6,6));
% FID(:,2,1)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,6,7,6));
% FID(:,2,2)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,7,7,6));

% %%Red boxes
% Scan 1
%  FID(:,1,1)=squeeze(V9259.LegCSI.RoemerEqualfid(:,5,9,6));
%  FID(:,1,2)=squeeze(V9259.LegCSI.RoemerEqualfid(:,6,9,6));
%  FID(:,2,1)=squeeze(V9259.LegCSI.RoemerEqualfid(:,7,8,6));
%  FID(:,2,2)=squeeze(V9259.LegCSI.RoemerEqualfid(:,8,8,6));

% Scan 2
% FID(:,1,1)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,5,6,6));
% FID(:,1,2)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,5,7,6));
% FID(:,2,1)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,5,8,6));
% FID(:,2,2)=squeeze(V9259.LegCSI2.RoemerEqualfid(:,6,9,6));
% % % 


% For DEPSI
% SPECTRA=flip(fftshift(fft(DEPSI.T2(:,:,:,8),[],1),1),1);
% ppm_axis=DEPSI.ppmaxis;

% Apodization and zerofilling
FID=FID.*repmat(exp(-2*pi*V9259.LegCSI.Parameters.time).',[1 size(FID,2) size(FID,3)]); % Apodization
FID(end+1:2*size(FID,1),:,:,:)=0;
SPECTRA=fftshift(fft(FID,[],1),1);


if size(SPECTRA,1)==1024
        ppm_axis=V9259.LegCSI.xaxis;

elseif size(SPECTRA,1)==2048
        ppm_axis=V9259.LegCSI.xaxiszerofill;

end

% SPECTRA=squeeze(specdenoisedfid3(:,4:7,3:6,7));
% SPECTRA=squeeze(fftshift(fft(V9245.Glucose8.LPfidRoemer(:,:,:,7),[],1),1));
% SPECTRA=squeeze(V9259.LegCSI2.combinedspectradata(:,2:4,2:4,7));
% SPECTRA=flip(permute(squeeze(V9259.LegCSI2.combinedspectradata(:,3:end,1:11,6)),[1 3 2]),3);
% SPECTRA=flip(permute(squeeze(V9259.LegCSI2.combinedspectradata(:,4:6,5:7,6)),[1 3 2]),3);

% SPECTRA=squeeze(V9247.HighRes.FinalspectraDenoised(:,:,:,7));

% SPECTRA=flip(permute(squeeze(V9247.LowRes.FinalspectraRoemer(:,8,3:7,4:10)),[1
% 3 2]),2); %Coronal view
% SPECTRA=squeeze(V9245.Glucose1.FinalspectraRoemer(:,:,:,7));
% SPECTRA=flip(permute(squeeze(V9245.Glucose11.FinalspectraRoemer(:,5,:,:)),[1 3 2]),2); %Coronal view
%
% SPECTRAthreshold=SPECTRA;
% cutoffvalue=6e+04;
% SPECTRA(find(SPECTRAthreshold > cutoffvalue))=cutoffvalue;
max_peak1=max(real(SPECTRA),[],'all');
specmax=max(real(SPECTRA),[],'all');
% specmax=500;% Over run scaling for scan1
% specmax=860;% Over run scaling for scan2

zlimitwindow=find(ppm_axis>0 & ppm_axis<10);
specmin=min(real(SPECTRA(zlimitwindow,:,:)),[],'all');
% min=min(real(SPECTRA),[],'all');
nkx=size(SPECTRA,2);
nky=size(SPECTRA,3);

% figure
% cnumber = 0;
% for idx_ix = 1:nkx
%     for idx_iy = 1:nky
%         cnumber = cnumber + 1;
%         subplot(nkx,nky,cnumber);
%         plot(ppm_axis, real(squeeze(SPECTRA(:,idx_ix,idx_iy))),'b');
%         set(gca,'XDir','reverse'); box off;
%         h = gca; h.YAxis.Visible = 'off'; h.XAxis.Visible = 'off';
%         xlim([0 10]);
%         ylim([min max]);
%     end
% end


figure('color','w','WindowState','maximized')
for idx_iy = 1:nky
    for idx_ix = 1:nkx
        %         im_spectra = double(real(SPECTRA(:,idx_ix, idx_iy).*exp(-1i * (2* pi * (ppm_axis-4.7).'*45.7 * 0.00184))));% With first order phase correction
        im_spectra = double(real(SPECTRA(:,idx_ix, idx_iy)));
        hAxe = axes(...
            'Parent'    , gcf                        , ...
            'Box'       , 'on'                       , ...
            'LineWidth' , 2                          , ...
            'Position'  , [(idx_iy-1)*1/nky, 1-(idx_ix*1/nkx), 1/nky, 1/nkx], ...
            'XDir'      , 'reverse'                  , ...
            'XTick'     , []                         , ...
            'YTick'     , []                         , ...
            'XTickLabel', {''; ''}                   , ...
            'YTickLabel', {''; ''}                   , ...
            'TickDir'   , 'out'                      , ...
            'Ylim'      , [specmin specmax*1.15]     , ...
            'XLim'      , [0 10]                     ,...
            'XColor'    , [1 0 0]                    ,...
            'YColor'    , [1 0 0]                    ,...
            'Color'     , [1 1 1]);
        hLine = line( ...
            'Parent', hAxe            , ...
            'XData' , ppm_axis        , ...
            'YData' , im_spectra      , ...
            'Linewidth', 2         , ...
            'Color',[0 0 0]);
    end
end
