function [idealhammingwindow,Correctionmask,filteredCSI]=Hammingfilterpostcorrection(fiddata,acqpattern,Parameters)

filteredCSI=zeros(size(fiddata));
grid=Parameters.CSIdims;
dims=size(fiddata);
for k=1:prod(grid)
    
    [x,y,z]= ind2sub(grid,k);
    filterfuncx=0.54+0.46*cos(2*pi*x/grid(1));
    filterfuncy=0.54+0.46*cos(2*pi*y/grid(2));
    filterfuncz=0.54+0.46*cos(2*pi*z/grid(3));
    hammingwindow(k)=filterfuncx*filterfuncy*filterfuncz;
end

hammingwindow=fftshift(reshape(hammingwindow,grid));

idealhammingwindow=circshift(hammingwindow,1,1);
idealhammingwindow=circshift(idealhammingwindow,1,2);
idealhammingwindow=circshift(idealhammingwindow,1,3);

Correctionmask=idealhammingwindow./acqpattern;
% In datasets acquired with old k-space weighting patch shutter leaves some
% points with no line acquired(acqpattern=0). This dividing by zero may lead to inf points in correction mask.
% To avoid that check if acqpattern has zeros, if so correctionmask(correctionmask==inf)=0
if numel(find(~acqpattern))>0
    disp('Non-acquired k-space points in the dataset.')
    Correctionmask(isinf(Correctionmask))=0;
end

dimmatchedmask=repmat(Correctionmask,[1 1 1 dims(1)]);
dimmatchedmask=permute(dimmatchedmask,[4 1 2 3]);
if  isempty(Parameters.Index.channelIndex)==1
    
    filteredCSI(:,:,:,:)=squeeze(fiddata(:,:,:,:)).*dimmatchedmask;
    
else
    for n=1:Parameters.dims(Parameters.Index.channelIndex)
        
        filteredCSI(:,n,:,:,:)=squeeze(fiddata(:,n,:,:,:)).*dimmatchedmask;

    end
    
end
disp('Hamming filter correction at post-processing applied.')