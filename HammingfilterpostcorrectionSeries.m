function [idealhammingwindow,Correctionmask,filteredCSI]=HammingfilterpostcorrectionSeries(fiddata,acqpattern,Parameters)

filteredCSI=zeros(size(fiddata));
grid=Parameters.CSIdims;

for k=1:prod(grid)
    
    [x,y,z]= ind2sub(grid,k);
    filterfuncx=0.54+0.46*cos(2*pi*x/grid(1));
    filterfuncy=0.54+0.46*cos(2*pi*y/grid(2));
    filterfuncz=0.54+0.46*cos(2*pi*z/grid(3));
    hammingwindow(k)=filterfuncx*filterfuncy*filterfuncz;
end
hammingwindow=fftshift(reshape(hammingwindow,grid));
idealhammingwindow=circshift(hammingwindow,[1 1 0]);

Correctionmask=idealhammingwindow./acqpattern;

dimmatchedmask=repmat(Correctionmask,[1 1 1 Parameters.dims(1)]);
dimmatchedmask=permute(dimmatchedmask,[4 1 2 3]);
for m=1:Parameters.dims(Parameters.Index.dynIndex)
    if  isempty(Parameters.Index.channelIndex)==1
        
        filteredCSI(:,:,:,:,m)=squeeze(fiddata(:,:,:,:,m)).*dimmatchedmask;
        
    else
        for n=1:Parameters.dims(Parameters.Index.channelIndex)
            
            filteredCSI(:,n,:,:,:,m)=squeeze(fiddata(:,n,:,:,:,m)).*dimmatchedmask;
            
        end
        
    end
end