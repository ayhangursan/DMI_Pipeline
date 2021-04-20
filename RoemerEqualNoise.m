function [RoemerEqualfid, Sensitivity_maps]=RoemerEqualNoise(fiddata,NCov,channelIndex)

Sensitivity_maps=(squeeze(mean((fiddata(2:5,:,:,:,:))))); % Mean value for second to fifth points of spectrum
dim=size(fiddata);
grid_dims=dim(channelIndex+1:end);
numberofloc=prod(size(Sensitivity_maps))/dim(channelIndex);
RoemerEqualfid=zeros([dim(1) grid_dims]);
for k=1:numberofloc
    S=Sensitivity_maps(:,k);
    U=pinv(sqrt(S'*pinv(NCov)*S))*S'*pinv(NCov);
    V=U*squeeze(fiddata(:,:,k)).';
    RoemerEqualfid(:,k)=V;
end
Sensitivity_maps=reshape(Sensitivity_maps,[dim(2:end)]);
end