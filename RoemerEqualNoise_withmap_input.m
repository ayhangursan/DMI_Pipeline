function [RoemerEqualfid, Sensitivity_maps]=RoemerEqualNoise_withmap_input(fiddata,Sensitivity_maps,NCov,channelIndex)
% Roemer equal noise combination with sensitivity map is given as an input

if Sensitivity_maps==0
Sensitivity_maps=(squeeze(mean((fiddata(2:5,:,:,:,:)))));disp('Sensitivity map generated') % Mean value for second to fifth points of spectrum
else
    disp('Coil combination made with input reference sensitivity map')

end

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
end