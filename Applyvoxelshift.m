function outputkspace=Applyvoxelshift(inputkspace,GridShift,Parameters)
%% Apply voxelshift
% Use before channel combinationcombined
kx1 = (Parameters.CSIdims(1)-1)/Parameters.CSIdims(1);
ky1 = (Parameters.CSIdims(2)-1)/Parameters.CSIdims(2);
kz1 = (Parameters.CSIdims(3)-1)/Parameters.CSIdims(3);
[kx,ky,kz] = ndgrid(linspace(-kx1,kx1,Parameters.CSIdims(1)),linspace(-ky1,ky1,Parameters.CSIdims(2)),linspace(-kz1,kz1,Parameters.CSIdims(3)));
FIDshifted=zeros(size(inputkspace));

if isempty(Parameters.Index.channelIndex)==1
    for npfid = 1:Parameters.NP
        FIDshifted(npfid,:,:,:) = squeeze(inputkspace(npfid,:,:,:)).*exp(-pi*1i*kx*GridShift(1) - pi*1i*ky*GridShift(2) - pi*1i*kz*GridShift(3));
    end
else
    for npfid = 1:Parameters.NP
        for nchan=1:size(inputkspace,Parameters.Index.channelIndex)
            FIDshifted(npfid,nchan,:,:,:) = squeeze(inputkspace(npfid,nchan,:,:,:)).*exp(-pi*1i*kx*GridShift(1) - pi*1i*ky*GridShift(2) - pi*1i*kz*GridShift(3));
        end
    end
end
outputkspace=FIDshifted;