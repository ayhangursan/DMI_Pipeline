function Waterlocation=waterpeakloc(spectradata, xaxiszero)
%% output values will be used in HSVD
waterwindow=find(xaxiszero>3.5 & xaxiszero<5.9); % A wide window is picked as it will be used in baseline scan

[Maxvalue, Location]=max(real(spectradata(waterwindow,:,:,:)),[],1);

Waterlocation=xaxiszero(waterwindow(1)-1+squeeze(Location));
end