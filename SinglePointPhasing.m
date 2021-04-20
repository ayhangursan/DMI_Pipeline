function PhasedSpec=SinglePointPhasing(spec)

NumPoints=size(spec,1);
Numlines=prod(size(spec))/NumPoints;

%Vectorize spectra
specvec=reshape(spec,[NumPoints Numlines]);

%Find maximum value
[Magmax Magmaxindex]=max(abs(specvec),[],1);

% Find phase at maximum value
phase=zeros(1,Numlines);
for m=1:Numlines
phase(m)=atan2(imag(specvec(Magmaxindex(m),m)),real(specvec(Magmaxindex(m),m)));
end
% Rephase signal so that real(max)=mag(max)
PhasedSpec=specvec.*(exp(-1i*phase));

PhasedSpec=reshape(PhasedSpec,size(spec));
