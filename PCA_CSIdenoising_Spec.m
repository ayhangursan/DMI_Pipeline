function denoisedCSI_spec=PCA_CSIdenoising_Spec(spec_data,patch_size,Parameters)
% PCA based denoising of CSI data
disp('Starting PCA based denoising')
disp(strcat('Patch size:',num2str(patch_size)))

NP=size(spec_data,1);
extended_data=repmat(spec_data,[1 patch_size patch_size patch_size]);

denoisedCSI_spec=spec_data;

%% Determine eigenvalues for noise level
noisepatch=zeros(patch_size^3,NP*2);
n=patch_size^3;


for m=1:n
    [row,col,slice]=ind2sub([patch_size patch_size patch_size],m);
    noisepatch(m,1:NP)=real(extended_data(:,row+1+(Parameters.CSIdims(1)-ceil(patch_size/2)),col+1+(Parameters.CSIdims(2)-ceil(patch_size/2)),slice+1+(Parameters.CSIdims(3)-ceil(patch_size/2))));
    noisepatch(m,NP+1:end)=imag(extended_data(:,row+1+(Parameters.CSIdims(1)-ceil(patch_size/2)),col+1+(Parameters.CSIdims(2)-ceil(patch_size/2)),slice+1+(Parameters.CSIdims(3)-ceil(patch_size/2))));
end
[U,S,V]=svd(noisepatch,'econ');
noiseeigv = flip((diag(S)).^2);

%%

for mm=1:prod(Parameters.CSIdims)
    [Grid_row,Grid_col,Grid_slice]=ind2sub([Parameters.CSIdims(1) Parameters.CSIdims(2) Parameters.CSIdims(3)],mm);
    signalpatch=zeros(patch_size^3,NP*2);
    
    for m=1:n
        [row,col,slice]=ind2sub([patch_size patch_size patch_size],m);
        signalpatch(m,1:NP)=real(extended_data(:,Grid_row+row+(Parameters.CSIdims(1)-ceil(patch_size/2)),Grid_col+col+(Parameters.CSIdims(2)-ceil(patch_size/2)),Grid_slice+slice+(Parameters.CSIdims(3)-ceil(patch_size/2))));
        signalpatch(m,NP+1:end)=imag(extended_data(:,Grid_row+row+(Parameters.CSIdims(1)-ceil(patch_size/2)),Grid_col+col+(Parameters.CSIdims(2)-ceil(patch_size/2)),Grid_slice+slice+(Parameters.CSIdims(3)-ceil(patch_size/2))));
    end
    [U,S,V]=svd(signalpatch,'econ');
    eigv = flip((diag(S)).^2);
    
    nonnoiseeigenvalues=zeros(numel(eigv),1);
    nonnoiseeigenvalues(eigv-noiseeigv>0)=eigv(eigv-noiseeigv>0);
    C=diag(sqrt(flip(nonnoiseeigenvalues)));
    denoisedpatch=(U*C*V').';
    denoisedsignal=denoisedpatch(1:NP,ceil(n/2))+1i*denoisedpatch(NP+1:end,ceil(n/2));
    denoisedCSI_spec(:,Grid_row,Grid_col,Grid_slice)=denoisedsignal;
    clear denoisedsignal signalpatch U S V;
end

end
