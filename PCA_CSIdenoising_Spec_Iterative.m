function denoisedCSI_spec=PCA_CSIdenoising_Spec_Iterative(spec_data,patch_size,Parameters)
%% PCA based denoising of CSI data
disp('Starting PCA based denoising.')
disp(strcat('Patch size:',num2str(patch_size)))

NP=size(spec_data,1);
extended_data=repmat(spec_data,[1 3 3 3]);

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

%% Denoise by elimination eigenvalues equal/below noise eigenvalues
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
    %   This part works but not iterative
        denoisedsignal=denoisedpatch(1:NP,ceil(n/2))+1i*denoisedpatch(NP+1:end,ceil(n/2));
        denoisedCSI_spec(:,Grid_row,Grid_col,Grid_slice)=denoisedsignal;
    
    denoisedpatch=reshape(denoisedpatch(1:NP,:),[NP patch_size patch_size patch_size])+1i*reshape(denoisedpatch(NP+1:end,:),[NP patch_size patch_size patch_size]);

    for m=1:n
        [row,col,slice]=ind2sub([patch_size patch_size patch_size],m);
        AP_ind=Grid_row+row+(Parameters.CSIdims(1)-ceil(patch_size/2));
        RL_ind=Grid_col+col+(Parameters.CSIdims(2)-ceil(patch_size/2));
        FH_ind=Grid_slice+slice+(Parameters.CSIdims(3)-ceil(patch_size/2));
        extended_data(:,AP_ind,RL_ind,FH_ind)=denoisedpatch(:,row,col,slice);
    end
    clear denoisedsignal signalpatch U S V;
    %%
end
cutoff_AP=[Parameters.CSIdims(1)+1: 2*Parameters.CSIdims(1)];
cutoff_RL=[Parameters.CSIdims(2)+1: 2*Parameters.CSIdims(2)];
cutoff_FH=[Parameters.CSIdims(3)+1: 2*Parameters.CSIdims(3)];
% denoisedCSI_spec=extended_data(:,cutoff_AP,cutoff_RL,cutoff_FH);
disp('Finished PCA based denoising.')

end
