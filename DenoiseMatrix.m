function OutMatrix=DenoiseMatrix(InpMatrix)
%% Debug
% InpMatrix=noisepatch;
% InpMatrix=signalpatch;
%
n=size(InpMatrix,1);
m=size(InpMatrix,2);
r=min(m,n);
[U,S,V] = svd(InpMatrix,'econ');
eigv = flip((diag(S)).^2);
lam_r = eigv(1) / n;
clam = 0;
sigma2 = NaN;

for p=1:r
    lam = eigv(p) / n;
    clam = clam+lam;
    gam = (m-r+p)/n;
    sigsq1 = clam/(p) / max(gam,1);
    sigsq2 = (lam-lam_r)/4/sqrt(gam);
    if(sigsq2 < sigsq1)
        sigma2 = sqrt(sigsq1);
        cutoff_p = p+1;
    end
end
cutoff_p = r-cutoff_p;
eigv = flip(eigv);

if(cutoff_p > 1)
    Snew = zeros(size(S));
    Sd = diag(sqrt(eigv(1:cutoff_p)));
    Snew(1:cutoff_p,1:cutoff_p) = Sd;
end

rebuilt_data = reshape(U*Snew*V',size(InpMatrix));
% noise=sigma2
% 
% figure
% subplot(1,3,1)
% histogram(S,50)
% subplot(1,3,2)
% histogram(eigv,50)
% subplot(1,3,3)
% histogram(Snew,50)
% 
% figure
% subplot(2,2,1)
% plot(V9349B.xaxis,real(complex(InpMatrix(62,1:256),InpMatrix(62,256+1:end))).')
% set(gca,'XDir','reverse')
% xlim([-20 10])
% 
% subplot(2,2,2)
% plot(V9349B.xaxis,real(complex(rebuilt_data(62,1:256),rebuilt_data(62,256+1:end))).')
% set(gca,'XDir','reverse')
% xlim([-20 10])
% 
% subplot(2,2,3:4)
% raw=real(complex(InpMatrix(62,1:256),InpMatrix(62,256+1:end))).';
% raw=raw./std(raw(V9349B.noisewindow));
% DN=real(complex(rebuilt_data(62,1:256),rebuilt_data(62,256+1:end))).';
% DN=DN./std(DN(V9349B.noisewindow));
% 
% plot(V9349B.xaxis,raw)
% hold on
% plot(V9349B.xaxis,DN)
% hold off
% set(gca,'XDir','reverse')
% xlim([-20 10])
% %%
OutMatrix=rebuilt_data;
%%