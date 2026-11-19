function [PCA_mode, Time_coef, k, cumulative_variance]=Extract_flow_PCA(U,V,W)
%function [D,V,B, nPCA,CumsumVariance]=Extract_flow_PCA(Xdata,Ydata,Zdata)

if 1 % the following code comes from ChatGPT.
% INPUT (U,V,W): velocity component matrices (np spatial points × nt time points)

% Stack components: each row = spatial mode, each column = time
X = [U; V; W];  % size: 3np × nt

% Step 2: Subtract temporal mean (mean over columns)
X_mean = mean(X, 2);
X_centered = X - X_mean;

% Step 3: PCA via SVD
[U_pca, S, V_pca] = svd(X_centered, 'econ');  % X_centered = U*S*V'

% Step 4: Determine how many modes to retain (99.99% cumulative variance)
singular_values = diag(S);
variance_explained = singular_values.^2 / sum(singular_values.^2);
cumulative_variance = cumsum(variance_explained);
k = find(cumulative_variance >= 0.9999, 1, 'first');  % number of modes to retain

% Step 5: Extract spatial modes and temporal coefficients
Phi = U_pca(:, 1:k);            % spatial modes (3np × k)
A = S(1:k,1:k) * V_pca(:,1:k)'; % temporal coefficients (k × nt)
% OUTPUT:
% Phi -- each column is a spatial PCA mode (length 3n)
% A -- each row is the time series (coefficients) for the corresponding PCA mode
% X_approx = Phi * A + X_mean;
PCA_mode = [X_mean Phi];
Time_coef = A';
else


%% Pre-process the input data

[np nt] = size(Zdata);
% new array with space 3D->1D column
data = zeros(3*np,nt);

    for i_nt=1:nt
      % reshape matrix for a column
      data(:,i_nt) = reshape([Xdata(:,i_nt),Ydata(:,i_nt),Zdata(:,i_nt)],[3*np,1]); 
    end 

 

% linear principal component analysis ====================================
psi = data';
nt=size(psi,1); % number of experimental observations (e.g. frames recorded) TIME
np3=size(psi,2); % number of features recorded with each observation SPACE


% mean of all n observations (represented as (n x m)-matrix)
psi0=repmat(mean(psi),nt,1); 

% mean-corrected (n x m)-observation matrix 
Delta=psi-psi0; 
%Delta=Delta./repmat(std(Delta),nt,1);


%% PCA method to calculate principle components

C=Delta'*Delta;
[V,D]=eig(C);
[D,ind]=sort(diag(D),'descend'); % sort eigenvalues (and corresponding eigenvectors) by decreasing magnitude
%V=V(:,ind)./norm(V(:,ind)); % normalization %==>这一步可能有错！

CumsumVariance = cumsum(D)./sum(D);
nPCA =find(CumsumVariance>=0.9999,1); 
CumsumVariance=CumsumVariance(1:nPCA);

V = [psi0(1,:)'  V(:,1:nPCA)];
B = nan(nt,nPCA);
for i_PCA = 2:nPCA %这里改正了！原来的从1开始数，错了！
    B_temp = V(:,i_PCA)\Delta';
    B(:,i_PCA) = B_temp;
end
    
%sum(D(1:8))/sum(D)

end