clear,clc;
addpath('tSVD','proxFunctions','solvers','twist','Datasets');
addpath('ClusteringMeasure', 'LRR', 'Nuclear_norm_l21_Algorithm', 'unlocbox');

load('my_NGs.mat');
fprintf('Clustering on NGs \n');
cls_num = length(unique(gt));
%%
X{1} = X1; X{2} = X2; X{3} = X3;
for v=1:3
    [X{v}]=NormalizeData(X{v});
end

K = length(X); N = size(X{1},2);
lambda = 0.01; 
dim_H = 150;

for k=1:K
    Z{k} = zeros(N,N);
    E1{k} = zeros(size(X{k},1),N);
    E2{k} = zeros(size(X{k},1),N);
    H = zeros(dim_H,N);
    P{k} = zeros(size(X{k},1),dim_H);
    Y1{k} = zeros(size(X{k},1),N);
    Y2{k} = zeros(size(X{k},1),N);
    C{k} = zeros(N,N);
    W{k} = zeros(N,N);
end
Z{K+1} = zeros(N,N);
E3 = zeros(dim_H,N);
Y3 = zeros(dim_H,N);
C{K+1} = zeros(N,N);
W{K+1} = zeros(N,N);

w = zeros(N*N*K,1);
c = zeros(N*N*K,1);
dim1 = N;dim2 = N;dim3 = K;
sX = [N, N, K+1];

parOP         =    false;
ABSTOL        =    1e-6;
RELTOL        =    1e-4;

Isconverg = 0;epson = 1e-7;
iter = 0;
mu = 10e-5; max_mu = 10e10; pho_mu = 2;
rho = 10e-5; max_rho = 10e12; pho_rho = 2;
tic;

while(Isconverg == 0)
    fprintf('----processing iter %d--------\n', iter+1);
    for k=1:K
        
        % update Zv{k}
        tmp = X{k}'*Y1{k} + mu*X{k}'*X{k} - mu*X{k}'*E1{k} - W{k} +  rho*C{k};
        Z{k} = (rho*eye(N,N)+ mu*X{k}'*X{k})\tmp;
        
        % update P{k}
        U_a = X{k}-E2{k} + Y2{k}/mu;
        U_b = H*U_a';
        [svd_U,~,svd_V] = svd(U_b,'econ');
        P{k} = svd_V*svd_U';
        
        %3 update Y1 and Y2
        Y1{k} = Y1{k} + mu*(X{k}-X{k}*Z{k}-E1{k});
        Y2{k} = Y2{k} + mu*(X{k}-P{k}*H - E2{k});
    end
    
    % update E{k}
    F = [];
    F = [X{1}-X{1}*Z{1}+Y1{1}/mu;X{2}-X{2}*Z{2}+Y1{2}/mu;X{3}-X{3}*Z{3}+Y1{3}/mu];
    F = [F;X{1}-P{1}*H+Y2{1}/mu;X{2}-P{2}*H+Y2{2}/mu;X{3}-P{3}*H+Y2{3}/mu];
    F = [F;H-H*Z{K+1}+Y3/mu];
    [Econcat] = solve_l1l2(F,lambda/mu);

    E1{1} = Econcat(1:size(X{1},1),:);
    E1{2} = Econcat(size(X{1},1)+1:size(X{1},1)+size(X{2},1),:);
    E1{3} = Econcat(size(X{1},1)+size(X{2},1)+1:size(X{1},1)+size(X{2},1)+size(X{3},1),:);
    tmp_index_E1 = size(X{1},1)+size(X{2},1)+size(X{3},1);
    E2{1} = Econcat(tmp_index_E1+1:tmp_index_E1+size(X{1},1),:);
    E2{2} = Econcat(tmp_index_E1+size(X{1},1)+1:tmp_index_E1+size(X{1},1)+size(X{2},1),:);
    E2{3} = Econcat(tmp_index_E1+size(X{1},1)+size(X{2},1)+1:tmp_index_E1+size(X{1},1)+size(X{2},1)+size(X{3},1),:);
    tmp_index_E2 = tmp_index_E1 + size(X{1},1)+size(X{2},1)+size(X{3},1);
    E3 = Econcat(tmp_index_E2+1:end,:);
    
    % update Zs
    tmp1 = H'*Y3 + mu*(H'*H-H'*E3) - W{K+1} + rho*C{K+1};
    Z{K+1} = (mu*(H'*H) + rho*eye(N,N))\tmp1;
    
    % update H
    PTP = zeros(dim_H,dim_H);
    H_c = zeros(dim_H,N);
    for inner_iter_K = 1:K
        PTP = PTP + P{inner_iter_K}'*P{inner_iter_K};
        H_c = H_c + P{inner_iter_K}'*Y2{inner_iter_K}+mu*P{inner_iter_K}'*X{inner_iter_K}-mu*P{inner_iter_K}'*E2{inner_iter_K};
    end
    H_left = mu*(PTP + eye(dim_H,dim_H));
    H_right = mu*(Z{K+1}*Z{K+1}'-Z{K+1}-Z{K+1}');
    H_all = H_c + Y3*Z{K+1}'-Y3 + mu*E3 - mu*E3*Z{K+1}';
    H = lyap(H_left,H_right,H_all);
    
    % update Y3
    Y3 = Y3 + mu*(H-H*Z{K+1}-E3);
    
    % update C
    Z_tensor = cat(3, Z{:,:});
    W_tensor = cat(3, W{:,:});
    z = Z_tensor(:);
    w = W_tensor(:);
    
    [c, ~] = wshrinkObj(z + 1/rho*w,1/rho,sX,0,3);
    C_tensor = reshape(c, sX);
    for inter_C = 1:K+1
        C{inter_C} = C_tensor(:,:,inter_C);
    end
    
    % update W
    w = w + rho*(z - c);
    W_tensor = reshape(w, sX);
    for inter_W = 1:K+1
        W{inter_W} = W_tensor(:,:,inter_W);
    end
    
    %% coverge condition
    Isconverg = 1;
    for k=1:K
        if (norm(X{k}-X{k}*Z{k}-E1{k},inf)>epson)
            history.norm_XXZE1 = norm(X{k}-X{k}*Z{k}-E1{k},inf);
            fprintf('norm_XXZE1 %7.10f\n', history.norm_XXZE1);
            Isconverg = 0;
        end
        
        if (norm(Z{k}-C{k},inf)>epson)
            history.norm_ZC = norm(Z{k}-C{k},inf);
            fprintf('norm_ZC %7.10f\n', history.norm_ZC);
            Isconverg = 0;
        end
    end
    
    if (norm(Z{K+1}-C{K+1},inf)>epson)
        history.norm_ZCs = norm(Z{K+1}-C{K+1},inf);
        fprintf('norm_ZCs %7.10f\n', history.norm_ZCs);
        Isconverg = 0;
    end

    if (iter>200)
        Isconverg  = 1;
    end
    
    iter = iter + 1;
    mu = min(mu*pho_mu, max_mu);
    rho = min(rho*pho_rho, max_rho);
end

S = 0;
for k=1:K+1
    S = S + abs(Z{k})+abs(Z{k}');
end

NMI_total = zeros(1,30);
ACC_total = zeros(1,30);
F_total = zeros(1,30);
AVG_total = zeros(1,30);
Precision_total = zeros(1,30);
RI_total = zeros(1,30);
AR_total = zeros(1,30);
Recall_total = zeros(1,30);
for i_ave = 1:30
    C = SpectralClustering(S,cls_num);
    [A, nmi, avgent] = compute_nmi(gt,C);
    
    ACC = Accuracy(C,double(gt));
    [f,p,r] = compute_f(gt,C);
    [AR,RI,MI,HI]=RandIndex(gt,C);
    toc;
    fprintf('NMI: %f, ACC: %f, F-Score: %f\n',nmi,ACC,f);
    NMI_total(i_ave) = nmi;
    ACC_total(i_ave) = ACC;
    F_total(i_ave) = f;
    AVG_total(i_ave) = avgent;
    Precision_total(i_ave) = p;
    RI_total(i_ave) = RI;
    AR_total(i_ave) = AR;
    Recall_total(i_ave) = r;
end
NMI_mean = mean(NMI_total); NMI_std = std(NMI_total);
ACC_mean = mean(ACC_total); ACC_std = std(ACC_total);
F_mean = mean(F_total); F_std = std(F_total);
AVG_mean = mean(AVG_total); AVG_std = std(AVG_total);
Precision_mean = mean(Precision_total); Precision_std = std(Precision_total);
RI_mean = mean(RI_total); RI_std = std(RI_total);
AR_mean = mean(AR_total); AR_std = std(AR_total);
Recall_mean = mean(Recall_total); Recall_std = std(Recall_total);

fprintf('NMI: %f(%f), ACC: %f(%f), F-Score: %f(%f), Precision: %f(%f), AVG: %f(%f), RI: %f(%f), AR: %f(%f), Recall: %f(%f)\n',...
    NMI_mean,NMI_std,ACC_mean,ACC_std,F_mean,F_std,Precision_mean,Precision_std,AVG_mean,AVG_std,RI_mean,RI_std,AR_mean,AR_std,Recall_mean,Recall_std);
