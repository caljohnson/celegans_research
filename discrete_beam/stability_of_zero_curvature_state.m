%Stability of Zero Curvature State

%linear stability analysis for the discretized beam PDE
%about K = 0

%size n
n = 10;
e = ones(n,1);

%second derivative operator
A = spdiags([e -2*e e], [0 1 2], n-2, n);

%zero curvature
K = zeros(n-2,1);
    
%Discretized PDE: y_t = A^T A y
%                A y_t = A (A^T A y)
%                   K = A A^T K

%stability of zero state
%det. by eigvals of A A^T
[V, D] = eig(full(-A*A'));
lambdas = diag(D);
lambdas([end, end-1, 1])
plot(V(:,end)); hold on;
plot(V(:,end-1)); 
plot(V(:,1));



