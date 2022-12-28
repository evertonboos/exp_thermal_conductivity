function TK = isox_generate_T(KK,T0,mydata)
% calculate temperature T at time step N given a vector K = K(:)
% containing the heat conductivity parameter, where k = k(x).
% Input is only k(x), i.e, a vector with n+1 entries.

%% Parameters
[n,N,t0,tt0,et,alpha,beta,D,x,y,P,hat_D,d0,dn,d0_e1,dn_en,X,Y] = mydata.cheby;
run('data.m')
run('data_adequation.m')

% create a conductivity matrix from the input vector
KK = KK(:);
K = KK*zeros(1,n+1);
for i = 1:n+1
    K(:,i) = KK(:);
end

%% Assembling coefficient matrix
D1 = D(:,2:end-1);
D2 = D(2:end-1,:);

% Computing matrix
F = sparse((n+1)^2,(n+1)^2);
G = F;

for j = 1:n+1
    int = (j-1)*(n+1)+1:j*(n+1);
    A = l1*h1(y(j))*d0_e1 + D1*diag(K(2:end-1,j))*D2 - l1*h2(y(j))*dn_en;
    B = l2*h3(x(j))*d0_e1 + D1*diag(K(j,2:end-1))*D2 - l2*h4(x(j))*dn_en;
    F(int,int) = [A(1,:); D2*diag(K(:,j))*D; A(end,:)];
    G(int,int) = [B(1,:); D2*diag(K(j,:))*D; B(end,:)];
end

M = (1/l1^2)*F + (1/l2^2)*P*G*P';

%% IVP
TK = zeros((n+1)^2,length(tt0));
TK(:,1) = T0;

C = c(X,Y);
%C = flip(C);
C = C(:);
C = sparse(1:(n+1)^2, 1:(n+1)^2, C); % 'C' matrix, diagonal sparse

for i = 1:length(tt0)-1 % Crank-Nicolson's method
    
    % t_{i+1}
    h = tt0(i+1)-tt0(i); % \Delta t
    
    Mm = C - (h/2)*M;
    Mp = C + (h/2)*M;
    %Mm = speye((n+1)^2) - (h/2)*M;
    %Mp = speye((n+1)^2) + (h/2)*M;
    
    % right-hand side (vector S(t))
    ua = zeros((n+1)^2,1); va = zeros((n+1)^2,1);
    up = zeros((n+1)^2,1); vp = zeros((n+1)^2,1);
    for j = 1:n+1
        int = (j-1)*(n+1)+1:j*(n+1);
        ba = l1*(-d0*h1(y(j))*f1(y(j),tt0(i)) + dn*h2(y(j))*f2(y(j),tt0(i)));
        da = l2*(-d0*h3(x(j))*f3(x(j),tt0(i)) + dn*h4(x(j))*f4(x(j),tt0(i)));
        bp = l1*(-d0*h1(y(j))*f1(y(j),tt0(i+1)) + dn*h2(y(j))*f2(y(j),tt0(i+1)));
        dp = l2*(-d0*h3(x(j))*f3(x(j),tt0(i+1)) + dn*h4(x(j))*f4(x(j),tt0(i+1)));

        ua(int) = [ba(1); zeros(n-1,1); ba(end)];
        va(int) = [da(1); zeros(n-1,1); da(end)];
        up(int) = [bp(1); zeros(n-1,1); bp(end)];
        vp(int) = [dp(1); zeros(n-1,1); dp(end)];
    end

    Sa = (1/l1^2)*ua + (1/l2^2)*P*va + g(X,Y,tt0(i)); % S(t_i) (anterior)
    Sp = (1/l1^2)*up + (1/l2^2)*P*vp + g(X,Y,tt0(i+1)); % S(t_{i+1}) (posterior)
    
    % Crank-Nicolson's method
    TK(:,i+1) = Mm\( Mp*TK(:,i) + (h/2)*(Sa+Sp) );
end

%TK = TK(2:end);
% select time steps
%TK = TK(:,(1:N)*(et+1));

end

