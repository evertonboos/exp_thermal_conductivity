function S = ort_coeff_vector(t,mydata)
% Assembles coefficient vectors related to the right hand side related to
% time 't' (scalar).

%% Parameters
[n,N,t0,tt0,et,alpha,beta,D,x,y,P,hat_D,d0,dn,d0_e1,dn_en,X,Y] = mydata.cheby;
run('data.m')
run('data_adequation.m')

%% Computing vectors
u = zeros((n+1)^2,1); v = zeros((n+1)^2,1);
for j = 1:n+1
    int = (j-1)*(n+1)+1:j*(n+1);
    b = l1*(-d0*h1(y(j))*f1(y(j),t) + dn*h2(y(j))*f2(y(j),t));
    d = l2*(-d0*h3(x(j))*f3(x(j),t) + dn*h4(x(j))*f4(x(j),t));
    %d = l2*(-d0*h3(x(end-j+1))*f3(x(end-j+1),t) + dn*h4(x(end-j+1))*f4(x(end-j+1),t));

    u(int) = [b(1); zeros(n-1,1); b(end)];
    v(int) = [d(1); zeros(n-1,1); d(end)];
end

S = (1/l1^2)*u + (1/l2^2)*P*v + g(X,Y,t); % S(t)

end

