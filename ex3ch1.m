%domain definition [ax,bx] E [0,1] imposing bc
% u'(0)=sigma, u(1)=beta: let us consider sigma=0, beta=3
%Mehtod-1 and Method-2 to handle the neuman BC
%Following the Randell Leveque Book on Numerical Methods
ax = 0;
bx = 1;
sigma = 0;
beta = 3;
% RHS function e^x
f = @(x) exp(x);  % right hand side function
utrue = @(x) exp(x) + (sigma-exp(ax))*(x - bx) + beta - exp(bx);  % true soln

%grid points
m = [10 20 40 80];
mesh = length(m);
E = zeros(mesh,1);
hval = zeros(mesh,1);
Ux = zeros(11,2);
for method = 1:2
    for i = 1:mesh
        n = m(i);
        x = linspace(ax,bx,n+1)';
        h = x(2)-x(1);
        hval(i) = h;
        % set up matrix A (using sparse matrix storage):
        e = ones(n+1,1);
        A = 1/h^2 * spdiags([e -2*e e], [-1 0 1], n+1, n+1);
        % fix up values in first and last row for BCs:
        A(1,1) = -1/h;
        A(1,2) = 1/h;
        A(n+1,n) = 0;
        A(n+1,n+1) = 1;
        
        % Right hand side:
        F = f(x);
        %F(1) = sigma;  % Approach 1 first order
        F(1) = sigma+ h/2 * F(1)*(method-1); %Approach 2 2nd order
        F(n+1) = beta;
        
        % solve linear system:
        U = A\F;
        
        
        uhat = utrue(x);
        err = abs(U - uhat);
        E(i,method) = max(abs(err));
       if i==1
           x1=x;
           Ux(:,method)=U;
       end
    end
end

% plot solutions:
figure(1)
clf
xfine = linspace(ax,bx,101);
ufine = utrue(xfine);
plot(xfine,ufine,'-k','LineWidth',1);
hold on
plot(x1,Ux(:,1),'+b','MarkerSize',12);
plot(x1,Ux(:,2),'or','MarkerSize',12);
hold off
axis([-.1 1.1 2.2 3.2])
set(gca,'fontsize',18)
xlabel('x');
ylabel('U(x)');
legend('Exact','Method-1','Method-2','Location','northwest');

% log-log plot of errors:
figure(2)
clf
loglog(hval,E(:,1),'b+-','MarkerSize',12);
hold on
loglog(hval,E(:,2),'ro-','MarkerSize',12);
hold off
axis([5e-3 .2 1e-5 1])
set(gca,'fontsize',18)
legend('Method-1','Method-2');
xlabel('h');
ylabel('E(h)');
