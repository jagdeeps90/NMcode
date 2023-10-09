%% simple fdm formula
ax = 0;
bx = 1;
alpha = 0;
beta = 0;
%RHS function
f = @(x_int) sin(x_int);
%true solution
utrue = @(x) -sin(x) + sin(1)*x;
%number of grid points

n = 32;
%grid spacing
 
x = linspace(ax,bx,n)';
h = x(2)-x(1);
x_int = linspace(ax+h,bx-h,n-2)';

%%
% set up matrix A (using sparse matrix storage):
U = zeros(n,1);
U(1)= alpha;
U(n) = beta;

A = zeros(n-2,n-2);
A(1,1) = -2/h^2;
for i =2:n-2
    A(i, i-1) = 1/h^2;
    A(i,i)    = -2/h^2;
end
for i = 1:n-3
    A(i,i+1) = 1/h^2;
end

%%
% RHS 
rhs = f(x_int);
rhs(1) = rhs(1) - alpha/h^2;
rhs(n-2) = rhs(n-2) - beta/h^2;

%%
%solve linear system
Uc = A\rhs;
%include the bc
for i = 2:n-1
    U(i) = Uc(i-1);
end
%%
%U = [alpha;U;beta];
uhat = utrue(x);
figure(1)
plot(x,U,'*k','LineWidth',2);
hold on;
plot(x,uhat,'--o')
xlabel('x');
ylabel('U');
set(gca,'fontsize',18)
legend('FDM','Exact','Location','southeast');
hold off;

figure(2)
err = abs(U) - abs(uhat);
plot(x,err);
xlabel('x');
ylabel('E(h)');
set(gca,'fontsize',18)
