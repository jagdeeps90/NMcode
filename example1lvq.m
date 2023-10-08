%% problem statement
% u(x) = sin(x) and \bar{x} = 1, % u'(x) = cos(1) = 0.5403023
%show the error for different values of h for first order and second order
%FDM formula on loglog scale to better visualize the erros

%code 

exact = cos(1); % we need this to evaluate the erros Du(x) -u'(x) here u'(x) = cos(1)

h_val = [1e-1 5e-2 1e-2 5e-3 1e-3];
disp(' ')
disp('       h             exact            Dpu            Dmu             D0u             D3u')
for i= 1:length(h_val)
    h = h_val(i);
    %forward
    Dpu = (sin(1+h)-sin(1))/h;
    Dmu = (sin(1) - sin(1-h))/h;
    D0u = (sin(1+h) - sin(1-h))/(2*h);
    D3u = (2*sin(1+h)+3*sin(1)-6*sin(1-h)+sin(1-2*h))/(6*h);

    %errors
       % errors:
   Epu(i) = Dpu - exact;
   Emu(i) = Dmu - exact;
   E0u(i) = D0u - exact;
   E3u(i) = D3u - exact;
   % print line of table:
   disp(sprintf('%13.4e  %13.4e  %13.4e   %13.4e   %13.4e  %13.4e',...
                 h, exact, Epu(i),Emu(i),E0u(i),E3u(i)))
 end

% plot errors:
figure(1)
loglog(h_val,abs(Epu),'-o','LineWidth',1);
axis([5e-4 .2 1e-12 1])
hold on
loglog(h_val,abs(Emu),'--s','LineWidth',2)
loglog(h_val,abs(E0u),'+-','LineWidth',2)
loglog(h_val,abs(E3u),'*-','LineWidth',2)
hold off

legend('Dpu','Dmu','D0u','D3u','Interpreter','Latex','Location','southeast')

xlabel('h','Interpreter','latex');
ylabel('E(h)','Interpreter','latex');

set(gca,'fontsize',18);
