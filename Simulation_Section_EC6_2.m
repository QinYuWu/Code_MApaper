%The quantile function of the supremum with Wasserstein uncertainty (Figure EC.2)

k=2;
epsilon=0.1;
alpha = 0.01:0.01:0.99;
n=length(alpha);
q_0=norminv(alpha,0,1);       %Quantile of standard normal distribution
q_1=qW_1(alpha,k,epsilon);    %Supremum quantile under MA1

for i=1:n                                                    
    q_2(i)=norminv(alpha(i),0,1)+(1-1/k)*(1-alpha(i))^(-1/k)*epsilon;   %Supremum quantile under MA2
end

plot(alpha,q_0,'k',alpha,q_1,'r',alpha,q_2,'b');
xlabel('s','Interpreter','latex')
ylabel('Quantile function','Interpreter','latex');
legend({'$\rm Standard~normal$','$F_{k,\epsilon|F_0}^1$','$F_{k,\epsilon|F_0}^2$'},'Interpreter','latex','Location','NorthWest' )
set(gcf, 'PaperPosition', [0 0 12 10]);
set(gcf, 'PaperSize', [12 10]);
saveas(gcf, 'quantile_W.pdf');






%The WR and MA approaches under Wasserstein uncertainty (Figure EC.3)

%%ES
k=2;
epsilon=0.1;
alpha = 0.01:0.01:0.99;
n=length(alpha);
for i=1:n
    ES_W1(i)=(1/(1-alpha(i)))*quadl(@(s) qW_1(s,k,epsilon),alpha(i),0.999);  %Robust ES under MA1
end


for i=1:n
    ES_W(i)=(1/(1-alpha(i)))*quadl(@(s) norminv(s,0,1),alpha(i),0.999)+epsilon/(1-alpha(i))^(1/k);  %Robust ES under WR and MA2
end

for i=1:n
    ES_norm(i)=(1/(1-alpha(i)))*quadl(@(s) norminv(s,0,1),alpha(i),0.999);   %Robust ES of standard normal distribution
end

plot(alpha,ES_norm,'k',alpha,ES_W,'-r',alpha,ES_W1,'b');
xlabel('$\alpha$','Interpreter','latex')
ylabel('Robust ${\rm ES}_\alpha$','Interpreter','latex');
legend({'$\rm Standard~normal$','${\rm ES}_\alpha^{\rm WR}$, ${\rm ES}_\alpha^{\rm MA_2}$','${\rm ES}_\alpha^{\rm MA_1}$'},'Interpreter','latex','Location','NorthWest')
set(gcf, 'PaperPosition', [0 0 12 10]);
set(gcf, 'PaperSize', [12 10]);
saveas(gcf, 'ES_W.pdf');

%%Power distortion 
k=2;
q=(1-1/k)^(-1);
epsilon=0.1;
alpha = linspace(1,10,100);
n=length(alpha);

for i=1:n
    PW_W(i)=quadl(@(s) norminv(s,0,1).*alpha(i).*s.^(alpha(i)-1),0.001,0.999)+epsilon*alpha(i)*(1/(q*(alpha(i)-1)+1))^(1/q); %Robust PD under WR
end

for i=1:n
    PW_W1(i)=quadl(@(s) qW_1(s,k,epsilon).*alpha(i).*s.^(alpha(i)-1),0.001,0.999);  %Robust PD under MA1
end


for i=1:n
   PW_W2(i)=quadl(@(s) norminv(s,0,1).*alpha(i).*s.^(alpha(i)-1),0.001,0.999)+alpha(i)*epsilon*(1-1/k)*beta(1-1/k,alpha(i)); %Robust PD under MA2
end

for i=1:n
    PW_norm(i)=quadl(@(s) norminv(s,0,1).*alpha(i).*s.^(alpha(i)-1),0.001,0.999);   %Robust PD of standard normal distribution
end

plot(alpha,PW_norm,'k',alpha,PW_W,'r',alpha,PW_W1,'b',alpha,PW_W2,'g');
axis([1 10 0 2.5]);
xlabel('k','Interpreter','latex')
ylabel('Robust ${\rm PD}_k$','Interpreter','latex');
legend({'$\rm Standard~normal$','${\rm PD}_k^{\rm WR}$','${\rm PD}_k^{\rm MA_1}$','${\rm PD}_k^{\rm MA_2}$'},'Interpreter','latex','Location','NorthWest')
set(gcf, 'PaperPosition', [0 0 12 10]);
set(gcf, 'PaperSize', [12 10]);
saveas(gcf, 'power_W.pdf');









%The WR and MA approaches under mean-variance uncertainty (Figure EC.4)

%%ES
mu = 0;
sigma = 1;

g=@(x) sqrt(x./(1-x));

alpha = linspace(0.001,0.999,1000);
ES0 = mu + sigma*sqrt(alpha./(1-alpha));                %Robust ES under WR and MA2
ES1 = zeros(1,1000);
count = 1;
for k=alpha
    ES1(count) = mu + sigma./(1-k)*integral(g,k,1);     %Robust ES under MA1
    count = count + 1;
end

plot(alpha,ES0,'-r',alpha,ES1,'b');
axis([0 1 0 40]);
xlabel('$\alpha$','Interpreter','latex')
ylabel('Robust ${\rm ES}_\alpha$','Interpreter','latex');
legend({'${\rm ES}_\alpha^{\rm WR}$, ${\rm ES}_\alpha^{\rm MA_2}$','${\rm ES}_\alpha^{\rm MA_1}$'},'Interpreter','latex','Location','NorthWest')
set(gcf, 'PaperPosition', [0 0 12 10]);
set(gcf, 'PaperSize', [12 10]);
saveas(gcf, 'ES.pdf');


%%Power distortion
mu = 0;
sigma = 1;

g=@(x,s) x.*(1-x).^(-0.5).*(1+x).^(s-1.5);

alpha = linspace(1,10,1000);
PD0 = mu + sigma*sqrt((alpha-1).^2./(2*alpha-1));             %Robust ES under WR
PD1 = mu + sigma*sqrt(pi)*gamma(alpha+1/2)./gamma(alpha);     %Robust ES under MA1
PD2 = zeros(1,1000);
count = 1;
for k=alpha
    PD2(count) = mu + sigma*k*0.5^k*integral(@(x)g(x,k),-1,1);  %Robust ES under MA2
    count = count + 1;
end

plot(alpha,PD0,'-r',alpha,PD1,'b',alpha,PD2,'g');
axis([1 10 0 6]);
xlabel('k','Interpreter','latex')
ylabel('Robust ${\rm PD}_k$','Interpreter','latex');
legend({'${\rm PD}_k^{\rm WR}$','${\rm PD}_k^{\rm MA_1}$','${\rm PD}_k^{\rm MA_2}$'},'Interpreter','latex','Location','NorthWest')
set(gcf, 'PaperPosition', [0 0 12 10]);
set(gcf, 'PaperSize', [12 10]);
saveas(gcf, 'power.pdf');

%%RVaR
mu = 0;
sigma = 1;
alpha = 0.5;

g=@(x) sqrt(x./(1-x));

bet = linspace(0.501,0.999,1000);
rvar0 = repmat(mu + sigma*sqrt(alpha./(1-alpha)),1,1000);        %Robust RVaR under WR
rvar1 = zeros(1,1000);
count = 1;
for k=bet
    rvar1(count) = mu + sigma./(k-alpha)*integral(g,alpha,k);    %Robust RVaR under MA1
    count = count + 1;
end

plot(bet,rvar0,'-r',bet,rvar1,'b');
axis([0.5 1 0.5 2.5]);
xlabel('$\beta$','Interpreter','latex')
ylabel('Robust ${\rm RVaR}_{\alpha,\beta}$','Interpreter','latex');
legend({'${\rm RVaR}_{0.5,\beta}^{\rm WR}$','${\rm RVaR}_{0.5,\beta}^{\rm MA_1}$'},'Interpreter','latex','Location','NorthWest','FontSize',11)
set(gcf, 'PaperPosition', [0 0 12 10]);
set(gcf, 'PaperSize', [12 10]);
saveas(gcf, 'RVaR.pdf');


%%Expectile
mu = 0;
sigma = 1;

g=@(x) sqrt(x./(1-x));

alpha = linspace(0.5,0.999,100);
c=alpha./(1-alpha);
b=(2.*alpha-1)./(1-alpha);

expectile0 = mu + sigma.*(c-1)./(2.*sqrt(c));     %Robust ex under WR and MA2

gg1 = @(s) atan((-s/(s - 1))^(1/2)) + (-s/(s - 1))^(1/2)/(s/(s - 1) - 1);
count = 1;
for k=b
    g1=@(x) x^3-(pi/2)*x^2+(k+1)*x-pi/2-k.*(1+x^2)*(pi/2-gg1(x^2/(1+x^2)));  %Robust ex under MA1
    a=double(solve(g1));
    expectile1(count) = mu + sigma*(a(a>0));
    count = count + 1;
end

plot(alpha,expectile0,'r',alpha,expectile1,'b');
axis([0.5 1 0 30]);
xlabel('$\alpha$','Interpreter','latex')
ylabel('Robust ${\rm ex}_\alpha$','Interpreter','latex');
legend({'${\rm ex}_\alpha^{\rm WR}$, ${\rm ex}_\alpha^{\rm MA_2}$','${\rm ex}_\alpha^{\rm MA_1}$'},'Interpreter','latex','Location','NorthWest')
set(gcf, 'PaperPosition', [0 0 12 10]);
set(gcf, 'PaperSize', [12 10]);
saveas(gcf, 'expectile.pdf');




