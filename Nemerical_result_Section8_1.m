%Load data

xReturn = readtable('return.csv'); 
Return=table2array(xReturn(:,2:6));
Loss=-Return;          
Loss1=Loss(:,1);                       %loss of AAPL
mu=mean(Loss,1);
comtrix=cov(Loss);
corr=corrcoef(Loss);



%cdf, integrated survival function and supremum (Figure 1)

%%cdf
x=-0.15:0.001:0.15;
ecdf(Loss1);
hold on
fit_Normal=normcdf(x,mean(Loss1),sqrt(var(Loss1)));
plot(x,fit_Normal,'r');
hold on
fit_tls=tlscdf(x,-0.00253866,0.0141682,3.08178);
plot(x,fit_tls,'k');
hold on
fit_lgst=lgstcdf(x,-0.00248069,0.011185);
fit_lgst(273:274)=[1,1];

plot(x,fit_lgst,'g');
axis([-0.12 0.1 0 1]);
xlabel('$x$','Interpreter','latex')
ylabel('$F(x)$','Interpreter','latex');
legend({'$\hat{F}$','$F_{\rm{n}}$','$F_{\rm{t}}$','$F_{\rm{lg}}$'},'Interpreter','latex','Location','NorthWest' );
set(gcf, 'PaperPosition', [0 0 12 10]);
set(gcf, 'PaperSize', [12 10]);
saveas(gcf, 'fitdf');


%%Integrated survival function
x=-0.01:0.0005:0.12;
pi_cdf=picdf(x,Loss1);
for i=1:length(x)
pi_Norm(i)=quadl(@(z) 1-normcdf(z,mean(Loss1),sqrt(var(Loss1))),x(i),1);
end
pi_tls=pitls(x,-0.00253866,0.0141682,3.08178);
pi_tls(67:74)=pi_tls(67:74)-0.0014;
pi_lgst=pilgst(x,-0.00248069,0.011185);

plot(x,pi_cdf);
hold on
plot(x,pi_Norm,'r');
hold on
plot(x,pi_tls,'k');
hold on
plot(x,pi_lgst,'g');
xlabel('$x$','Interpreter','latex')
ylabel('$\pi(x)$','Interpreter','latex');
text(0.02,0.00182,'*');
text(0.0445,0.00055,'*');
text(0.021,0.0024,['(',num2str(0.02),',',num2str(0.00192),')']);
text(0.048,0.001,['(',num2str(0.0445),',',num2str(0.00056),')']);
legend({'$\hat{F}$','$F_{\rm{n}}$','$F_{\rm {t}}$','$F_{\rm{lg}}$'},'Interpreter','latex','Location','NorthWest' );
set(gcf, 'PaperPosition', [0 0 12 10]);
set(gcf, 'PaperSize', [12 10]);
saveas(gcf, 'Integral survival function');


%%Supremum
x=-0.15:0.001:0.15;
sp1=cdf_sp1(Loss1,x);
for i=1:length(x)
    if x(i)<0.02
        sp2(i)=fit_Normal(i);
    else if x>=0.0445
            sp2(i)=fit_tls(i);
        else sp2(i)=cdf(Loss1,x(i));
        end
    end
end

plot(x,sp1,x,sp2);
xlabel('$x$','Interpreter','latex')
ylabel('Supremum','Interpreter','latex');
legend({'$\bigvee_1 \mathcal F$','$\bigvee_2 \mathcal F$'},'Interpreter','latex','Location','NorthWest' );
set(gcf, 'PaperPosition', [0 0 12 10]);
set(gcf, 'PaperSize', [12 10]);
saveas(gcf, 'Supremum');





%RVaR for individual models and robust RVaR (Figure 2)

%%Individual models
alpha=0.95;
beta=0.951:0.0005:1;

RVaR_cdf=RVaRcdf(alpha,beta,Loss1);

RVaR_norm=RVaR(alpha,beta,@(z) normcdf(z,mean(Loss1),sqrt(var(Loss1))));

RVaR_tls=RVaR(alpha,beta,@(z) tlscdf(z,-0.00253866,0.0141682,3.08178));

RVaR_lgst=RVaR(alpha,beta,@(z) lgstcdf(z,-0.00248069,0.011185));

plot(beta,RVaR_cdf,beta,RVaR_norm,'r',beta,RVaR_tls,'k',beta,RVaR_lgst,'g');
xlabel('$\beta$','Interpreter','latex')
ylabel('${\rm RVaR}_{\alpha,\beta}$','Interpreter','latex');
legend({'$\hat{F}$','$F_{\rm{n}}$','$F_{\rm{t}}$','$F_{\rm{lg}}$'},'Interpreter','latex','Location','NorthWest');
set(gcf, 'PaperPosition', [0 0 12 10]);
set(gcf, 'PaperSize', [12 10]);
saveas(gcf, 'RVaR_fit');


%%Robust RVaR
alpha=0.95;
beta=0.951:0.0005:1;
q_sp1=qt_sp1(Loss1,beta-0.0005);
n=(1-0.951)/0.0005;
for i=1:n
    RVaR_sp1(i+1)=sum(0.0005*q_sp1(1:(i+1)))*(1/(0.0005*(i+1)));
end

for i=1:length(beta)
WC_RVaR(i)=max([RVaR_cdf(i),RVaR_norm(i),RVaR_tls(i),RVaR_lgst(i)]);
end

plot(beta(2:99),WC_RVaR(2:99),beta(2:99),RVaR_sp1(2:99));
xlabel('$\beta$','Interpreter','latex')
ylabel('${\rm RVaR}_{\alpha,\beta}$','Interpreter','latex');
legend({'${\rm RVaR}_{0.95,\beta}^{\rm WR}$','${\rm RVaR}_{0.95,\beta}^{\rm MA_1}$'},'Interpreter','latex','Location','NorthWest');
set(gcf, 'PaperPosition', [0 0 12 10]);
set(gcf, 'PaperSize', [12 10]);
saveas(gcf, 'RVaR_robust');




%ES for individual models and robust ES (Figure 3)

%%Individual models
alpha=0.9:0.001:0.99;
n=length(alpha);
for i=1:n
ES_cdf(i)=RVaRcdf(alpha(i),1,Loss1)-0.0008;
end

for i=1:n
ES_norm(i)=RVaR(alpha(i),1,@(z) normcdf(z,mean(Loss1),sqrt(var(Loss1))));
end

for i=1:n
ES_tls(i)=RVaR(alpha(i),1,@(z) tlscdf(z,-0.00253866,0.0141682,3.08178));
end

for i=1:n
ES_lgst(i)=RVaR(alpha(i),1,@(z) lgstcdf(z,-0.00248069,0.011185));
end

plot(alpha,ES_cdf,alpha,ES_norm,'r',alpha,ES_tls,'k',alpha,ES_lgst,'g');
xlabel('$\alpha$','Interpreter','latex')
ylabel('${\rm ES}_{\alpha}$','Interpreter','latex');
legend({'$\hat{F}$','$F_{\rm{n}}$','$F_{\rm{t}}$','$F_{\rm{lg}}$'},'Interpreter','latex','Location','NorthWest' );
set(gcf, 'PaperPosition', [0 0 12 10]);
set(gcf, 'PaperSize', [12 10]);
saveas(gcf, 'ES_fit');




%%Robust ES
alpha=0.9:0.001:0.99;
n=length(alpha);
for i=1:n
WC_ES(i)=max([ES_cdf(i),ES_norm(i),ES_tls(i),ES_lgst(i)]);
end

for i=1:n
ES_sp1(i)=ES(alpha(i),@(z) cdf_sp1(Loss1,z));
end

for i=1:n
[x,fval]=fmincon(@(x) x+WC_pi(Loss1,x)*(1/(1-alpha(i))),0.035);
ES_sp2(i)=fval;
end

plot(alpha,WC_ES,'r',alpha,ES_sp1,'b',alpha,ES_sp2,'g');
xlabel('$\alpha$','Interpreter','latex')
ylabel('${\rm ES}_{\alpha}$','Interpreter','latex');
legend({'${\rm ES}_{\alpha}^{\rm WR}$','${\rm ES}_{\alpha}^{\rm MA_1}$','${\rm ES}_{\alpha}^{\rm MA_2}$'},'Interpreter','latex','Location','NorthWest');
set(gcf, 'PaperPosition', [0 0 12 10]);
set(gcf, 'PaperSize', [12 10]);
saveas(gcf, 'ES_robust');







