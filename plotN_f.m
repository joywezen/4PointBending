%plot simulation
function []=plotN_f(N_f,loadLevelFactor)
run('/home/SERVER-zhang/Matlab/OneDBarFEM/SFRSCC05data.m')
run('/home/SERVER-zhang/Matlab/external_functions/environment_setting.m')
figure
subplot(1,2,1)
%plot Oh 1986 model
hOh=loglog(concrete.meanweibull,loadLevelFactor,'-s','linewidth',1.5,'color',colorsetred(4,:));
hold on
%plot experiment data
for j=1:size(concrete.data,2)
    hdata=loglog(cell2mat(concrete.data(:,j)),repmat(loadLevelFactor(j),size(cell2mat(concrete.data(:,j)))),'^','color',colorsetred(4,:));
    hold on
end
%plot simulation
hSim=loglog(N_f,loadLevelFactor,'-.o','linewidth',1.5,'color',colorsetblue(6,:));

set(gca,'XScale','log')
grid on
box off
legend([hOh,hdata(1),hSim],{'Oh 1986','Goel 2014','Simu'})
xlabel('$N$','interpreter','latex')
ylabel('load level','interpreter','latex')
subplot(1,2,2)
%relative error
err=abs(concrete.meanweibull-N_f)./N_f*100;
bar(loadLevelFactor,err,'FaceColor',colorsetgreen(5,:))
set ( gca, 'xdir', 'reverse' )
grid on
box off
xlabel('load level','interpreter','latex')
ylabel('$\epsilon\%$ ','interpreter','latex')
legend('relative error')
end