lgnd=legend('Initial','Estimated','Exact');
set(lgnd,'FontSize',18);

lgnd=legend('Exact','Filter','smoother','Back sampling');
set(lgnd,'FontSize',18);
title(['mode ( ', num2str(ix),' , ', num2str(iy), ' )'],'FontSize',16)

lgnd=legend('Filter','Smoother','Back sampling');
set(lgnd,'FontSize',18);
