%% Function to check the neighbourhood 
function checkFamily(x,familyMat,ii)
family = familyMat(ii,:);
ii_neigh = family(family~=0);
plotname = "Family of the " + int2str(ii) + "-th node";
% Plot
figure
scatter(x(:,1),x(:,2),'b','filled')
hold on
scatter(x(ii_neigh,1),x(ii_neigh,2),'r','filled')
scatter(x(ii,1),x(ii,2),'g','filled')
xlabel('x [m]')
ylabel('y [m]')
title(plotname)
set(gca,'FontSize',15)
end