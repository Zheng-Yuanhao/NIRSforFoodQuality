%load('contourmap_scatterplot-mango_684-1990_10243.mat')
load('contourmap_scatterplot-pharmaceuticaltable-602-600-1622.mat')
contourf(X,Y,Z); hold on
colorbar
x=scatterplot_xyz(:,1);
y=scatterplot_xyz(:,2);
z=scatterplot_xyz(:,3);
plot3(x(1:2),y(1:2),z(1:2),'-ow','LineWidth',1,'MarkerFaceColor','w','MarkerSize',0.5); hold on
plot3(x(99:100),y(99:100),z(99:100),'-ow','LineWidth',1,'MarkerFaceColor','w','MarkerSize',0.5); hold on
plot3(x(2:99),y(2:99),z(2:99),'-ow','LineWidth',1,'MarkerFaceColor','w','MarkerSize',2.5); hold on
plot3(x(1),y(1),z(1),'ro','MarkerFaceColor','r','MarkerSize',6); hold on
plot3(x(size(x,1)),y(size(y,1)),z(size(z,1)),'-ro','LineWidth',1,'MarkerFaceColor','r','MarkerSize',6)
