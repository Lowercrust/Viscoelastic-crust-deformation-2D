function plotshear(shear,ca,model)
% addpath()
figure;
plus=shear;
minus=shear;
plus(plus<0)=0;
minus(minus>0)=0;
plus=-log10(plus);
minus=log10(abs(minus));
% range=[14,ca(2)];
plus(plus>ca(2))=ca(2);
plus=ca(2)-plus;
minus(minus<-ca(2))=-ca(2);
minus=-(minus+ca(2));
log=plus+minus;
% plus(inf)=0;
% minus(inf)=0;
pdeplot(model,'xydata',log) %,'contour','on','levels',10
colormap(brewermap([],'BrBG'))
% ca=range(2)-range(1);
caxis([ca(1)-ca(2) ca(2)-ca(1)]);
axis ij
axis equal tight
end