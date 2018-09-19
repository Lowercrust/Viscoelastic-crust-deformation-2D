function plotfaultstress(yq1,savefaultstress,a,b)
for i=a:b
plot(yq1,savefaultstress{i})
hold on
end
end