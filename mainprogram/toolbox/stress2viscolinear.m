function [esyxv,esyzv]=stress2viscolinear(stressyx,stressyz,etaeff)
esyxv=2*stressyx./etaeff; % viscous shear strain rate yx component
esyzv=2*stressyz./etaeff; % viscous shear strain rate yz component
end