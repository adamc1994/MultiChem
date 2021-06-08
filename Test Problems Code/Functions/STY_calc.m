function output = STY_calc(yield,tR,conc,MW)

% calculates STY. 
% tR = residence time. 
% conc = concentration of SM in reactor.
% MW = molecular weight of product.

STY = ((yield.*conc.*MW)./tR)*60; % kg/m3h 
output=STY;