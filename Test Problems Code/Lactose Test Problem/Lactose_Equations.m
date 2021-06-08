% CSTR in series model for the lactose isomerisation
function [dydt] = Lactose_Equations( ~ , y )
global ko Ea discrete_v n flowrate
conc = reshape ( y(1:n * 3), n, [] ); % 3 = number of chemical species
T = y(n*3 + 1 : n*3 + n-1); % 3 = number of chemical species

% Calculate rate constants, different reactor temperatures allowed
ksq = [ko;ko;ko;ko]; % 4 row vectors for each CSTR
Easq = [Ea;Ea;Ea;Ea]; % 4 row vectors for each CSTR
Tsq = [T T T]; % 3 T columns for each reaction step
k = ksq .* exp(-Easq ./ (8.314*(Tsq+273.15)));

% Lactose isomerisation reaction
% Conc at each discretised point 
A = conc(:,1);
B = conc(:,2);
C = conc(:,3);

% Calculate the concentration change, isothermal operation assumed
dAdt = [0; % inlet concentration doesn't change
    -flowrate * diff(A) ./ diff(discrete_v) ... % CSTR flowrate adjustment 
    - (k(:,1) .* A(2:end) + k(:,2) .* A(2:end))]; % rate of reaction for A

dBdt = [0; % inlet concentration doesn't change
    -flowrate * diff(B) ./ diff(discrete_v) ... % CSTR flowrate adjustment 
    - (k(:,3) .* B(2:end) - k(:,1) .* A(2:end))]; % rate of reaction for B

dCdt = [0; % inlet concentration doesn't change
    -flowrate * diff(C) ./ diff(discrete_v)... % CSTR flowrate adjustment 
    + (k(:,2) .* A(2:end) + k(:,3) .* B(2:end))]; % rate of reaction for C
    
dydt = [reshape([dAdt dBdt dCdt],1,[])'; [0 0 0 0]'];  
end