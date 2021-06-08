% CSTR in series model for the Paal-Knorr reaction
function [dydt] = PK_Equations( ~ , y )
global ko Ea discrete_v n flowrate
conc = reshape ( y(1:n * 4), n, [] );
T = y(n*4 + 1 : n*4 + 4);

% Calculate rate constants, different reactor temperatures allowed
ksq = [ko;ko;ko;ko]; % 4 row vectors for each CSTR
Easq = [Ea;Ea;Ea;Ea]; % 4 row vectors for each CSTR
Tsq = [T T]; % 2 T columns for each reaction step
k = ksq .* exp(-Easq ./ (8.314*(Tsq+273.15)));

% Paal-Knorr reaction
% Conc at each discretised point 
A = conc(:,1);
B = conc(:,2);
C = conc(:,3);
D = conc(:,4);

% Calculate the concentration change, isothermal operation assumed
dAdt = [0; % inlet concentration doesn't change
    -flowrate * diff(A) ./ diff(discrete_v) ... % CSTR flowrate adjustment 
    - (k(:,1) .* A(2:end) .* B(2:end))]; % rate of reaction for A

dBdt = [0; % inlet concentration doesn't change
    -flowrate * diff(B) ./ diff(discrete_v) ... % CSTR flowrate adjustment 
    - (k(:,1) .* A(2:end) .* B(2:end))]; % rate of reaction for B

dCdt = [0; % inlet concentration doesn't change
    -flowrate * diff(C) ./ diff(discrete_v)... % CSTR flowrate adjustment 
    + (k(:,1) .* A(2:end) .* B(2:end) - k(:,2) .* C(2:end))]; % rate of reaction for C

dDdt = [0; % inlet concentration doesn't change
    -flowrate * diff(D) ./ diff(discrete_v) ... % CSTR flowrate adjustment 
    + (k(:,2) .* C(2:end))]; % rate of reaction for D
    
dydt = [reshape([dAdt dBdt dCdt dDdt ],1,[])'; [0 0 0 0]'];  
end