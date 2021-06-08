% CSTR in series model for the SnAr reaction
function [dydt] = SnAr_Equations( ~ , y )
global ko Ea discrete_v n flowrate
conc = reshape ( y(1:n * 5), n, [] ); % 5 = number of chemical species
T = y(n*5 + 1 : n*5 + n-1); % 5 = number of chemical species

% Calculate rate constants, different reactor temperatures allowed
ksq = [ko;ko;ko;ko]; % 4 row vectors for each CSTR
Easq = [Ea;Ea;Ea;Ea]; % 4 row vectors for each CSTR
Tsq = [T T T T]; % 4 T columns for each reaction step
k = ksq .* exp(-Easq ./ (8.314*(Tsq+273.15)));

% SnAr reaction
% Conc at each discretised point 
A = conc(:,1);
B = conc(:,2);
C = conc(:,3);
D = conc(:,4);
E = conc(:,5);

% Calculate the concentration change, isothermal operation assumed
dAdt = [0; % inlet concentration doesn't change
    -flowrate * diff(A) ./ diff(discrete_v) ... % CSTR flowrate adjustment 
    - (k(:,1) .* A(2:end) .* B(2:end) + k(:,2) .* A(2:end) .* B(2:end))]; % rate of reaction for A

dBdt = [0; % inlet concentration doesn't change
    -flowrate * diff(B) ./ diff(discrete_v) ... % CSTR flowrate adjustment 
    - (k(:,1) .* A(2:end) .* B(2:end) + k(:,2) .* A(2:end) .* B(2:end) + k(:,3) .* B(2:end) .* C(2:end) + k(:,4) .* B(2:end) .* D(2:end))]; % rate of reaction for B

dCdt = [0; % inlet concentration doesn't change
    -flowrate * diff(C) ./ diff(discrete_v)... % CSTR flowrate adjustment 
    + (k(:,1) .* A(2:end) .* B(2:end) - k(:,3) .* B(2:end) .* C(2:end))]; % rate of reaction for C

dDdt = [0; % inlet concentration doesn't change
    -flowrate * diff(D) ./ diff(discrete_v) ... % CSTR flowrate adjustment 
    + (k(:,2) .* A(2:end) .* B(2:end) - k(:,4) .* B(2:end) .* D(2:end))]; % rate of reaction for D

dEdt = [0; % inlet concentration doesn't change
    -flowrate * diff(E) ./ diff(discrete_v) ... % CSTR flowrate adjustment 
    + (k(:,3) .* B(2:end) .* C(2:end) + k(:,4) .* B(2:end) .* D(2:end))]; % rate of reaction for E
    
dydt = [reshape([dAdt dBdt dCdt dDdt dEdt],1,[])'; [0 0 0 0]'];  
end