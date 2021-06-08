function [Cfinal, Yield, TimeData, ConcData] = PK_Constants(ConcA, ConcB, T, res_time)
%{
function [Cfinal, Yield, Efactor, STY, TimeData, ConcData] = series_cstr(ConcA, T, res_time)
Determines the final concentration and yields of the exiting
reactants/products, also determines the E factor and space time yield for
the reaction at the user specified conditions
---------------------------------------------------------------------------
INPUT
ConcA       Scalar input for the feed concentration of reactant A entering
            the first CSTR
T           Vector input for the temperatures of each CSTR reactor, this is
            used in calculating the rates of reaction for each reactor. Reactors
            assumed to operate isothermally
res_time    a

OUTPUT
Cfinal      a

Yield       a

Efactor     a

STY         a

TimeData    a

ConcData    a
---------------------------------------------------------------------------
    Set constants for the PRF reactor
---------------------------------------------------------------------------
    %}
global Ea ko discrete_v flowrate n
n = 4 + 1;
V = 8;
% discretise volume elements for each segment
discrete_v = linspace(0, V, n)';
Co = [ConcA ConcB zeros(1,2)]; % initial conc
Co = [Co; zeros(n-1,4)];
Co = reshape (Co, 1, []);


Ea = [12.2e3 20.0e3]; % J/mol
ko = [15.40 405.19]; % min-1 for first order, M-1 min-1 for second order

flowrate = V / res_time;
yo = [Co T];
%{
---------------------------------------------------------------------------
    Perform numerical intergration
---------------------------------------------------------------------------
    %}
options = odeset('NonNegative', 1:length(Co) );
 
options = odeset(options, 'RelTol', 1e-3);

[TimeData,ConcData] = ...
    ode45(@PK_Equations,[0 res_time*4],yo, options);
%{
---------------------------------------------------------------------------
    Calculate performance parameters
---------------------------------------------------------------------------
    %}

%                 A                 B                    C                D
Cfinal = [ConcData(end,n) ConcData(end,2*n) ConcData(end,3*n) ConcData(end,4*n)];

Yield = Cfinal / ConcData(1,1) + 1e-6; % correcting to ensure non zero value when using TS-EMO
