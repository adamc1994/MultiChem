function [Cfinal, Yield, TimeData, ConcData] = SnAr_Constants(ConcA, ConcB, T, res_time)
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
Co = [ConcA ConcB zeros(1,3)]; % initial conc at reactor inlet Ca0 Cb0 Cc0...
Co = [Co; zeros(n-1,5)];
Co = reshape (Co, 1, []);


Ea = [4.32e4 3.53e4 4.08e4 6.89e4]; % J/mol
ko = [1.5597e6 13.9049e3 10.4046e3 370.3652e6]; % min-1 for first order, M-1 min-1 for second order

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
    ode45(@SnAr_Equations,[0 res_time*4],yo, options);
%{
    ConcData comes out in the form of:
    
    Columns 1 : n are the concentrations of A at each discrete section
    n+1 : 2*n are the concentrations of B at each discrete section
    Etc...
    Columns no_of_reactants * n : no_of_reactants * n + (n - 1) are the
    temperatures of each CSTR reactor
    %}

%{
---------------------------------------------------------------------------
    Calculate performance parameters
---------------------------------------------------------------------------
    %}

%                 A                 B                    C                D            E
Cfinal = [ConcData(end,n) ConcData(end,2*n) ConcData(end,3*n) ConcData(end,4*n) ConcData(end,5*n)];

Yield = Cfinal / ConcData(1,1) + 1e-6; % correcting to ensure non zero value when using TSEMO
