function [Rsfcl, E, J, Tsc, state] = fcn(current, Ta, Ts)

% superconductor dimensions and thermal properties
sc_length = 50.0;
sc_diameter = 0.004;                            % superconductor diameter (m)
sc_area = pi * (sc_diameter / 2)^2;             % superconductor cross-sectional area (m^2)
Cv = 1.0e6;
k = 1.5e3;
Csc = sc_length * sc_area * Cv;
Rsc = 1 / (k * sc_length * sc_diameter * pi);   % includes the total surface area (minus the "ends")
%nominal_resistance = rho_Tc * sc_length / sc_area;  % SFCL resistance at temperature=Tc

% superconductor characteristics
Jc_77K=1.5e7;               % critical current density at 77K (A/m^2)
E0=0.1;                     % E-field for transition from superconducting state to flux-flow state (V/m)
Ec=1e-6 * 100;              % definition of E-field required for Jc (V/m)
alpha_77K=6.0;              % exponent value during superconducting state (dimensionless)
beta=3;                     % exponent value during flux-flow state (dimensionless)
Tc=95;                      % critical temperature (K)
rho_Tc=1.0e-6;              % superconductor resistivity at Tc (ohm-m^2)

% inputs from previous iterations
persistent T;               % the instantaneous temperature of the superconductor
persistent prevE;           % store the E-field value from the previous iteration
persistent RsfclNew;        % history of Rsfcl values for averaging
persistent Qremoved;
persistent Qsc;
persistent time;            % elapsed simulation time

if isempty(T)
    T = Ta;                 % Ta is the ambient temperature of the coolant; it is assumed that this does not change
end
if isempty(prevE)
    prevE = 0.0;
end
if isempty(time)
    time = 0.0;
end
if isempty(RsfclNew)
    RsfclNew = [1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6];
end
if isempty(Qremoved)
    Qremoved = 0.0;
end
if isempty(Qsc)
    Qsc = 0.0;
end

RATE_LIMITING = 0;              % set to 1 to enable
AVERAGING_FILTER = 0;           % a value >0 sets the number of samples used for the averaging filter; a value of 0 disables; max value is 10
time = time + 1;
TEMPERATURE_TIMESTEPS = 1;      % defines number of timesteps needed for each thermal calculation; value of 1 means every timestep


J = abs(current) / sc_area;     % inputs from power system

if T < Tc
    if (prevE < E0)
        % superconducting
        state = 0;
    else
        % flux-flow
        state = 1;
    end
    
    Jc = Jc_77K * ((Tc - T) / (Tc - 77));
else
    % normal conducting
    state = 2;
    
    Jc = Jc_77K;    % value of Jc does not matter if T >= Tc
end

alpha_alt = log10(E0 / Ec) / log10( (Jc_77K / Jc)^(1 - (1/beta)) * (E0 / Ec)^(1 / alpha_77K) );
alpha = max(beta, alpha_alt);

% todo: prevent state from changing every timestep
if state == 0
    % superconducting
    E=Ec * (J / Jc)^alpha;
elseif state == 1
    % flux-flow
    E=E0 * (Ec / E0)^(beta/alpha_77K) * (Jc_77K / Jc) * (J / Jc_77K)^beta;
else
    % normal conducting
    E=rho_Tc * (T / Tc) * J;
end

prevE = E;

RsfclNew = [0, RsfclNew(1), RsfclNew(2), RsfclNew(3), RsfclNew(4), RsfclNew(5), RsfclNew(6), RsfclNew(7), RsfclNew(8), RsfclNew(9)];    % todo: convert to array of any length

if abs(current) > 0
    RsfclNew(1) = E * sc_length / abs(current);
else
    RsfclNew(1) = RsfclNew(2);
end

if (RATE_LIMITING == 1)
    %rate limiting; hold Rsfcl value at max rate, if rate is outwith limit
    RATE_LIMIT = 1.5;
    RATE_COMPARE = RsfclNew(2);
    if RsfclNew(1) > RATE_LIMIT * RATE_COMPARE
        RsfclNew(1) = RATE_LIMIT * RATE_COMPARE;
    elseif RsfclNew(1) < RATE_COMPARE / RATE_LIMIT
        RsfclNew(1) = RATE_COMPARE / RATE_LIMIT;
    end
end

if (AVERAGING_FILTER > 0)
    RsfclNew(1) = mean(RsfclNew(1:AVERAGING_FILTER));     % help smooth-out variations in Rsfcl
end
Rsfcl = RsfclNew(1);

% assume that if the flux-flow resistance >= the equivalent normal state resitance, then the superconductor has quenched
if Rsfcl >= (rho_Tc * (T / Tc) * J * sc_length / abs(current))
    state = 2;
    Rsfcl = (rho_Tc * (T / Tc) * J * sc_length / abs(current));
end

%QscDiff = current^2 * RsfclNew(1) * Ts;
QscDiff = current^2 * Rsfcl * Ts;
QremovedDiff = ((T - Ta) / Rsc) * Ts;
Qsc = Qsc + QscDiff;
Qremoved = Qremoved + QremovedDiff;

if (mod(time, TEMPERATURE_TIMESTEPS) == 0)
    T = T + (Qsc - Qremoved) / Csc;         % slight difference from [Langston2005]: Ta replaced by T
    Qsc = 0.0;
    Qremoved = 0.0;
end

Tsc = T;
