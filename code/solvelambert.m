% -------------------------------------------------------------------------

function [V1, V2] = solvelambert(R1, R2, t, proRetro, mu)

% This function solves Lambert’s problem.

%code from lamberto dell'elce
% 
% mu         - gravitational parameter (kmˆ3/sˆ2)
% R1, R2     - initial and final position vectors (km)
% r1, r2     - magnitudes of R1 and R2
% t          - the time of flight from R1 to R2 (a constant) (s)
% V1, V2     - initial and final velocity vectors (km/s)
% c12        - cross product of R1 into R2
% theta      - angle between R1 and R2
% proRetro   - 1 if the orbit is prograde, -1 if the orbit is
%              retrograde
% A          - a constant given by Equation 5.35
% z          - alpha*xˆ2, where alpha is the reciprocal of the
%              semimajor axis and x is the universal anomaly
% y(z)       - a function of z given by Equation 5.38
% F(z,t)     - a function of the variable z and constant t, given by
%              Equation 5.40
% dFdz(z)    - the derivative of F(z,t), given by Equation 5.43
% ratio      - F/dFdz
% tol        - tolerance on precision of convergence
% nmax       - maximum number of iterations of Newton’s procedure
% f, g       - Lagrange coefficients
% gdot       - time derivative of g
% C(z), S(z) - Stumpff functions
% dum        - a dummy variable

N = length(t);


% Magnitudes of R1 and R2:
r1 = sqrt(sum(R1.^2, 2));
r2 = sqrt(sum(R2.^2, 2));
c12 = vect(R1, R2);
theta = acos(sum(R1 .* R2, 2) ./ r1 ./ r2);

% Determine whether the orbit is prograde or retrograde:
aux = (c12(:, 3) < 0 & proRetro == 1);
theta(aux) = 2 * pi - theta(aux);

aux = (c12(:, 3) >= 0 & proRetro == - 1);
theta(aux) = 2 * pi - theta(aux);

% Equation 5.35:
A = sin(theta) .* sqrt(r1 .* r2 ./ (1 - cos(theta)));

% Determine approximately where F(z,t) changes sign, and
% use that value of z as the starting value for Equation 5.45:
z = - 1000 * ones(N, 1);
F = - 1 * ones(N, 1);
idxNeg = (F < 0);
while(any(idxNeg))
    z(idxNeg) = z(idxNeg) + 0.1;
    [~, F] = equations(z, t, r1, r2, mu, A);
    idxNeg = (F < 0 & ~isnan(F));
end

% Set an error tolerance and a limit on the number of iterations:
tol = 1.e-8;
nmax = 5000;

% Iterate on Equation 5.45 until z is determined to within the error
% tolerance:
ratio = 1;
n = 0;

while (max(abs(ratio)) > tol) && (n <= nmax),
    n = n + 1;
    [y, F, dFdz] = equations(z, t, r1, r2, mu, A);
    F;
    ratio = F ./ dFdz;
    z = z - ratio;
    size(z);
    max(abs(ratio));
end


% Report if the maximum number of iterations is exceeded:
if n >= nmax
    fprintf('\n\n  Number of iterations exceeds %g \n\n ', nmax)
end

% Equation 5.46a:
f = 1 - y ./ r1;

% Equation 5.46b:
g = A .* sqrt(y / mu);

% Equation 5.46d:
gdot = 1 - y ./ r2;

% Equation 5.28:
V1 = (R2 - repmat(f, 1, 3) .* R1) ./ repmat(g, 1, 3);

% Equation 5.29:
V2 = (repmat(gdot, 1, 3) .* R2 - R1) ./ repmat(g, 1, 3);

V1 = real(V1); V2 = real(V2);

return

% -------------------------------------------------------------------------

function [y, F, dFdz] = equations(z, t, r1, r2, mu, A)

S = stumpS(z);
C = stumpC(z);


y = r1 + r2 + A .* (z .* S - 1) ./ sqrt(C);
F = (y ./ C).^1.5 .* S + A .* sqrt(y) - sqrt(mu) .* t;

dFdz = zeros(size(z));

idx0 = (z == 0);
dFdz(idx0) = sqrt(2) / 40 * y(idx0).^1.5 + ...
    A(idx0) / 8 .* (sqrt(y(idx0)) + A(idx0) .* sqrt(0.5 ./ y(idx0)));

dFdz(~idx0) = (y(~idx0) ./ C(~idx0)).^1.5 .* (0.5 ./ z(~idx0) .* ...
    (C(~idx0) - 1.5 * S(~idx0) ./ C(~idx0)) + ...
    3 / 4 * S(~idx0).^2 ./ C(~idx0)) ...
    + A(~idx0) / 8 .* (3 * S(~idx0) ./ C(~idx0) .* sqrt(y(~idx0)) ...
    + A(~idx0) .* sqrt(C(~idx0) ./ y(~idx0)));

return

% -------------------------------------------------------------------------

function s = stumpS(z)

s = ones(size(z)) / 6;

idxPos = (z > 0);
idxNeg = (z < 0);

s(idxPos) = (sqrt(z(idxPos)) - sin(sqrt(z(idxPos)))) ./ ...
    (sqrt(z(idxPos))).^3;
s(idxNeg) = (sinh(sqrt(- z(idxNeg))) - sqrt(- z(idxNeg))) ./ ...
    (sqrt(- z(idxNeg))).^3;

return

% -------------------------------------------------------------------------

function c = stumpC(z)

c = ones(size(z)) / 2;

idxPos = (z > 0);
idxNeg = (z < 0);

c(idxPos) = (1 - cos(sqrt(z(idxPos)))) ./ z(idxPos);
c(idxNeg) = (cosh(sqrt(- z(idxNeg))) - 1) ./ (- z(idxNeg));

return

% -------------------------------------------------------------------------

function c = vect(a, b)

c = zeros(size(a));
c(:, 1) = a(:, 2) .* b(:, 3) - a(:, 3) .* b(:, 2);
c(:, 2) = a(:, 3) .* b(:, 1) - a(:, 1) .* b(:, 3);
c(:, 3) = a(:, 1) .* b(:, 2) - a(:, 2) .* b(:, 1);

return

% -------------------------------------------------------------------------

function JD = gregorian2julian(YYYY, MM, DD, hh, mm, ss)

JD = 367 * YYYY ...
    - floor(7 / 4 * (YYYY + floor((MM + 9) / 12))) ...
    + floor(275 * MM / 9) ...
    + DD + 1721013.5 ...
    + (hh + (mm + ss / 60) / 60) / 24;
return

% -------------------------------------------------------------------------

function [year, month, day, h, m, s, dategreg] = julian2gregorian(JD)

I = floor(JD + 0.5);
Fr = abs(I - (JD + 0.5));	 

B = I;
aux = I >= 2299160;
A = floor((I(aux) - 1867216.25 ) / 36524.25);
a4 = floor(A / 4);
B(aux) = I(aux) + 1 + A - a4;
C = B + 1524;
D = floor((C - 122.1) / 365.25);
E = floor(365.25 * D);
G = floor((C - E) / 30.6001);

day = floor(C - E + Fr - floor(30.6001 * G));

month = G - 1;
month(G > 13.5) = month(G > 13.5) - 12;

year = D - 4716;
year(month <= 2.5) = year(month <= 2.5) + 1;

h = floor(Fr * 24);

m = floor(abs(h - (Fr * 24)) * 60);

minufrac = (abs(h - (Fr * 24)) * 60); 
s = ceil(abs(m - minufrac) * 60);

dategreg = [year, month, day, h, m, s];

return

% -------------------------------------------------------------------------

function [r, v] = keplerian2rv(KP, mu)

% 
% This function computes the initial position and velocity (from the input
% keplerian parameters).
% Input row vector of N elements.
% Outputs matrices of 3 * N elements
% 
% Lamberto Dell'Elce
% 

a = KP(:, 1)';
e = KP(:, 2)';
i = KP(:, 3)';
omega = KP(:, 4)';
RAAN = KP(:, 5)';
theta = KP(:, 6)';

N = length(a);

sTheta = sin(theta);
cTheta = cos(theta);

r_orb = ([1; 1; 1] * (a .* (1 - e.^2) ./ (1 + e .* cTheta))) .* ...
    [cTheta; sTheta; zeros(1, N)];
v_orb = ([1; 1; 1] * sqrt(mu ./ a ./ (1 - e.^2))) .* ...
    [- sTheta; e + cTheta; zeros(1, N)];

r = zeros(3, N);
v = zeros(3, N);
sRAAN_vec = sin(RAAN);
cRAAN_vec = cos(RAAN);
sI_vec = sin(i);
cI_vec = cos(i);
sOmega_vec = sin(omega);
cOmega_vec = cos(omega);
for j = 1 : N,
    sRAAN = sRAAN_vec(j);
    cRAAN = cRAAN_vec(j);
    sI = sI_vec(j);
    cI = cI_vec(j);
    sOmega = sOmega_vec(j);
    cOmega = cOmega_vec(j);
    
    ARAAN = [cRAAN, sRAAN, 0; - sRAAN, cRAAN, 0; 0, 0, 1];
    Ai = [1, 0, 0; 0, cI, sI; 0, - sI, cI];
    Aomega = [cOmega, sOmega, 0; - sOmega, cOmega, 0; 0, 0, 1];
    
    Atot = ARAAN' * Ai' * Aomega';
    
    r(:, j) = Atot * r_orb(:, j); % [m]
    v(:, j) = Atot * v_orb(:, j); % [m / s]
end

r = r';
v = v';

return

% -------------------------------------------------------------------------

function [coe, r, v] = planet_elements_and_sv(planet_id, jd)

% This function calculates the orbital elements and the state
% vector of a planet from the date (year, month, day)
% and universal time (hour, minute, second).
% 
% coe  - vector of heliocentric orbital elements
%        [a e RA incl w TA h w_hat L M E]
%        angles in deg, distances in km
% 
% User M-functions required: kepler_E,

mu = 1.32712440018e+11;

deg = pi/180;

% Obtain the data for the selected planet from Table 8.1:
[J2000_coe, rates] = planetary_elements(planet_id);

% Equation 8.104a:
t0 = (jd - 2451545) / 36525;

% Equation 8.104b:
elements = ones(size(t0)) * J2000_coe + t0 * rates;
a = elements(:, 1);
e = elements(:, 2);

% Equation 2.61:
h = sqrt(mu * a .* (1 - e.^2));

% Reduce the angular elements to within the range 0 - 360 degrees:
incl = elements(:, 3);
RA = zero_to_360(elements(:, 4));
w_hat = zero_to_360(elements(:, 5));
L = zero_to_360(elements(:, 6));
w = zero_to_360(w_hat - RA);
M = zero_to_360((L - w_hat));

% Algorithm 3.1 (for which M must be in radians)
E = kepler_E(e, M*deg);

% Equation 3.10 (converting the result to degrees):
TA = zero_to_360(2 * atan(sqrt((1 + e) ./(1 - e)) .* tan(E / 2)) / deg);
coe = [a, e, RA, incl, w, TA, h, w_hat, L, M, E / deg];

% Algorithm 4.2 (for which all angles must be in radians):
[r, v] = keplerian2rv([a e incl*deg w*deg RA*deg TA*deg], mu);

return

% -------------------------------------------------------------------------

function [J2000_coe, rates] = planetary_elements(planet_id)

aux = [ 0.38709893, 0.20563069, 7.00487, 48.33167, 77.45450, 252.25084, 0.00000066, 0.00002527, -23.51, -446.30, 573.57, 538101628.29;
0.72333199, 0.00677323, 3.39471, 76.68069, 131.53298, 181.97973, 0.00000092, -0.00004938, -2.86, -996.89, -108.80, 210664136.06;
1.00000011, 0.01671022, 0.00005, -11.26064, 102.94719, 100.46435, -0.00000005, -0.00003804, -46.94, -18228.25, 1198.28, 129597740.63;
1.52366231, 0.09341233, 1.85061, 49.57854, 336.04084, 355.45332, -0.00007221, 0.00011902, -25.47, -1020.19, 1560.78, 68905103.78;
5.20336301, 0.04839266, 1.30530, 100.55615, 14.75385, 34.40438, 0.00060737, -0.00012880, -4.15, 1217.17, 839.93, 10925078.35;
9.53707032, 0.05415060, 2.48446, 113.71504, 92.43194, 49.94432, -0.00301530, -0.00036762, 6.11, -1591.05, -1948.89, 4401052.95;
19.19126393, 0.04716771, 0.76986, 74.22988, 170.96424, 313.23218, 0.00152025, -0.00019150, -2.09, -1681.40, 1312.56, 1542547.79;
30.06896348, 0.00858587, 1.76917, 131.72169, 44.97135, 304.88003, -0.00125196, 0.00002514, -3.64, -151.25, -844.43, 786449.21;
39.48168677, 0.24880766, 17.14175, 110.30347, 224.06676, 238.92881, -0.00076912, 0.00006465, 11.07, -37.33, -132.25, 522747.90];

J2000_elements = aux(:, 1 : 6);
cent_rates = aux(:, 7 : end);

J2000_coe = J2000_elements(planet_id,:);
rates = cent_rates(planet_id,:);

% Convert from AU to km:
au = 149597871;
J2000_coe(1) = J2000_coe(1) * au;
rates(1) = rates(1) * au;

% Convert from arcseconds to fractions of a degree:
rates(3:6) = rates(3:6)/3600;

return

% -------------------------------------------------------------------------

function y = zero_to_360(x)

y = x;
aux = (x >= 360);
y(aux) = x(aux) - fix(x(aux) / 360) * 360;

aux = (x < 0);
y(aux) = x(aux) - (fix(x(aux) / 360) - 1) * 360;

return

% -------------------------------------------------------------------------

function E = kepler_E(e, M)

% Set an error tolerance:
error = 1.e-8;

%...Select a starting value for E:
E = M + e / 2;
E(M >= pi) = M(M >= pi) - e(M >= pi) / 2;

% Iterate on Equation 3.14 until E is determined to within
% the error tolerance:
ratio = 1;

while max(abs(ratio)) > error,
    ratio = (E - e .* sin(E) - M) ./ (1 - e .* cos(E));
    E = E - ratio;
end

return


