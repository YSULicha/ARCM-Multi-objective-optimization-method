%% ARCM: Multi-Objective Optimization Method for Energy Systems
%% DOI: 10.1109/TPWRS.2026.3682532

%% 1. Setup & Parameters
clc; 
clear;
close all;

T = 24; % Time periods
cf = 0.025; % Natural gas price
cb = [0.05, 0.05, 0.04, 0.04, 0.05, 0.06, 0.08, 0.12, 0.15, 0.16, 0.15, 0.14, 0.14, 0.12, 0.10, 0.12, 0.15, 0.18, 0.19, 0.18, 0.16, 0.12, 0.08, 0.06]; % Grid electricity buy price
cs = [0.03, 0.03, 0.03, 0.03, 0.04, 0.04, 0.05, 0.05, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05, 0.04, 0.04, 0.05, 0.06, 0.07, 0.07, 0.06, 0.05, 0.04, 0.04]; % Grid electricity sell price
clc = 0.18; % Load curtailment compensation
eg = [0.55, 0.53, 0.52, 0.50, 0.51, 0.55, 0.60, 0.65, 0.68, 0.65, 0.60, 0.55, 0.45, 0.40, 0.38, 0.35, 0.40, 0.50, 0.65, 0.75, 0.72, 0.68, 0.60, 0.58]; % Grid carbon emission factors
ef = [0.20, 0.20, 0.20, 0.21, 0.21, 0.22, 0.22, 0.23, 0.24, 0.25, 0.25, 0.25, 0.24, 0.24, 0.23, 0.23, 0.22, 0.21, 0.21, 0.20, 0.20, 0.20, 0.20, 0.20]; % Natural gas carbon emission factors
cbc = 0.03; % BESS degradation cost

Le = [45, 42, 40, 41, 43, 50, 65, 88, 95, 92, 90, 88, 85, 86, 89, 93, 105, 120, 135, 125, 110, 90, 70, 55] * 10; % Electrical load
Lc = [5, 4, 3, 3, 4, 8, 15, 25, 35, 45, 58, 65, 70, 72, 70, 65, 55, 40, 30, 20, 15, 10, 8, 6] * 10; % Cooling load
Ppv  = [0, 0, 0, 0, 0, 3, 12, 25, 40, 48, 55, 58, 55, 49, 41, 28, 14, 4, 0, 0, 0, 0, 0, 0] * 10; % PV generation
Pwt  = [25, 28, 32, 30, 26, 22, 18, 15, 19, 23, 28, 35, 38, 36, 31, 27, 24, 29, 33, 37, 40, 36, 32, 28] * 10; % Wind generation
Tout = [2, 1, 0, 0, 1, 3, 6, 9, 12, 15, 18, 20, 21, 22, 21, 19, 16, 13, 10, 8, 6, 4, 3, 2]; % Outdoor temperature

Lh = [70, 68, 65, 66, 70, 78, 85, 90, 82, 75, 68, 65, 64, 66, 70, 75, 88, 100, 110, 105, 95, 85, 78, 72] * 10; % Total heat load
Lhp = 0.9 * Lh; % Process heat load

lsr = 0.2; % Load shifting ratio
lcr = 0.05; % Load curtailment ratio

PcM = 800; Pcm = 200; % Max/Min CHP power
HgM = 1000; % Max gas boiler power
PgM = 1500; % Max grid power
Sm = 2000; Sn = 400; S0 = 1000; % Max/Min/Init BESS capacity
PbcM = 500; PbdM = 500; % Max BESS charge/discharge power
PgtM = 400; % Max GT power

qa = 0.004; qb = 1.9; qc = 6; % CHP quadratic cost coefficients
Ppts = linspace(Pcm, PcM, 7); % CHP PWL power points
Fpts = qa * (Ppts.^2) + qb * Ppts + qc; % CHP PWL fuel points
Nseg = length(Ppts) - 1; % CHP PWL segments

rhp = 1.2; % CHP heat-to-power ratio
ehg = 0.9; % Gas boiler efficiency
ec = 0.95; ed = 0.95; % BESS charge/discharge efficiency
cop = 3.5; % EC coefficient of performance
egt = 0.35; % GT efficiency
TiM = 24; Tim = 19.5; Ti0 = 22; % Max/Min/Init indoor temperature
Ca = 50; Ra = 4.0; % Thermal capacitance and resistance
Af = exp(-1/(Ra * Ca)); Bf = 1 - Af; % Thermal dynamic coefficients

%% 2. Optimization Variables & Constraints
yalmip('clear');
Hc = sdpvar(T,1); Hg = sdpvar(T,1); Pb = sdpvar(T,1); Ps = sdpvar(T,1);
Pbc = sdpvar(T,1); Pbd = sdpvar(T,1); S = sdpvar(T,1); Pe = sdpvar(T,1);
Pld = sdpvar(T,1); Plu = sdpvar(T,1); Plc = sdpvar(T,1); Pgt = sdpvar(T,1);
Tin = sdpvar(T,1); Hs = sdpvar(T,1); Hps = sdpvar(T,1); 
lam = sdpvar(length(Ppts), T, 'full');
zc = binvar(Nseg, T, 'full');
zs = binvar(T,1); 

Fg = Hg / ehg;
Fcp = Fpts * lam; 
Fgt = Pgt / egt; 

f1 = sum(cf * Fcp) + sum(cf * Fg) + cb * Pb - cs * Ps + sum(clc * Plc) + sum(cbc * (Pbc + Pbd)) + sum(cf * Fgt);
f2 = sum(ef .* Fcp) + sum(ef' .* Fg) + eg * Pb + sum(ef' .* Fgt);
          
cons = [];
cons = [cons, S(1) == S0 + Pbc(1)*ec - Pbd(1)/ed];
cons = [cons, Tin(1) == Ti0 * Af + (Hs(1) * Ra + Tout(1)) * Bf];
for t = 2:T
    cons = [cons, Tin(t) == Tin(t-1) * Af + (Hs(t) * Ra + Tout(t)) * Bf];
end
cons = [cons, Tim <= Tin <= TiM];
Pcp = Ppts * lam; 
for t = 1:T
    L_act = Le(t) - Pld(t) + Plu(t) - Plc(t);
    cons = [cons, Pcp(t) + Ppv(t) + Pwt(t) + Pb(t) + Pbd(t) + Pgt(t) == L_act + Ps(t) + Pbc(t) + Pe(t)];
    cons = [cons, Hc(t) + Hg(t) == Hs(t) + Hps(t)];
    cons = [cons, Hps(t) == Lhp(t)];
    cons = [cons, Pe(t) * cop == Lc(t)];
    cons = [cons, Hc(t) <= rhp * Pcp(t)];
    if t > 1
        cons = [cons, S(t) == S(t-1) + Pbc(t)*ec - Pbd(t)/ed];
    end
end
cons = [cons, 0 <= Hg <= HgM, 0 <= Pb <= PgM, 0 <= Ps <= PgM];
cons = [cons, 0 <= Pbc <= PbcM, 0 <= Pbd <= PbdM];
cons = [cons, Sn <= S <= Sm, S(T) == S0];
cons = [cons, sum(Pld) == sum(Plu)];
cons = [cons, 0 <= Pld <= lsr * Le', 0 <= Plu <= lsr * Le', 0 <= Plc <= lcr * Le'];

Ms = 5000; 
cons = [cons, Pld <= Ms * zs, Plu <= Ms * (1 - zs)]; 
cons = [cons, 0 <= lam <= 1];
for t = 1:T
    cons = [cons, sum(lam(:,t)) == 1, sum(zc(:,t)) == 1];
    cons = [cons, lam(1,t) <= zc(1,t)];
    for i = 2:Nseg
        cons = [cons, lam(i,t) <= zc(i-1,t) + zc(i,t)];
    end
    cons = [cons, lam(end,t) <= zc(end,t)];
end
cons = [cons, 0 <= Pgt <= PgtM, 0 <= Hc, Hs >= 0];

%% 3. ARCM Solver
ops = sdpsettings('solver', 'gurobi', 'verbose', 0, 'gurobi.MIPGap', 1e-4);

sol1 = optimize(cons, f1, ops);
f1_min = value(f1); f2_1 = value(f2);

sol2 = optimize(cons, f2, ops);
f2_min = value(f2); f1_2 = value(f1);

Ns = 20; 
pts = [f1_2, f2_min; f1_min, f2_1];
df1 = f1_2 - f1_min; df2 = f2_1 - f2_min;
grid = linspace(0, 1, Ns + 2);

for i = 1:Ns
    if mod(i, 2) == 1 
        c1 = f1_min + grid(i+1) * df1;
        sA = optimize([cons, f1 <= c1], f2, ops);
        if sA.problem == 0
            sB = optimize([cons, f2 <= value(f2)], f1, ops);
            if sB.problem == 0, pts = [pts; value(f1), value(f2)]; end
        end
    else 
        c2 = f2_1 - grid(i+1) * df2;
        sA = optimize([cons, f2 <= c2], f1, ops);
        if sA.problem == 0
            sB = optimize([cons, f1 <= value(f1)], f2, ops);
            if sB.problem == 0, pts = [pts; value(f1), value(f2)]; end
        end
    end
end

u_pts = unique(round(pts, 4), 'rows');
[~, idx] = sort(u_pts(:,1));
s_pts = u_pts(idx, :);

%% 4. Command-Line Results Display
fprintf('--- "What-if" Scenario Analysis: Cost-Emission Trade-off ---\n');
fprintf(repmat('=', 1, 110)); fprintf('\n');
fprintf('%-8s %-18s %-22s %-18s %-22s %-25s\n', 'Point', 'Total Cost ($)', 'Total Emissions (kg)', 'Cost Increase ($)', 'Emission Decrease (kg)', 'Marginal Cost ($/kg CO2)');
fprintf(repmat('-', 1, 110)); fprintf('\n');
fprintf('%-8d %-18.2f %-22.2f %-18s %-22s %-25s\n', 1, s_pts(1,1), s_pts(1,2), 'N/A (Baseline)', 'N/A (Baseline)', 'N/A');

for i = 2:size(s_pts, 1)
    dc = s_pts(i,1) - s_pts(i-1,1);
    de = s_pts(i-1,2) - s_pts(i,2); 
    if de > 1e-6
        mr = sprintf('%.3f', dc / de);
    else
        mr = 'Inf';
    end
    fprintf('%-8d %-18.2f %-22.2f %-18.2f %-22.2f %-25s\n', i, s_pts(i,1), s_pts(i,2), dc, de, mr);
end
fprintf(repmat('=', 1, 110)); fprintf('\n\n');

%% 5. Pareto Frontier Plotting
figure('Name', 'IES Pareto Frontier (ARCM)', 'Position', [100, 100, 700, 550]);
hold on; box on;
plot(s_pts(:,1), s_pts(:,2), 'o-', 'LineWidth', 2.5, 'Color', '#4198AC', 'MarkerSize', 8, 'MarkerFaceColor', '#7BC0CD', 'DisplayName', 'Pareto Optimal Solutions');
plot(f1_min, f2_1, 'p', 'MarkerSize', 16, 'Color', '#ED8D5A', 'MarkerFaceColor', '#ED8D5A', 'DisplayName', 'Minimum Cost Point');
plot(f1_2, f2_min, 's', 'MarkerSize', 14, 'Color', '#51999F', 'MarkerFaceColor', '#51999F', 'DisplayName', 'Minimum Emissions Point');
xlabel('Total Operating Cost ($)', 'FontSize', 12); 
ylabel('Total Carbon Emissions (kg CO2)', 'FontSize', 12);
legend('show', 'Location', 'best');
grid on;
