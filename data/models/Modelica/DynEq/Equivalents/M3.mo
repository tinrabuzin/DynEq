within DynEq.Equivalents;

model M3
  extends Modelica.Icons.Example;
  // General Parameters
  parameter Real P_0 = 18.549341834292655e3;
  parameter Real Q_0 = 6.258853e3;
  parameter Real v_0 = 1.0003343 * 20e3;
  parameter Real angle_0 = 0;
  parameter Real omega_0 = 0.9997503;
  // Tunable parameters
  parameter Real K_P = 2 "Needs to be > 1";
  parameter Real K_Q = 2 "Needs to be > 1";
  parameter Real FP_PV = 0.11 "Needs to be 0 <= F_P <= 1";
  parameter Real FQ_PV = 0.11 "Needs to be 0 <= F_Q <= 1";
  parameter Real Mpv = 0.1;
  parameter Real MPV = 0.1;
  // Components
  inner DynEq.Essentials.SystemBase SysData(SBase = 1e+6, fBase = 50) annotation(
    Placement(visible = true, transformation(origin = {-70, 90}, extent = {{-30, -10}, {30, 10}}, rotation = 0)));
  DynEq.Elements.ElmVac source(P_0 = -P_0, Q_0 = -Q_0, VBase = 20e3, angle_0 = angle_0, v_0 = v_0) annotation(
    Placement(visible = true, transformation(origin = {-48, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.Line ln_pvlarge(R = 0.01, X = 0.01, V_bA = 20e3, V_bB = 20e3) annotation(
    Placement(visible = true, transformation(origin = {12, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.Line trafo(R = 0.070547, X = 0.000511, V_bA = 20e3, V_bB = 20e3) annotation(
    Placement(visible = true, transformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.Line ln_pvsmall(R = 0.2, V_bA = 20e3, V_bB = 20e3, X = 0.1) annotation(
    Placement(visible = true, transformation(origin = {12, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.Line ln_ld(R = 0.01, X = 0.01, V_bA = 20e3, V_bB = 20e3) annotation(
    Placement(visible = true, transformation(origin = {12, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput u1 annotation(
    Placement(visible = true, transformation(origin = {-80, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(extent = {{-120, 40}, {-80, 80}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput u2 annotation(
    Placement(visible = true, transformation(origin = {-80, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(extent = {{-120, -80}, {-80, -40}}, rotation = 0)));
  DynEq.Elements.Solar.LiegeModel.PVSmall PVSmall(M_b = -PVSmall.P_0 / Mpv, P_0(fixed = false), PqFlag = false, Q_0(fixed = false), VBase (displayUnit = "V") = 20e3, angle_0(fixed = false), v_0(fixed = false)) annotation(
    Placement(visible = true, transformation(origin = {40, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  DynEq.Elements.Solar.LiegeModel.PVLarge PVLarge(M_b = -PVLarge.P_0 / MPV, P_0(fixed = false), PqFlag = false, Q_0(fixed = false), VBase (displayUnit = "V") = 20e3, Vmin = 0.8, Vr = 0.95, angle_0(fixed = false), v1 = 2 - PVLarge.v0, v_0(fixed = false), Qmx=100, Qmn=-100) annotation(
    Placement(visible = true, transformation(origin = {40, -30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  DynEq.Elements.CompositeLoad composite_load(P_0(fixed = false), Q_0(fixed = false), angle_0(fixed = false), v_0(fixed = false), VBase = 20e3, omega_0 = omega_0)  annotation(
    Placement(visible = true, transformation(origin = {40, -60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  // Line and trafo parameters
  final parameter Real xt = trafo.Z.im;
  final parameter Real rt = trafo.Z.re;
  final parameter Real gt = rt / (rt ^ 2 + xt ^ 2);
  final parameter Real bt = -xt / (rt ^ 2 + xt ^ 2);
  final parameter Real gPV = ln_pvlarge.Z.re / (ln_pvlarge.Z.re ^ 2 + ln_pvlarge.Z.im ^ 2);
  final parameter Real bPV = -ln_pvlarge.Z.im / (ln_pvlarge.Z.re ^ 2 + ln_pvlarge.Z.im ^ 2);
  final parameter Real gpv = ln_pvsmall.Z.re / (ln_pvsmall.Z.re ^ 2 + ln_pvsmall.Z.im ^ 2);
  final parameter Real bpv = -ln_pvsmall.Z.im / (ln_pvsmall.Z.re ^ 2 + ln_pvsmall.Z.im ^ 2);
  final parameter Real gld = ln_ld.Z.re / (ln_ld.Z.re ^ 2 + ln_ld.Z.im ^ 2);
  final parameter Real bld = -ln_ld.Z.im / (ln_ld.Z.re ^ 2 + ln_ld.Z.im ^ 2);
  // Algebraic variables at the trafo secondary
  final parameter Types.ActivePower Pprod(fixed = false);
  final parameter Types.ReactivePower Qprod(fixed = false);
  final parameter Types.ActivePower P1(fixed = false);
  final parameter Types.ReactivePower Q1(fixed = false);
  final parameter Types.Voltage V1(fixed = false, start = 20e3);
  final parameter Types.Angle d1(fixed = false, start = 0);
  // Large PV algebraic variables
  final parameter Types.ActivePower PPV1(fixed = false);
  final parameter Types.ReactivePower QPV1(fixed = false);
  final parameter Types.ActivePower PPV2(fixed = false);
  final parameter Types.ReactivePower QPV2(fixed = false);
  final parameter Types.Voltage VPV(fixed = false, start = 20e3);
  final parameter Types.Angle dPV(fixed = false, start = 0);
  // Small PV algebraic variables
  final parameter Types.ActivePower Ppv1(fixed = false);
  final parameter Types.ReactivePower Qpv1(fixed = false);
  final parameter Types.ActivePower Ppv2(fixed = false);
  final parameter Types.ReactivePower Qpv2(fixed = false);
  final parameter Types.Voltage Vpv(fixed = false, start = 20e3);
  final parameter Types.Angle dpv(fixed = false, start = 0);
  // Load algebraic variables
  final parameter Types.ActivePower Pld1(fixed = false);
  final parameter Types.ReactivePower Qld1(fixed = false);
  final parameter Types.ActivePower Pld2(fixed = false);
  final parameter Types.ReactivePower Qld2(fixed = false);
  final parameter Types.Voltage Vld(fixed = false, start = 20e3);
  final parameter Types.Angle dld(fixed = false, start = 0);
initial equation
/* =============== POWER FLOW EQUATIONS ============ */
// Trafo primary - trafo secondary
  P_0 = gt * v_0 ^ 2 - v_0 * V1 * (gt * cos(angle_0 - d1) + bt * sin(angle_0 - d1));
  Q_0 = (-bt * v_0 ^ 2) - v_0 * V1 * (gt * sin(angle_0 - d1) - bt * cos(angle_0 - d1));
  -P1 = gt * V1 ^ 2 - v_0 * V1 * (gt * cos(d1 - angle_0) + bt * sin(d1 - angle_0));
  -Q1 = (-bt * V1 ^ 2) - v_0 * V1 * (gt * sin(d1 - angle_0) - bt * cos(d1 - angle_0));
// Production values
  -Pprod = P1 / (K_P - 1);
  -Qprod = Q1 / (K_Q - 1);
  PPV1 = FP_PV * Pprod;
  QPV1 = FQ_PV * Qprod;
  Ppv1 = (1 - FP_PV) * Pprod;
  Qpv1 = (1 - FQ_PV) * Qprod;
// Consumption values
  Pld1 = K_P / (K_P - 1) * P1;
  Qld1 = K_Q / (K_Q - 1) * Q1;
// Trafo secondary - composite_load
  Pld1 = gld * V1 ^ 2 - V1 * Vld * (gld * cos(d1 - dld) + bld * sin(d1 - dld));
  Qld1 = (-bld * V1 ^ 2) - V1 * Vld * (gld * sin(d1 - dld) - bld * cos(d1 - dld));
  -Pld2 = gld * Vld ^ 2 - V1 * Vld * (gld * cos(dld - d1) + bld * sin(dld - d1));
  -Qld2 = (-bld * Vld ^ 2) - V1 * Vld * (gld * sin(dld - d1) - bld * cos(dld - d1));
// Trafo secondary - large PV
  PPV1 = gPV * V1 ^ 2 - V1 * VPV * (gPV * cos(d1 - dPV) + bPV * sin(d1 - dPV));
  QPV1 = (-bPV * V1 ^ 2) - V1 * VPV * (gPV * sin(d1 - dPV) - bPV * cos(d1 - dPV));
  -PPV2 = gPV * VPV ^ 2 - V1 * VPV * (gPV * cos(dPV - d1) + bPV * sin(dPV - d1));
  -QPV2 = (-bPV * VPV ^ 2) - V1 * VPV * (gPV * sin(dPV - d1) - bPV * cos(dPV - d1));
// Trafo secondary - small PV
  Ppv1 = gpv * V1 ^ 2 - V1 * Vpv * (gpv * cos(d1 - dpv) + bpv * sin(d1 - dpv));
  Qpv1 = (-bpv * V1 ^ 2) - V1 * Vpv * (gpv * sin(d1 - dpv) - bpv * cos(d1 - dpv));
  -Ppv2 = gpv * Vpv ^ 2 - V1 * Vpv * (gpv * cos(dpv - d1) + bpv * sin(dpv - d1));
  -Qpv2 = (-bpv * Vpv ^ 2) - V1 * Vpv * (gpv * sin(dpv - d1) - bpv * cos(dpv - d1));
/* =============== END OF POWER FLOW EQUATIONS ============ */
// Assigning the pf results
  composite_load.P_0 = Pld2;
  composite_load.Q_0 = Qld2;
  composite_load.v_0 = Vld;
  composite_load.angle_0 = dld;
  PVLarge.P_0 = PPV2;
  PVLarge.Q_0 = QPV2;
  PVLarge.v_0 = VPV;
  PVLarge.angle_0 = dPV;
  PVSmall.P_0 = Ppv2;
  PVSmall.Q_0 = Qpv2;
  PVSmall.v_0 = Vpv;
  PVSmall.angle_0 = dpv;
equation
  connect(u2, source.fPu) annotation(
    Line(points = {{-80, 10}, {-66, 10}, {-66, 6}, {-58, 6}}, color = {0, 0, 127}));
  connect(u1, source.vPu) annotation(
    Line(points = {{-80, -10}, {-66, -10}, {-66, -4}, {-58, -4}}, color = {0, 0, 127}));
  connect(trafo.termB, ln_pvsmall.termA) annotation(
    Line(points = {{-9, 0}, {2, 0}}, color = {0, 85, 0}));
  connect(ln_pvlarge.termA, trafo.termB) annotation(
    Line(points = {{2, -30}, {-2, -30}, {-2, 0}, {-9, 0}}, color = {0, 85, 0}));
  connect(ln_ld.termA, trafo.termB) annotation(
    Line(points = {{2, -60}, {-2, -60}, {-2, 0}, {-9, 0}}, color = {0, 85, 0}));
  connect(source.p, trafo.termA) annotation(
    Line(points = {{-48, 0}, {-31, 0}}, color = {0, 85, 0}));
  connect(ln_pvsmall.termB, PVSmall.p) annotation(
    Line(points = {{22, 0}, {40, 0}}, color = {0, 85, 0}));
  connect(PVSmall.omega, u2) annotation(
    Line(points = {{50, 6}, {60, 6}, {60, 20}, {-80, 20}, {-80, 10}}, color = {0, 0, 127}));
  connect(ln_pvlarge.termB, PVLarge.p) annotation(
    Line(points = {{22, -30}, {40, -30}}, color = {0, 85, 0}));
  connect(PVLarge.omega, u2) annotation(
    Line(points = {{50, -24}, {60, -24}, {60, 20}, {-80, 20}, {-80, 10}}, color = {0, 0, 127}));
  connect(ln_ld.termB, composite_load.p) annotation(
    Line(points = {{22, -60}, {40, -60}}, color = {0, 85, 0}));
  connect(composite_load.win, u2) annotation(
    Line(points = {{50, -60}, {60, -60}, {60, 20}, {-80, 20}, {-80, 10}}, color = {0, 0, 127}));
  annotation(
    __OpenModelica_simulationFlags(csvInput = "C:/Users/trabu/dev/DynEq/data/input_signals/EID021_om.csv", lv = "LOG_STATS", s = "rungekutta"),
    experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-6, Interval = 0.01),
  __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian --replaceHomotopy=actual");
end M3;
