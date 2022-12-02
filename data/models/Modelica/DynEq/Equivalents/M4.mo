within DynEq.Equivalents;

model M4
  extends Modelica.Icons.Example;
  // General Parameters
  parameter Real P_0 = 17.593473230782624e3;
  parameter Real Q_0 = 6.44994877688119e3;
  parameter Real v_0 = 1.0003343*20e3;
  parameter Real angle_0 = 0;
  parameter Real omega_0 = 0.9997503;
   // Tunable parameters
  parameter Real KPLP = 0.24549457830753876;
  parameter Real KPLQ = 1;
  parameter Real Mpv = 0.5127903688798994;
  // Components
  inner DynEq.Essentials.SystemBase SysData(SBase = 1e+6, fBase = 50) annotation(
    Placement(visible = true, transformation(origin = {-70, 90}, extent = {{-30, -10}, {30, 10}}, rotation = 0)));
  DynEq.Elements.ElmVac  source(P_0 = -P_0, Q_0=-Q_0 ,VBase(displayUnit = "V") = 20e3,angle_0 = angle_0, v_0 = v_0) annotation(
    Placement(visible = true, transformation(origin = {-48, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.Solar.WECC.PlantPVD1 solar( Ddn = 0, Ft0 = 0, Ft1 = 0.0001, Ft2 = 200, Ft3 = 201, Imax = 1.1487749575056005, M_b = -solar.P_0 / Mpv , P_0(fixed = false), PqFlag = false, Q_0(fixed = false), Qmn = -100, Qmx = 100, Tg = 0.029999999959649696, VBase(displayUnit = "V") = 20e3, Vt0 = 0.02518511722573208, Vt1 = 0.5209725867935494, Vt2 = 200, Vt3 = 201, Xc = 0, angle_0(fixed = false), dqdv = 4.062342889612549, fr_recov = 0, v0 = 0.7912697529915187, v1 = 1.2269022561147775, v_0(fixed = false), vr_recov = 0) annotation(
    Placement(visible = true, transformation(origin = {40, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  DynEq.Elements.CompositeLoad composite_load( Kpm=0.14378444928586107, Mlf=0.1, P_0(fixed = false), Q_0(fixed = false), VBase(displayUnit = "V") = 20e3, angle_0(fixed = false), omega_0 = omega_0, v_0(fixed = false)) annotation(
    Placement(visible = true, transformation(origin = {40, -30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput u2 annotation(
    Placement(visible = true, transformation(origin = {-80, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(extent = {{-120, -80}, {-80, -40}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput u1 annotation(
    Placement(visible = true, transformation(origin = {-80, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(extent = {{-120, 40}, {-80, 80}}, rotation = 0)));

  final parameter Real g = line.Z.re/(line.Z.re^2+line.Z.im^2);
  final parameter Real b = -line.Z.im/(line.Z.re^2+line.Z.im^2);
  
  final parameter Types.ActivePower P1(fixed=false);
  final parameter Types.ReactivePower Q1(fixed=false);
  final parameter Types.Voltage V1(fixed=false, start=20e3);
  final parameter Types.Angle d1(fixed=false, start=0);

  DynEq.Elements.Line line(PStartA=P_0, PStartB(fixed = false), QStartA=Q_0, QStartB(fixed = false), R = 10, SNomA(displayUnit = "V.A") = 10e3, SNomB(displayUnit = "V.A") = 10e3, VBaseA(displayUnit = "V") = 20e3, VBaseB(displayUnit = "V") = 20e3, VStartA = v_0, VStartB(fixed = false), X = 10, angleStartA = angle_0, angleStartB(fixed = false)) annotation(
    Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

initial equation


// Production values
  solar.P_0 = - (KPLP-1) * P1;
  solar.Q_0 = - (KPLQ-1) * Q1;
// Consumption values
  composite_load.P_0 = P1 * KPLP;
  composite_load.Q_0 = Q1 * KPLQ;

// Trafo secondary - composite_load
  P_0 = g * v_0 ^ 2 - v_0 * V1 * (g * cos(angle_0 - d1) + b * sin(angle_0 - d1));
  Q_0 = -b*v_0^2 - v_0*V1*(g*sin(angle_0-d1) - b*cos(angle_0-d1));
  -P1 = g*V1^2 - v_0*V1*(g*cos(d1-angle_0) + b*sin(d1-angle_0));
  -Q1 = -b*V1^2 - v_0*V1*(g*sin(d1-angle_0) - b*cos(d1-angle_0));


//
  line.angleStartB = d1;
  line.VStartB = V1;
  line.PStartB = -P1;
  line.QStartB = -Q1;
//
  composite_load.v_0 = V1;
  composite_load.angle_0 = d1;
  solar.v_0 = V1;
  solar.angle_0 = d1;
  
  equation
  connect(u2, source.fPu) annotation(
    Line(points = {{-80, 10}, {-66, 10}, {-66, 6}, {-58, 6}}, color = {0, 0, 127}));
  connect(u1, source.vPu) annotation(
    Line(points = {{-80, -10}, {-66, -10}, {-66, -4}, {-58, -4}}, color = {0, 0, 127}));
  connect(solar.omega, u2) annotation(
    Line(points = {{50, 6}, {58, 6}, {58, 20}, {-74, 20}, {-74, 10}, {-80, 10}}, color = {0, 0, 127}));
  connect(composite_load.win, u2) annotation(
    Line(points = {{50, -30}, {58, -30}, {58, 20}, {-80, 20}, {-80, 10}}, color = {0, 0, 127}));
  connect(source.p, line.termA) annotation(
    Line(points = {{-48, 0}, {-11, 0}}, color = {0, 85, 0}));
  connect(line.termB, composite_load.p) annotation(
    Line(points = {{11, 0}, {26, 0}, {26, -30}, {30, -30}}, color = {0, 85, 0}));
  connect(solar.p, line.termB) annotation(
    Line(points = {{40, 0}, {10, 0}}, color = {0, 85, 0}));
  annotation(
    __OpenModelica_simulationFlags(csvInput = "C:/Users/trabu/dev/DynEq/data/experiments/exp1/inputs_om.csv", lv = "LOG_STATS", s = "rungekutta"),
    experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-6, Interval = 0.01));
end M4;