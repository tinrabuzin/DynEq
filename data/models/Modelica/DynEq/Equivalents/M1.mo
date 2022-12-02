within DynEq.Equivalents;

model M1
  extends Modelica.Icons.Example;
  // General Parameters
  parameter Real P_0 = 26.446944539574236e3;
  parameter Real Q_0 = 11.741284438756121e3;
  parameter Real v_0 = 1.0003343 * 20e3;
  parameter Real angle_0 = 1;
  parameter Real omega_0 = 0.9997503;
  // Variables and parameters for powerflow calculations


  // Components
  inner DynEq.Essentials.SystemBase SysData(SBase = 1e+6, fBase = 50) annotation(
    Placement(visible = true, transformation(origin = {-64, 90}, extent = {{-30, -10}, {30, 10}}, rotation = 0)));
  DynEq.Elements.ElmVac source(P_0 = -P_0, Q_0 = -Q_0, VBase(displayUnit = "V") = 20e3, angle_0 = angle_0, v_0 = v_0) annotation(
    Placement(visible = true, transformation(origin = {-50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.CompositeLoad composite_load(VBase = 20e3, P_0(fixed=false), Q_0(fixed=false),  v_0(fixed=false), angle_0(fixed=false), omega_0 = omega_0)  annotation(
    Placement(visible = true, transformation(origin = {30, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput u1 annotation(
    Placement(visible = true, transformation(origin = {-90, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(extent = {{-120, 40}, {-80, 80}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput u2 annotation(
    Placement(visible = true, transformation(origin = {-90, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(extent = {{-120, -80}, {-80, -40}}, rotation = 0)));   
    Elements.Line line(PStartA = P_0, PStartB(fixed = false), QStartA = Q_0, QStartB(fixed = false), SNomA(displayUnit = "V.A") = 10e3, SNomB(displayUnit = "V.A") = 10e3, VBaseA(displayUnit = "V") = 20e3, VBaseB(displayUnit = "V") = 20e3, VStartA = v_0, VStartB(fixed = false), X = 36.81172403721348, angleStartA = angle_0, angleStartB(fixed = false)) annotation(
    Placement(visible = true, transformation(origin = {-10, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

// Parameters for powerflow calculations
  final parameter Real b = -line.Z.im/(line.Z.re^2+line.Z.im^2);
  
  final parameter Types.ActivePower P1(fixed=false);
  final parameter Types.ReactivePower Q1(fixed=false);
  final parameter Types.Voltage V1(fixed=false, start=20e3);
  final parameter Types.Angle d1(fixed=false, start=angle_0);
initial equation
// Consumption values
  P_0 =- v_0 * V1 *  b * sin(angle_0 - d1);
  Q_0 = -b * v_0^2   + v_0 * V1 *   b * cos(angle_0 - d1);
  -P1 = - v_0 * V1 * b * sin(d1 - angle_0);
  -Q1 = -b * V1^2    + v_0 * V1 *  b * cos(d1 - angle_0);
  
  composite_load.P_0 = P1;
  composite_load.Q_0 = Q1;
  composite_load.v_0 = V1;
  composite_load.angle_0 = d1;
//
  line.angleStartB = d1;
  line.VStartB = V1;
  line.PStartB = -P1;
  line.QStartB = -Q1;
equation
  connect(composite_load.win, u2) annotation(
    Line(points = {{40, 0}, {50, 0}, {50, 30}, {-90, 30}}, color = {0, 0, 127}));
  connect(u2, source.fPu) annotation(
    Line(points = {{-90, 30}, {-70, 30}, {-70, 6}, {-60, 6}}, color = {0, 0, 127}));
  connect(u1, source.vPu) annotation(
    Line(points = {{-90, -30}, {-70, -30}, {-70, -4}, {-60, -4}}, color = {0, 0, 127}));
  connect(line.termB, composite_load.p) annotation(
    Line(points = {{0, 0}, {30, 0}}, color = {0, 85, 0}));
  connect(source.p, line.termA) annotation(
    Line(points = {{-50, 0}, {-20, 0}}, color = {0, 85, 0}));
  annotation(
    __OpenModelica_simulationFlags(csvInput = "C:/Users/trabu/dev/DynEq/data/experiments/exp1/inputs_om.csv", lv = "LOG_STATS", s = "rungekutta"),
    experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-6, Interval = 0.01),
  __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian --replaceHomotopy=actual");
end M1;