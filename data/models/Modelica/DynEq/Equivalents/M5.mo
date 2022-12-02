within DynEq.Equivalents;

model M5
  extends Modelica.Icons.Example;
  // General Parameters
  parameter Real P_0 = 17.48294575617976e3 ;
  parameter Real Q_0 = 6.3518083493054105e3;
  parameter Real v_0 = 1.0003343 * 20e3;
  parameter Real angle_0 = 1;
  parameter Real omega_0 = 0.9997503;
  parameter Real Mpv = 0.6908858660550528;
  parameter Real MPV = 0.1389568509275564;
  parameter Real K_P=1.524174276070246;
  parameter Real K_Q=1;
  parameter Real FP_PV=0.3256234010565161;
  parameter Real FQ_PV=0.8854910129881104;

  // Components
  inner DynEq.Essentials.SystemBase SysData(SBase (displayUnit = "V.A") = 1e4, fBase = 50) annotation(
    Placement(visible = true, transformation(origin = {-70, 90}, extent = {{-30, -10}, {30, 10}}, rotation = 0)));
  DynEq.Elements.ElmVac source(P_0 = -P_0, Q_0 = -Q_0, VBase(displayUnit = "V") = 20e3, angle_0 = angle_0, v_0 = v_0) annotation(
    Placement(visible = true, transformation(origin = {-48, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.Line line(R=0.5687412799846374, X=36.81172403721348,
  PStartA=P_0,  PStartB(fixed=false), QStartA=Q_0, QStartB(fixed=false), SNomA(displayUnit = "V.A") = 10e3, SNomB(displayUnit = "V.A") = 10e3, VBaseA (displayUnit = "V") = 20e3, VBaseB (displayUnit = "V") = 20e3,
  VStartA=v_0,  VStartB(fixed=false), angleStartA=angle_0, angleStartB(fixed=false)) annotation(
    Placement(visible = true, transformation(origin = {-10, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput u1 annotation(
    Placement(visible = true, transformation(origin = {-80, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(extent = {{-120, 40}, {-80, 80}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput u2 annotation(
    Placement(visible = true, transformation(origin = {-80, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(extent = {{-120, -80}, {-80, -40}}, rotation = 0)));
  DynEq.Elements.Solar.WECC.PlantPVD1 PVSmall(M_b = -PVSmall.P_0 / Mpv, P_0(fixed=false), PqFlag = false, Q_0(fixed=false), VBase = 20e3, angle_0(fixed = false), v_0(fixed = false), Qmx=100, Qmn=-100) annotation(
    Placement(visible = true, transformation(origin = {40, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  DynEq.Elements.Solar.WECC.PlantPVD1 PVLarge(M_b = -PVLarge.P_0 / MPV, P_0(fixed=false), PqFlag = false, Q_0(fixed=false), Qmn=-100, Qmx=100, VBase (displayUnit = "V") = 20e3, angle_0(fixed = false), v1 = 0.9999, v_0(fixed = false)) annotation(
    Placement(visible = true, transformation(origin = {40, -30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  DynEq.Elements.CompositeLoad composite_load(
    Rs=0.022385714633135736,
    Xs=0.04705445502458098,
    Xm=4.9907315720133925,
    Rr=0.0057682147961691,
    Xr=0.1796502249194955,
    H=1.184018860699618,
    expT=2.0093181296438734,
    expP=1.340374676870228,
    expQ=1.020485789709415,
    Kpm=0.20869547329584562, Mlf = 0.95, P_0(fixed=false), Q_0(fixed=false), angle_0(fixed = false), v_0(fixed = false), VBase = 20e3, omega_0 = omega_0)  annotation(
    Placement(visible = true, transformation(origin = {40, -60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  // Line parameters
  //final parameter Real g = line.Z.re/(line.Z.re^2+line.Z.im^2);
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
  
  PVLarge.P_0 = -FP_PV * (K_P - 1) * P1;
  PVLarge.Q_0 = -FQ_PV * (K_Q - 1)  * Q1;
  PVSmall.P_0 = -(1 - FP_PV) * (K_P - 1) * P1;
  PVSmall.Q_0 = -(1 - FQ_PV) * (K_Q - 1) * Q1;
  composite_load.P_0 = P1 * K_P;
  composite_load.Q_0 = Q1 * K_Q;
//
  line.angleStartB = d1;
  line.VStartB = V1;
  line.PStartB = -P1;
  line.QStartB = -Q1;
//
  composite_load.v_0 = V1;
  composite_load.angle_0 = d1;
  PVLarge.v_0 = V1;
  PVLarge.angle_0 = d1;
  PVSmall.v_0 = V1;
  PVSmall.angle_0 = d1;
  
equation
  connect(u2, source.fPu) annotation(
    Line(points = {{-80, 10}, {-66, 10}, {-66, 6}, {-58, 6}}, color = {0, 0, 127}));
  connect(u1, source.vPu) annotation(
    Line(points = {{-80, -10}, {-66, -10}, {-66, -4}, {-58, -4}}, color = {0, 0, 127}));
  connect(line.termB, PVSmall.p) annotation(
    Line(points = {{1, 0}, {40, 0}}, color = {0, 85, 0}));
  connect(PVSmall.omega, u2) annotation(
    Line(points = {{50, 6}, {60, 6}, {60, 20}, {-80, 20}, {-80, 10}}, color = {0, 0, 127}));
  connect(PVLarge.omega, u2) annotation(
    Line(points = {{50, -24}, {60, -24}, {60, 20}, {-80, 20}, {-80, 10}}, color = {0, 0, 127}));
  connect(composite_load.win, u2) annotation(
    Line(points = {{50, -60}, {60, -60}, {60, 20}, {-80, 20}, {-80, 10}}, color = {0, 0, 127}));
  connect(line.termA, source.p) annotation(
    Line(points = {{-21, 0}, {-48, 0}}, color = {0, 85, 0}));
  connect(line.termB, PVLarge.p) annotation(
    Line(points = {{1, 0}, {1, -30}, {40, -30}}, color = {0, 85, 0}));
  connect(line.termB, composite_load.p) annotation(
    Line(points = {{1, 0}, {1, -60}, {40, -60}}, color = {0, 85, 0}));
  annotation(
    __OpenModelica_simulationFlags(csvInput = "C:/Users/trabu/dev/DynEq/data/experiments/exp1/inputs_om.csv", lv = "LOG_STATS", s = "rungekutta"),
    experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-6, Interval = 0.01),
  __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian --replaceHomotopy=actual",
  Diagram(coordinateSystem(extent = {{-100, 100}, {60, -80}})));
end M5;