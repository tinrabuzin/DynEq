within DynEq.Equivalents;

model M2
  extends Modelica.Icons.Example;
  // General Parameters
  parameter Real P_0 = 18.549341834292655e3 ;
  parameter Real Q_0 = 6.258853e3;
  parameter Real v_0 = 1.0003343*20e3;
  parameter Real angle_0 = 0;
  parameter Real omega_0 = 0.9997503;
   // Tunable parameters
  parameter Real KPLP = 0.2;
  parameter Real KPLQ = 0.2;
  parameter Real Mpv = 0.01;
  // Components
  inner DynEq.Essentials.SystemBase SysData(SBase = 1e+6, fBase = 50) annotation(
    Placement(visible = true, transformation(origin = {-70, 90}, extent = {{-30, -10}, {30, 10}}, rotation = 0)));
  DynEq.Elements.ElmVac  source(P_0 = -P_0, Q_0=-Q_0 ,VBase(displayUnit = "V") = 20e3,angle_0 = angle_0, v_0 = v_0) annotation(
    Placement(visible = true, transformation(origin = {-48, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.Solar.WECC.PlantPVD1 solar( Ddn = 0, Ft0 = 0, Ft1 = 0.0001, Ft2 = 200, Ft3 = 201, M_b = -solar.P_0 / Mpv , P_0(fixed = false), PqFlag = false, Q_0(fixed = false), Qmn = -100, Qmx = 100, VBase(displayUnit = "V") = 20e3, Vt0 = 0.1, Vt1 = 0.95, Vt2 = 200, Vt3 = 201, Xc = 0, angle_0(fixed = false), fr_recov = 0, v1 = 1 + 1 - solar.v0, v_0(fixed = false), vr_recov = 0) annotation(
    Placement(visible = true, transformation(origin = {40, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  DynEq.Elements.CompositeLoad composite_load( Kpm = 0.2, Mlf = 0.5, P_0(fixed = false), Q_0(fixed = false), VBase(displayUnit = "V") = 20e3, angle_0(fixed = false), omega_0 = omega_0, v_0(fixed = false)) annotation(
    Placement(visible = true, transformation(origin = {40, -30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  DynEq.Elements.Line ln_ld( 
  PStartB(fixed=false), QStartB(fixed=false),R = 1, SNomA(displayUnit = "V.A") = 10e3, SNomB(displayUnit = "V.A") = 10e3, VBaseA (displayUnit = "V") = 20e3, VBaseB(displayUnit = "V") =20e3, VStartB(fixed=false), X = 1, angleStartB(fixed=false),
  VStartA(fixed=false), angleStartA(fixed=false),
  PStartA(fixed=false), QStartA(fixed=false)) annotation(
    Placement(visible = true, transformation(origin = {12, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.Line trafo(
  PStartB(fixed=false), QStartB(fixed=false),R = 0.070547, SNomA(displayUnit = "V.A") = 10e3, SNomB(displayUnit = "V.A") = 10e3,  VBaseA (displayUnit = "V") = 20e3, VBaseB(displayUnit = "V") =20e3, VStartB(fixed=false), X = 0.000511, angleStartB(fixed=false)) annotation(
    Placement(visible = true, transformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.Line ln_prod(
  PStartB(fixed=false), QStartB(fixed=false),R = 1, SNomA(displayUnit = "V.A") = 10e3, SNomB(displayUnit = "V.A") = 10e3,  VBaseA (displayUnit = "V") = 20e3, VBaseB(displayUnit = "V") =20e3, VStartB(fixed=false), X = 1, angleStartB(fixed=false),
  VStartA(fixed=false), angleStartA(fixed=false),
  PStartA(fixed=false), QStartA(fixed=false))  annotation(
    Placement(visible = true, transformation(origin = {12, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput u2 annotation(
    Placement(visible = true, transformation(origin = {-80, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(extent = {{-120, -80}, {-80, -40}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput u1 annotation(
    Placement(visible = true, transformation(origin = {-80, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(extent = {{-120, 40}, {-80, 80}}, rotation = 0)));

  final parameter Real xt = trafo.Z.im;
  final parameter Real rt = trafo.Z.re;
  final parameter Real gt = rt/(rt^2+xt^2);
  final parameter Real bt = -xt/(rt^2+xt^2);
  final parameter Real gsol = ln_prod.Z.re/(ln_prod.Z.re^2+ln_prod.Z.im^2);
  final parameter Real bsol = -ln_prod.Z.im/(ln_prod.Z.re^2+ln_prod.Z.im^2);
  final parameter Real gld = ln_ld.Z.re/(ln_ld.Z.re^2+ln_ld.Z.im^2);
  final parameter Real bld = -ln_ld.Z.im/(ln_ld.Z.re^2+ln_ld.Z.im^2);
  
  final parameter Types.ActivePower P1(fixed=false);
  final parameter Types.ReactivePower Q1(fixed=false);
  final parameter Types.Voltage V1(fixed=false, start=20e3);
  final parameter Types.Angle d1(fixed=false, start=0);
  
  final parameter Types.ActivePower Psol1(fixed=false);
  final parameter Types.ReactivePower Qsol1(fixed=false);
  final parameter Types.ActivePower Pld1(fixed=false);
  final parameter Types.ReactivePower Qld1(fixed=false);
  
  final parameter Types.ActivePower Psol2(fixed=false);
  final parameter Types.ReactivePower Qsol2(fixed=false);
  final parameter Types.Voltage Vsol(fixed=false, start=20e3);
  final parameter Types.Angle dsol(fixed=false, start=0);
  final parameter Types.ActivePower Pld2(fixed=false);
  final parameter Types.ReactivePower Qld2(fixed=false);
  final parameter Types.Voltage Vld(fixed=false,start=20e3);
  final parameter Types.Angle dld(fixed=false,start=0);

initial equation
// Trafo primary - trafo secondary
  P_0 = gt * v_0 ^ 2 - v_0 * V1 * (gt * cos(angle_0 - d1) + bt * sin(angle_0 - d1));
  Q_0 = -bt*v_0^2 - v_0*V1*(gt*sin(angle_0-d1) - bt*cos(angle_0-d1));
  -P1 = gt*V1^2 - v_0*V1*(gt*cos(d1-angle_0) + bt*sin(d1-angle_0));
  -Q1 = -bt*V1^2 - v_0*V1*(gt*sin(d1-angle_0) - bt*cos(d1-angle_0));
// Production values
  -Psol1 = KPLP/(1-KPLP)*P1;
  -Qsol1 = KPLQ/(1-KPLQ)*Q1;
// Consumption values
  Pld1 =  P1/(1-KPLP);
  Qld1 = Q1/(1-KPLQ);
// Trafo secondary - composite_load
  Pld1 = gld * V1 ^ 2 - V1 * Vld * (gld * cos(d1 - dld) + bld * sin(d1 - dld));
  Qld1 = -bld*V1^2 - V1*Vld*(gld*sin(d1-dld) - bld*cos(d1-dld));
  -Pld2 = gld*Vld^2 - V1*Vld*(gld*cos(dld-d1) + bld*sin(dld-d1));
  -Qld2 = -bld*Vld^2 - V1*Vld*(gld*sin(dld-d1) - bld*cos(dld-d1));
// Trafo secondary - solar
  Psol1 =   gsol * V1^2   - V1 * Vsol * (gsol * cos(d1 - dsol) + bsol * sin(d1 - dsol));
  Qsol1 =  -bsol * V1^2   - V1 * Vsol * (gsol * sin(d1 - dsol) - bsol * cos(d1 - dsol));
  -Psol2 =  gsol * Vsol^2 - V1 * Vsol * (gsol * cos(dsol - d1) + bsol * sin(dsol - d1));
  -Qsol2 = -bsol * Vsol^2 - V1 * Vsol * (gsol * sin(dsol - d1) - bsol * cos(dsol - d1));
// Assigning the pf results
  //
  trafo.angleStartB = d1;
  trafo.VStartB = V1;
  trafo.PStartB = -P1;
  trafo.QStartB = -Q1;
  trafo.VStartA = v_0;
  trafo.angleStartA = angle_0;
  trafo.PStartA = P_0;
  trafo.QStartA = Q_0;
  //
  ln_ld.angleStartB = dld;
  ln_ld.VStartB = Vld;
  ln_ld.PStartB = -Pld2;
  ln_ld.QStartB = -Qld2;
  ln_ld.angleStartA = d1;
  ln_ld.VStartA = V1;
  ln_ld.PStartA = Pld1;
  ln_ld.QStartA = Qld1;
  //
  ln_prod.angleStartB = dsol;
  ln_prod.VStartB = Vsol;
  ln_prod.PStartB = Psol2;
  ln_prod.QStartB = Qsol2;
  ln_prod.angleStartA = d1;
  ln_prod.VStartA = V1;
  ln_prod.PStartA = Psol1;
  ln_prod.QStartA = Qsol1;
  //
  composite_load.P_0 = Pld2;
  composite_load.Q_0 = Qld2;
  composite_load.v_0 = Vld;
  composite_load.angle_0 = dld;
  solar.P_0 = Psol2;
  solar.Q_0 = Qsol2;
  solar.v_0 = Vsol;
  solar.angle_0 = dsol;
  equation
  connect(u2, source.fPu) annotation(
    Line(points = {{-80, 10}, {-66, 10}, {-66, 6}, {-58, 6}}, color = {0, 0, 127}));
  connect(u1, source.vPu) annotation(
    Line(points = {{-80, -10}, {-66, -10}, {-66, -4}, {-58, -4}}, color = {0, 0, 127}));
  connect(source.p, trafo.termA) annotation(
    Line(points = {{-36, 0}, {-30, 0}}, color = {0, 85, 0}));
  connect(trafo.termB, ln_prod.termA) annotation(
    Line(points = {{-10, 0}, {2, 0}}, color = {0, 85, 0}));
  connect(ln_ld.termA, trafo.termB) annotation(
    Line(points = {{2, -30}, {-4, -30}, {-4, 0}, {-10, 0}}, color = {0, 85, 0}));
  connect(solar.omega, u2) annotation(
    Line(points = {{50, 6}, {58, 6}, {58, 20}, {-74, 20}, {-74, 10}, {-80, 10}}, color = {0, 0, 127}));
  connect(composite_load.win, u2) annotation(
    Line(points = {{50, -30}, {58, -30}, {58, 20}, {-80, 20}, {-80, 10}}, color = {0, 0, 127}));
  connect(ln_prod.termB, solar.p) annotation(
    Line(points = {{22, 0}, {30, 0}}, color = {0, 85, 0}));
  connect(ln_ld.termB, composite_load.p) annotation(
    Line(points = {{22, -30}, {30, -30}}, color = {0, 85, 0}));
  annotation(
    __OpenModelica_simulationFlags(csvInput = "C:/Users/trabu/dev/DynEq/data/input_signals/EID021_om.csv", lv = "LOG_STATS", s = "rungekutta"),
    experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-6, Interval = 0.01),
  __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian --replaceHomotopy=actual");
end M2;