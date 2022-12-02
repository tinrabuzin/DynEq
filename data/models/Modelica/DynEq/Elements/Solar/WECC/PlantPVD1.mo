within DynEq.Elements.Solar.WECC;

model PlantPVD1
  extends DynEq.Elements.BaseClasses.OnePort(p(i(re(nominal=M_b/VBase),im(nominal=M_b/VBase))));
  parameter Types.ApparentPower M_b "Base power of the PV generation" annotation (Dialog(group="Power flow data"));
  parameter SI.PerUnit Imax = 1.1 "Maximum allowable total converter current"annotation(Dialog(group = "PVD1 Model Parameters"));
  parameter Boolean PqFlag "Priority on current limit flag: 1=P prio.; 0 = Q prio."annotation(Dialog(group = "PVD1 Model Parameters"));
  parameter SI.Time Tg = 0.02 "Inverter current regulator time constat"annotation(Dialog(group = "PVD1 Model Parameters"));
  parameter SI.PerUnit Xc = 0 "Line drop compensation reactance"annotation(Dialog(group = "PVD1 Model Parameters"));
  parameter SI.PerUnit Qmx = 0.328 "Maximum reactive power"annotation(Dialog(group = "PVD1 Model Parameters"));
  parameter SI.PerUnit Qmn = -0.328 "Minimum reactive power"annotation(Dialog(group = "PVD1 Model Parameters"));
  parameter SI.PerUnit v0 = 0.9 "Low voltage threshold for Volt/Var Control"annotation(Dialog(group = "PVD1 Model Parameters"));
  parameter SI.PerUnit v1 = 1.1 "High voltage threshold for Volt/Var Control"annotation(Dialog(group = "PVD1 Model Parameters"));
  parameter SI.PerUnit dqdv = 0 "Voltage/Var droop compensation"annotation(Dialog(group = "PVD1 Model Parameters"));
  parameter SI.PerUnit fdbd = -99 "Frequency deadband over frequency response" annotation(Dialog(group = "PVD1 Model Parameters"));
  parameter Real Ddn = 0 "Down regulation droop" annotation(Dialog(group = "PVD1 Model Parameters"));

  parameter Real vr_recov = 1 "Amount of generation to reconnect after voltage disconnection" annotation(Dialog(group = "PVD1 Model Parameters"));
  parameter Real fr_recov = 1 "Amount of generation to reconnect after frequency disconnection" annotation(Dialog(group = "PVD1 Model Parameters"));
  parameter Real Ft0 = 0.99 annotation(Dialog(group = "PVD1 Model Parameters"));
  parameter Real Ft1 = 0.995 annotation(Dialog(group = "PVD1 Model Parameters"));
  parameter Real Ft2 = 1.005 annotation(Dialog(group = "PVD1 Model Parameters"));
  parameter Real Ft3 = 1.01 annotation(Dialog(group = "PVD1 Model Parameters"));
  parameter Real Vt0 = 0.88 annotation(Dialog(group = "PVD1 Model Parameters"));
  parameter Real Vt1 = 0.9 annotation(Dialog(group = "PVD1 Model Parameters"));
  parameter Real Vt2 = 100 annotation(Dialog(group = "PVD1 Model Parameters"));
  parameter Real Vt3 = 200 annotation(Dialog(group = "PVD1 Model Parameters"));
  ElmGenstat static_generator(P_0 = P_0, Q_0 = Q_0, M_b = M_b, VBase = VBase, angle_0 = angle_0, v_0 = v_0)  annotation(
    Placement(visible = true, transformation(origin = {50, -3.55271e-15}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput omega annotation(
    Placement(visible = true, transformation(origin = {100, 60}, extent = {{10, -10}, {-10, 10}}, rotation = 0), iconTransformation(origin = {-100, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  DynEq.Elements.Solar.WECC.PVD1 pvd1( Ddn = Ddn, Ft0 = Ft0, Ft1 = Ft1, Ft2 = Ft2, Ft3 = Ft3, Imax = Imax, PqFlag = PqFlag,Pref = -P_0 / M_b, Qmn = Qmn, Qmx = Qmx, Qref = -Q_0 / M_b, Tg = Tg, Vt0 = Vt0, Vt1 = Vt1, Vt2 = Vt2, Vt3 = Vt3, Xc = Xc, dqdv = dqdv, fdbd = fdbd, fr_recov = fr_recov, u_0 = v_0 / VBase, v0 = v0, v1 = v1, vr_recov = vr_recov)  annotation(
    Placement(visible = true, transformation(origin = {10, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(static_generator.p, p) annotation(
    Line(points = {{60, 0}, {110, 0}}, color = {0, 0, 255}));
  connect(omega, static_generator.omega) annotation(
    Line(points = {{100, 60}, {50, 60}, {50, 10}}, color = {0, 0, 127}));
  connect(pvd1.Ip, static_generator.id_ref) annotation(
    Line(points = {{22, 6}, {40, 6}}, color = {0, 0, 127}));
  connect(pvd1.Iq, static_generator.iq_ref) annotation(
    Line(points = {{22, -6}, {40, -6}}, color = {0, 0, 127}));
  connect(static_generator.v, pvd1.Vt) annotation(
    Line(points = {{62, -6}, {76, -6}, {76, 20}, {-8, 20}, {-8, 6}, {2, 6}}, color = {0, 0, 127}));
  connect(static_generator.i, pvd1.It) annotation(
    Line(points = {{62, 6}, {68, 6}, {68, 14}, {-4, 14}, {-4, 0}, {2, 0}}, color = {0, 0, 127}));
  connect(omega, pvd1.freq) annotation(
    Line(points = {{100, 60}, {-12, 60}, {-12, -6}, {2, -6}}, color = {0, 0, 127}));
protected
  annotation (
    Icon(graphics = {Ellipse(origin = {50, -48}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-150, 148}, {50, -52}}, endAngle = 360), Line(origin = {-0.0673837, 60.7896}, points = {{0, -83}, {0, 25}}, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 10), Text(origin = {-1, -52}, extent = {{-65, 12}, {65, -12}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)),
    __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
end PlantPVD1;
