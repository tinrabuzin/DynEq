within DynEq.Elements.Solar.LiegeModel;

model PVLarge
  extends DynEq.Elements.BaseClasses.OnePort(p(i(re(nominal=M_b/VBase),im(nominal=M_b/VBase))));
  parameter SI.ApparentPower M_b "Base power of the PV generation" annotation(
    Dialog(group = "Power flow data"));
  parameter SI.PerUnit Imax = 1.1 "Maximum allowable total converter current" annotation(
    Dialog(group = "PVD1 Model Parameters"));
  parameter Boolean PqFlag "Priority on current limit flag: 1=P prio.; 0 = Q prio." annotation(
    Dialog(group = "PVD1 Model Parameters"));
  parameter SI.Time Tg = 0.02 "Inverter current regulator time constat" annotation(
    Dialog(group = "PVD1 Model Parameters"));
  parameter SI.PerUnit Xc = 0 "Line drop compensation reactance" annotation(
    Dialog(group = "PVD1 Model Parameters"));
  parameter SI.PerUnit Qmx = 0.328 "Maximum reactive power" annotation(
    Dialog(group = "PVD1 Model Parameters"));
  parameter SI.PerUnit Qmn = -0.328 "Minimum reactive power" annotation(
    Dialog(group = "PVD1 Model Parameters"));
  parameter SI.PerUnit v0 = 0.9 "Low voltage threshold for Volt/Var Control" annotation(
    Dialog(group = "PVD1 Model Parameters"));
  parameter SI.PerUnit v1 = 1.1 "High voltage threshold for Volt/Var Control" annotation(
    Dialog(group = "PVD1 Model Parameters"));
  parameter SI.PerUnit dqdv = 0 "Voltage/Var droop compensation" annotation(
    Dialog(group = "PVD1 Model Parameters"));
  parameter Real Vr = 0.9 annotation(
    Dialog(group = "Partial Tripping Parameters"));
  parameter Real Vmin = 0.2 annotation(
    Dialog(group = "Partial Tripping Parameters"));
  parameter Real c = 0.3 annotation(
    Dialog(group = "Partial Tripping Parameters"));
  parameter Real d = 0.2 annotation(
    Dialog(group = "Partial Tripping Parameters"));
  parameter Real u = 1 annotation(
    Dialog(group = "Partial Tripping Parameters"));
  parameter Real T1 = 0.3 annotation(
    Dialog(group = "Partial Tripping Parameters"));
  parameter Real T2 = 0.7 annotation(
    Dialog(group = "Partial Tripping Parameters"));
  DynEq.Elements.Solar.ElmGenstat static_generator(M_b = M_b, P_0 = P_0, Q_0 = Q_0, VBase = VBase, angle_0 = angle_0, v_0 = v_0) annotation(
    Placement(visible = true, transformation(origin = {50, -3.55271e-15}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput omega annotation(
    Placement(visible = true, transformation(origin = {100, 70}, extent = {{20, -20}, {-20, 20}}, rotation = 0), iconTransformation(origin = {-100, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  DynEq.Elements.Solar.LiegeModel.InverterControlLarge control_large(Imax = Imax, PqFlag = PqFlag, Pref = -P_0 / M_b, Qmn = Qmn, Qmx = Qmx, Qref = -Q_0 / M_b, T1 = T1, T2 = T2, Tg = Tg, Vmin = Vmin, Vr = Vr, Xc = Xc, c = c, d = d, dqdv = dqdv, u = u, u_0 = v_0 / VBase, v0 = v0, v1 = v1)  annotation(
    Placement(visible = true, transformation(origin = {10, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(static_generator.p, p) annotation(
    Line(points = {{60, 0}, {110, 0}}, color = {0, 0, 255}));
  connect(omega, static_generator.omega) annotation(
    Line(points = {{100, 70}, {50, 70}, {50, 10}, {50, 10}}, color = {0, 0, 127}));
  connect(control_large.Ip, static_generator.id_ref) annotation(
    Line(points = {{20, 6}, {40, 6}, {40, 4}, {40, 4}}, color = {0, 0, 127}));
  connect(control_large.Iq, static_generator.iq_ref) annotation(
    Line(points = {{20, -4}, {40, -4}, {40, -6}, {40, -6}}, color = {0, 0, 127}));
  connect(static_generator.i, control_large.It) annotation(
    Line(points = {{62, 6}, {64, 6}, {64, 14}, {-4, 14}, {-4, 4}, {0, 4}, {0, 4}}, color = {0, 0, 127}));
  connect(static_generator.v, control_large.Vt) annotation(
    Line(points = {{62, -8}, {66, -8}, {66, 16}, {-6, 16}, {-6, 8}, {0, 8}, {0, 8}}, color = {0, 0, 127}));
protected
  annotation(
    Icon(graphics = {Ellipse(origin = {50, -48}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-150, 148}, {50, -52}}, endAngle = 360), Line(origin = {-0.0673837, 60.7896}, points = {{0, -83}, {0, 25}}, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 10), Text(origin = {-1, -52}, extent = {{-65, 12}, {65, -12}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)),
    __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
end PVLarge;
