within DynEq.Elements.Solar.LiegeModel;

model PVSmall
  extends DynEq.Elements.BaseClasses.OnePort(p(i(re(nominal=M_b/VBase),im(nominal=M_b/VBase))));
  parameter SI.ApparentPower M_b "Base power of the PV generation" annotation(
    Dialog(group = "Power flow data"));
  parameter SI.PerUnit Imax = 1.1 "Maximum allowable total converter current" annotation(
    Dialog(group = "PVD1 Model Parameters"));
  parameter Boolean PqFlag "Priority on current limit flag: 1=P prio.; 0 = Q prio." annotation(
    Dialog(group = "PVD1 Model Parameters"));
  parameter SI.Time Tg = 0.02 "Inverter current regulator time constat" annotation(
    Dialog(group = "PVD1 Model Parameters"));
  parameter Real Vmin = 0.2 annotation(
    Dialog(group = "Partial Tripping Parameters"));
  parameter Real c = 0.3 annotation(
    Dialog(group = "Partial Tripping Parameters"));
  parameter Real d = 0.2 annotation(
    Dialog(group = "Partial Tripping Parameters"));
  DynEq.Elements.Solar.ElmGenstat static_generator(M_b = M_b, P_0 = P_0, Q_0 = Q_0, VBase = VBase, angle_0 = angle_0, v_0 = v_0) annotation(
    Placement(visible = true, transformation(origin = {50, -3.55271e-15}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput omega annotation(
    Placement(visible = true, transformation(origin = {100, 70}, extent = {{20, -20}, {-20, 20}}, rotation = 0), iconTransformation(origin = {-100, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  DynEq.Elements.Solar.LiegeModel.InverterControlSmall control_small(Imax = Imax, PqFlag = PqFlag, Pref = -P_0 / M_b, Qref = -Q_0 / M_b, Tg = Tg, Vmin = Vmin, c = c, d = d, u_0 = v_0 / VBase)  annotation(
    Placement(visible = true, transformation(origin = {10, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(static_generator.p, p) annotation(
    Line(points = {{60, 0}, {110, 0}}, color = {0, 0, 255}));
  connect(omega, static_generator.omega) annotation(
    Line(points = {{100, 70}, {50, 70}, {50, 10}, {50, 10}}, color = {0, 0, 127}));
  connect(control_small.Ip, static_generator.id_ref) annotation(
    Line(points = {{22, 6}, {40, 6}}, color = {0, 0, 127}));
  connect(control_small.Iq, static_generator.iq_ref) annotation(
    Line(points = {{22, -6}, {40, -6}}, color = {0, 0, 127}));
  connect(static_generator.v, control_small.Vt) annotation(
    Line(points = {{62, -6}, {66, -6}, {66, -20}, {-6, -20}, {-6, 0}, {110, 0}}, color = {0, 0, 127}));
protected
  annotation(
    Icon(graphics = {Ellipse(origin = {50, -48}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-150, 148}, {50, -52}}, endAngle = 360), Line(origin = {-0.0673837, 60.7896}, points = {{0, -83}, {0, 25}}, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 10), Text(origin = {-1, -52}, extent = {{-65, 12}, {65, -12}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)),
    __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
end PVSmall;
