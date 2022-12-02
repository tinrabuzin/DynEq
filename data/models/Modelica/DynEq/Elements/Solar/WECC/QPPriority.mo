within DynEq.Elements.Solar.WECC;
model QPPriority
  parameter Boolean PqFlag annotation(Evaluate = true);
  parameter SI.PerUnit Imax;
  Modelica.Blocks.Interfaces.RealInput Ip annotation (
    Placement(visible = true, transformation(origin = {-100, 50}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, 50}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput Iq annotation (
    Placement(visible = true, transformation(origin = {-100, -50}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, -50}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput Iqcmd annotation (
    Placement(visible = true, transformation(origin = {110, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput Ipcmd annotation (
    Placement(visible = true, transformation(origin = {110, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.VariableLimiter IpLimiter annotation (
    Placement(visible = true, transformation(origin = {0, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.VariableLimiter IqLimiter annotation (
    Placement(visible = true, transformation(origin = {0, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput Iqmin = -Iqmax;
  Modelica.Blocks.Interfaces.RealOutput Iqmax = homotopy(if PqFlag then sqrt(Imax^2 - Ipcmd^2) else Imax, Imax);
  Modelica.Blocks.Interfaces.RealOutput Ipmax = homotopy(if PqFlag then Imax else sqrt(Imax^2 - Iqcmd^2), Imax);
  Modelica.Blocks.Sources.Constant Ipmin(k = 0)  annotation(
    Placement(visible = true, transformation(origin = {-50, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation

  connect(IqLimiter.y, Iqcmd) annotation (
    Line(points = {{11, -50}, {110, -50}}, color = {0, 0, 127}));
  connect(IpLimiter.y, Ipcmd) annotation (
    Line(points = {{11, 50}, {110, 50}}, color = {0, 0, 127}));
  connect(Ip, IpLimiter.u) annotation (
    Line(points = {{-100, 50}, {-12, 50}}, color = {0, 0, 127}));
  connect(Iq, IqLimiter.u) annotation (
    Line(points = {{-100, -50}, {-12, -50}}, color = {0, 0, 127}));
  connect(Ipmax, IpLimiter.limit1);
  connect(Iqmin, IqLimiter.limit2);
  connect(Iqmax, IqLimiter.limit1);
  connect(Ipmin.y, IpLimiter.limit2) annotation(
    Line(points = {{-38, 10}, {-20, 10}, {-20, 42}, {-12, 42}}, color = {0, 0, 127}));
  annotation (
    Icon(graphics={  Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {20, 100}, extent = {{-120, 0}, {80, -20}}, textString = "PQ Priority"), Text(origin = {44, 60}, extent = {{-120, 0}, {-78, -20}}, textString = "Ipref"), Text(origin = {44, -40}, extent = {{-120, 0}, {-78, -20}}, textString = "Iqref"), Text(origin = {172, 52}, extent = {{-120, 0}, {-78, -20}}, textString = "Ipcmd"), Text(origin = {172, -28}, extent = {{-120, 0}, {-78, -20}}, textString = "Iqcmd")}));
end QPPriority;
