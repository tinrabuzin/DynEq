within DynEq.Elements.Solar.LiegeModel;

model TrippingSmall
  parameter Real Vmin = 0.2;
  parameter Real c = 0.3;
  parameter Real d = 0.2;
  Modelica.Blocks.Interfaces.RealOutput trip annotation(
    Placement(visible = true, transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput v annotation(
    Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  DynEq.Elements.Solar.LiegeModel.Minimum min_f1(y_start = 1)  annotation(
    Placement(visible = true, transformation(origin = {50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.Limiter limiter_low(uMax = 1, uMin = 0) annotation(
    Placement(visible = true, transformation(origin = {16, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
protected
  Modelica.Blocks.Interfaces.RealOutput f1 annotation(
    Placement(visible = true, transformation(origin = {-10, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
// Stage 1
  f1 = if v < Vmin then c / (1 - d) / Vmin * v - c * d / (1 - d) else 1;
  connect(min_f1.y, trip) annotation(
    Line(points = {{62, 0}, {110, 0}}, color = {0, 0, 127}));
  connect(f1, limiter_low.u) annotation(
    Line(points = {{-10, 0}, {4, 0}}, color = {0, 0, 127}));
  connect(limiter_low.y, min_f1.u) annotation(
    Line(points = {{28, 0}, {38, 0}}, color = {0, 0, 127}));
  annotation(
    __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
    Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {-8, 46}, extent = {{-68, 18}, {68, -18}}, textString = "Partial Tripping")}, coordinateSystem(initialScale = 0.1)));
end TrippingSmall;
