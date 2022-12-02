within DynEq.Elements.Solar.WECC;

model GenerationTripping
  parameter Real Lv0;
  parameter Real Lv1;
  parameter Real Lv2;
  parameter Real Lv3;
  Modelica.Blocks.Interfaces.RealInput u annotation(
    Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput TrpLow annotation(
    Placement(visible = true, transformation(origin = {110, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput TrpHigh annotation(
    Placement(visible = true, transformation(origin = {110, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.Solar.LiegeModel.Minimum compute_minimum(y_start=1) annotation(
    Placement(visible = true, transformation(origin = {70, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput TrpLowIntern annotation(
    Placement(visible = true, transformation(origin = {-10, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {120, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.Limiter limiter_low(uMax = 1, uMin = 0)  annotation(
    Placement(visible = true, transformation(origin = {30, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
// Tripping for low values of the input signal
  TrpLowIntern = (u - Lv0) / (Lv1 - Lv0);
  TrpHigh =1;
  connect(compute_minimum.y, TrpLow) annotation(
    Line(points = {{82, 50}, {110, 50}}, color = {0, 0, 127}));
  connect(limiter_low.u, TrpLowIntern) annotation(
    Line(points = {{18, 50}, {-10, 50}}, color = {0, 0, 127}));
  connect(limiter_low.y, compute_minimum.u) annotation(
    Line(points = {{42, 50}, {58, 50}}, color = {0, 0, 127}));
  annotation(
    Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {32, 93}, extent = {{-132, 7}, {68, -13}}, textString = "Generation Tripping"), Text(origin = {144, 53}, extent = {{-132, 7}, {-44, -13}}, textString = "TripLow"), Text(origin = {144, -47}, extent = {{-132, 7}, {-44, -13}}, textString = "TripHigh"), Text(origin = {36, 3}, extent = {{-132, 7}, {-44, -13}}, textString = "Input")}));
end GenerationTripping;
