within DynEq.Elements.Solar.LiegeModel;

model TrippingLarge
  parameter Real Vr = 0.9;
  parameter Real Vmin = 0.2;
  parameter Real c = 0.3;
  parameter Real d = 0.2;
  parameter Real u = 1;
  parameter Real T1 = 0.3;
  parameter Real T2 = 3;
  Modelica.Blocks.Interfaces.RealOutput trip annotation(
    Placement(visible = true, transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput v annotation(
    Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  DynEq.Elements.Solar.LiegeModel.Minimum min_f1(y_start = 1)  annotation(
    Placement(visible = true, transformation(origin = {-10, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Product prod_f1_f2 annotation(
    Placement(visible = true, transformation(origin = {30, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Logical.Timer timer annotation(
    Placement(visible = true, transformation(origin = {30, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Logical.LessEqual lessEqual annotation(
    Placement(visible = true, transformation(origin = {-10, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant limit(k = Vr) annotation(
    Placement(visible = true, transformation(origin = {-50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.Solar.LiegeModel.Minimum min_f3(y_start = 1)  annotation(
    Placement(visible = true, transformation(origin = {-10, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Product prod_f1_f2_f3 annotation(
    Placement(visible = true, transformation(origin = {70, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const(k = 1) annotation(
    Placement(visible = true, transformation(origin = {-80, -66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.Limiter limiter_low(uMax = 1, uMin = 0) annotation(
    Placement(visible = true, transformation(origin = {-46, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.Limiter limiter(uMax = 1, uMin = 0) annotation(
    Placement(visible = true, transformation(origin = {-46, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
protected
  Modelica.Blocks.Interfaces.RealOutput f3 annotation(
    Placement(visible = true, transformation(origin = {-70, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-40, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput f1 annotation(
    Placement(visible = true, transformation(origin = {-70, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-40, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  if (timer.y <= T1 and timer.y > 0) then
    f1 = if v < Vmin then c / (1 - d) / Vmin * v - c * d / (1 - d) else 1;
    f3 = 1;
  elseif (timer.y >= T1) then
      f3 = (-1 / (u * T2 - T1) * (timer.y - T1)) + 1;
      f1 = 1;
  else
    f1 = 1;
    f3 = 1;
  end if;
  connect(lessEqual.y, timer.u) annotation(
    Line(points = {{2, 70}, {16, 70}, {16, 70}, {18, 70}}, color = {255, 0, 255}));
  connect(v, lessEqual.u1) annotation(
    Line(points = {{-100, 0}, {-70, 0}, {-70, 70}, {-22, 70}}, color = {0, 0, 127}));
  connect(limit.y, lessEqual.u2) annotation(
    Line(points = {{-38, 30}, {-30, 30}, {-30, 62}, {-22, 62}, {-22, 62}, {-22, 62}}, color = {0, 0, 127}));
  connect(prod_f1_f2_f3.y, trip) annotation(
    Line(points = {{82, 0}, {102, 0}, {102, 0}, {110, 0}}, color = {0, 0, 127}));
  connect(prod_f1_f2.y, prod_f1_f2_f3.u2) annotation(
    Line(points = {{42, -60}, {48, -60}, {48, -6}, {58, -6}, {58, -6}, {58, -6}}, color = {0, 0, 127}));
  connect(prod_f1_f2.u2, const.y) annotation(
    Line(points = {{18, -66}, {-70, -66}, {-70, -66}, {-68, -66}}, color = {0, 0, 127}));
  connect(min_f3.y, prod_f1_f2_f3.u1) annotation(
    Line(points = {{2, -10}, {40, -10}, {40, 6}, {58, 6}}, color = {0, 0, 127}));
  connect(min_f1.y, prod_f1_f2.u1) annotation(
    Line(points = {{2, -40}, {8, -40}, {8, -54}, {18, -54}}, color = {0, 0, 127}));
  connect(limiter_low.y, min_f1.u) annotation(
    Line(points = {{-34, -40}, {-22, -40}}, color = {0, 0, 127}));
  connect(limiter_low.u, f1) annotation(
    Line(points = {{-58, -40}, {-70, -40}}, color = {0, 0, 127}));
  connect(f3, limiter.u) annotation(
    Line(points = {{-70, -10}, {-58, -10}}, color = {0, 0, 127}));
  connect(limiter.y, min_f3.u) annotation(
    Line(points = {{-34, -10}, {-22, -10}}, color = {0, 0, 127}));
  annotation(
    __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
    Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {-8, 46}, extent = {{-68, 18}, {68, -18}}, textString = "Partial Tripping")}, coordinateSystem(initialScale = 0.1)));
end TrippingLarge;
