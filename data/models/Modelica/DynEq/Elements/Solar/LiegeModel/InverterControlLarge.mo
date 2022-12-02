within DynEq.Elements.Solar.LiegeModel;

model InverterControlLarge
  parameter SI.PerUnit Imax = 1.1 "Maximum allowable total converter current";
  parameter Boolean PqFlag "Priority on current limit flag: 1=P prio.; 0 = Q prio.";
  parameter SI.Time Tg = 0.02 "Inverter current regulator time constat";
  parameter SI.PerUnit Xc = 0 "Line drop compensation reactance";
  parameter SI.PerUnit Qmx = 0.328 "Maximum reactive power";
  parameter SI.PerUnit Qmn = -0.328 "Minimum reactive power";
  parameter SI.PerUnit v0 = 0.9 "Low voltage threshold for Volt/Var Control";
  parameter SI.PerUnit v1 = 1.1 "High voltage threshold for Volt/Var Control";
  parameter SI.PerUnit dqdv = 0 "Voltage/Var droop compensation";

  parameter Real Vr = 0.9;
  parameter Real Vmin = 0.2;
  parameter Real c = 0.3;
  parameter Real d = 0.2;
  parameter Real u = 1;
  parameter Real T1 = 0.2;
  parameter Real T2 = 0.7;
  
  Modelica.Blocks.Interfaces.RealInput Vt annotation(
    Placement(visible = true, transformation(origin = {-200, -10}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-200, 140}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput It annotation(
    Placement(visible = true, transformation(origin = {-200, -70}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-200, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Math.Gain compensation(k = Xc) annotation(
    Placement(visible = true, transformation(origin = {-150, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add add annotation(
    Placement(visible = true, transformation(origin = {-110, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.Limiter numerical_limit( uMax = Modelica.Constants.inf, uMin = 0.01) annotation(
    Placement(visible = true, transformation(origin = {-70, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Division division annotation(
    Placement(visible = true, transformation(origin = {90, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.Solar.WECC.QPPriority qppriority(Imax = Imax, PqFlag = PqFlag) annotation(
    Placement(visible = true, transformation(origin = {130, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Division division1 annotation(
    Placement(visible = true, transformation(origin = {90, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Continuous.FirstOrder PCurrentController(T = Tg, initType = Modelica.Blocks.Types.Init.InitialOutput, k = 1, y_start = Pref / u_0) annotation(
    Placement(visible = true, transformation(origin = {182, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Continuous.FirstOrder QCurrentController(T = Tg, initType = Modelica.Blocks.Types.Init.InitialOutput, k = -1, y_start = -Qref / u_0) annotation(
    Placement(visible = true, transformation(origin = {182, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant active_power_reference(k = Pref) annotation(
    Placement(visible = true, transformation(origin = {-10, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.DeadZone deadband_voltage(uMax = v1, uMin = v0) annotation(
    Placement(visible = true, transformation(origin = {-70, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain voltage_droop(k = dqdv) annotation(
    Placement(visible = true, transformation(origin = {-30, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.Limiter limiter(uMax = Qmx, uMin = Qmn) annotation(
    Placement(visible = true, transformation(origin = {50, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add add1 annotation(
    Placement(visible = true, transformation(origin = {14, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant reactive_power_reference(k = Qref) annotation(
    Placement(visible = true, transformation(origin = {-70, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput Ip annotation(
    Placement(visible = true, transformation(origin = {210, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {210, 110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput Iq annotation(
    Placement(visible = true, transformation(origin = {210, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {210, -90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  parameter SI.PerUnit Qref;
  parameter SI.PerUnit Pref;
  parameter SI.PerUnit u_0;
  Modelica.Blocks.Math.Product product1 annotation(
    Placement(visible = true, transformation(origin = {150, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Product product4 annotation(
    Placement(visible = true, transformation(origin = {150, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.Solar.LiegeModel.TrippingLarge trippingLarge(T1 = T1, T2 = T2, Vmin = Vmin, Vr = Vr, c = c, d = d, u = u)  annotation(
    Placement(visible = true, transformation(origin = {70, 150}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(It, compensation.u) annotation(
    Line(points = {{-200, -70}, {-162, -70}}, color = {0, 0, 127}));
  connect(compensation.y, add.u2) annotation(
    Line(points = {{-139, -70}, {-130.5, -70}, {-130.5, -76}, {-122, -76}}, color = {0, 0, 127}));
  connect(Vt, add.u1) annotation(
    Line(points = {{-200, -10}, {-128, -10}, {-128, -63}, {-122, -63}, {-122, -64}}, color = {0, 0, 127}));
  connect(numerical_limit.u, Vt) annotation(
    Line(points = {{-82, -10}, {-200, -10}}, color = {0, 0, 127}));
  connect(numerical_limit.y, division.u2) annotation(
    Line(points = {{-59, -10}, {-8, -10}, {-8, -36}, {78, -36}}, color = {0, 0, 127}));
  connect(division.y, qppriority.Iq) annotation(
    Line(points = {{101, -30}, {106, -30}, {106, -5}, {120, -5}}, color = {0, 0, 127}));
  connect(division1.u2, numerical_limit.y) annotation(
    Line(points = {{78, 24}, {-8, 24}, {-8, -10}, {-59, -10}}, color = {0, 0, 127}));
  connect(division1.y, qppriority.Ip) annotation(
    Line(points = {{101, 30}, {106, 30}, {106, 5}, {120, 5}}, color = {0, 0, 127}));
  connect(add.y, deadband_voltage.u) annotation(
    Line(points = {{-99, -70}, {-82, -70}}, color = {0, 0, 127}));
  connect(voltage_droop.u, deadband_voltage.y) annotation(
    Line(points = {{-42, -70}, {-59, -70}}, color = {0, 0, 127}));
  connect(add1.y, limiter.u) annotation(
    Line(points = {{25, -70}, {38, -70}}, color = {0, 0, 127}));
  connect(voltage_droop.y, add1.u1) annotation(
    Line(points = {{-19, -70}, {-10, -70}, {-10, -69}, {2, -69}, {2, -64}}, color = {0, 0, 127}));
  connect(reactive_power_reference.y, add1.u2) annotation(
    Line(points = {{-59, -110}, {-10, -110}, {-10, -95}, {2, -95}, {2, -76}}, color = {0, 0, 127}));
  connect(limiter.y, division.u1) annotation(
    Line(points = {{61, -70}, {62.5, -70}, {62.5, -24}, {78, -24}}, color = {0, 0, 127}));
  connect(QCurrentController.y, Iq) annotation(
    Line(points = {{193, -70}, {210, -70}}, color = {0, 0, 127}));
  connect(PCurrentController.y, Ip) annotation(
    Line(points = {{193, 90}, {210, 90}}, color = {0, 0, 127}));
  connect(product1.y, PCurrentController.u) annotation(
    Line(points = {{162, 90}, {170, 90}, {170, 90}, {170, 90}}, color = {0, 0, 127}));
  connect(product1.u2, qppriority.Ipcmd) annotation(
    Line(points = {{138, 84}, {126, 84}, {126, 20}, {144, 20}, {144, 6}, {142, 6}, {142, 6}}, color = {0, 0, 127}));
  connect(product4.y, QCurrentController.u) annotation(
    Line(points = {{162, -70}, {168, -70}, {168, -70}, {170, -70}}, color = {0, 0, 127}));
  connect(qppriority.Iqcmd, product4.u1) annotation(
    Line(points = {{142, -4}, {144, -4}, {144, -52}, {130, -52}, {130, -64}, {138, -64}, {138, -64}}, color = {0, 0, 127}));
  connect(trippingLarge.v, Vt) annotation(
    Line(points = {{60, 150}, {20, 150}, {20, 80}, {-186, 80}, {-186, -10}, {-200, -10}}, color = {0, 0, 127}));
  connect(trippingLarge.trip, product1.u1) annotation(
    Line(points = {{82, 150}, {112, 150}, {112, 96}, {138, 96}}, color = {0, 0, 127}));
  connect(product4.u2, trippingLarge.trip) annotation(
    Line(points = {{138, -76}, {112, -76}, {112, 150}, {80, 150}, {80, 150}, {82, 150}}, color = {0, 0, 127}));
  connect(active_power_reference.y, division1.u1) annotation(
    Line(points = {{2, 50}, {60, 50}, {60, 36}, {78, 36}, {78, 36}}, color = {0, 0, 127}));
  annotation(
    Diagram(coordinateSystem(extent = {{-200, -200}, {200, 200}}, initialScale = 0.05), graphics = {Rectangle(origin = {-5, -71}, lineColor = {0, 170, 0}, pattern = LinePattern.Dash, lineThickness = 1, extent = {{-85, 21}, {79, -57}}), Text(origin = {-120, -41}, lineColor = {0, 170, 0}, lineThickness = 1, extent = {{30, -1}, {76, -9}}, textString = "Volt/Var Control", horizontalAlignment = TextAlignment.Left)}),
    Icon(coordinateSystem(extent = {{-200, -200}, {200, 200}}, initialScale = 0.05), graphics = {Rectangle(origin = {-52, 50}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-150, 150}, {252, -250}}), Text(origin = {26, 194}, extent = {{-228, -34}, {174, 6}}, textString = "PVLarge")}),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
    __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
end InverterControlLarge;
