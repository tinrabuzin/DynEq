within DynEq.Elements.Solar.WECC;

model PVD1
  parameter SI.PerUnit Imax = 1.1 "Maximum allowable total converter current";
  parameter Boolean PqFlag "Priority on current limit flag: 1=P prio.; 0 = Q prio.";
  parameter SI.Time Tg = 0.02 "Inverter current regulator time constat";
  parameter SI.PerUnit Xc = 0 "Line drop compensation reactance";
  parameter SI.PerUnit Qmx = 0.328 "Maximum reactive power";
  parameter SI.PerUnit Qmn = -0.328 "Minimum reactive power";
  parameter SI.PerUnit v0 = 0.9 "Low voltage threshold for Volt/Var Control";
  parameter SI.PerUnit v1 = 1.1 "High voltage threshold for Volt/Var Control";
  parameter SI.PerUnit dqdv = 0 "Voltage/Var droop compensation";
  parameter SI.PerUnit fdbd = -99 "Frequency deadband over frequency response";
  parameter Real Ddn = 0 "Down regulation droop";
  parameter Real vr_recov = 1 "Amount of generation to reconnect after voltage disconnection";
  parameter Real fr_recov = 1 "Amount of generation to reconnect after frequency disconnection";
  parameter Real Ft0 = 0.99;
  parameter Real Ft1 = 0.995;
  parameter Real Ft2 = 1.005;
  parameter Real Ft3 = 1.01;
  parameter Real Vt0 = 0.88;
  parameter Real Vt1 = 0.9;
  parameter Real Vt2 = 1.1;
  parameter Real Vt3 = 1.2;
  Modelica.Blocks.Interfaces.RealInput Vt annotation(
    Placement(visible = true, transformation(origin = {-200, -10}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-90, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput It annotation(
    Placement(visible = true, transformation(origin = {-200, -70}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain compensation(k = Xc) annotation(
    Placement(visible = true, transformation(origin = {-150, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add add annotation(
    Placement(visible = true, transformation(origin = {-110, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.Limiter numerical_limit(uMax = Modelica.Constants.inf, uMin = 0.01) annotation(
    Placement(visible = true, transformation(origin = {-70, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Division division annotation(
    Placement(visible = true, transformation(origin = {90, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.Solar.WECC.QPPriority qppriority(Imax = Imax, PqFlag = PqFlag) annotation(
    Placement(visible = true, transformation(origin = {130, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Division division1 annotation(
    Placement(visible = true, transformation(origin = {90, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Continuous.FirstOrder PCurrentController(T = Tg, initType = Modelica.Blocks.Types.Init.SteadyState, k = 1, y_start = Pref / u_0) annotation(
    Placement(visible = true, transformation(origin = {182, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Continuous.FirstOrder QCurrentController(T = Tg, initType = Modelica.Blocks.Types.Init.SteadyState, k = -1, y_start = -Qref / u_0) annotation(
    Placement(visible = true, transformation(origin = {182, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput freq annotation(
    Placement(visible = true, transformation(origin = {-200, 150}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-90, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant freq_ref(k = 1) annotation(
    Placement(visible = true, transformation(origin = {50, 150}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant active_power_reference(k = Pref) annotation(
    Placement(visible = true, transformation(origin = {-50, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.DeadZone deadband_voltage(uMax = v1, uMin = v0) annotation(
    Placement(visible = true, transformation(origin = {-70, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain voltage_droop(k = -dqdv) annotation(
    Placement(visible = true, transformation(origin = {-30, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.Limiter limiter(uMax = Qmx, uMin = Qmn) annotation(
    Placement(visible = true, transformation(origin = {50, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add add1 annotation(
    Placement(visible = true, transformation(origin = {14, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant reactive_power_reference(k = Qref) annotation(
    Placement(visible = true, transformation(origin = {-70, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput Ip annotation(
    Placement(visible = true, transformation(origin = {210, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput Iq annotation(
    Placement(visible = true, transformation(origin = {210, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  parameter SI.PerUnit Qref;
  parameter SI.PerUnit Pref;
  parameter SI.PerUnit u_0;
  Modelica.Blocks.Math.Product product1 annotation(
    Placement(visible = true, transformation(origin = {150, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.Solar.WECC.GenerationTripping  voltage_tripping(Lv0 = Vt0, Lv1 = Vt1, Lv2 = Vt2, Lv3 = Vt3) annotation(
    Placement(visible = true, transformation(origin = {30, 110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Product product2 annotation(
    Placement(visible = true, transformation(origin = {70, 110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Product product3 annotation(
    Placement(visible = true, transformation(origin = {110, 130}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Product product4 annotation(
    Placement(visible = true, transformation(origin = {150, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
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
  connect(voltage_tripping.u, Vt) annotation(
    Line(points = {{20, 110}, {-128, 110}, {-128, -10}, {-200, -10}}, color = {0, 0, 127}));
  connect(voltage_tripping.TrpLow, product2.u1) annotation(
    Line(points = {{42, 116}, {58, 116}, {58, 116}, {58, 116}}, color = {0, 0, 127}));
  connect(voltage_tripping.TrpHigh, product2.u2) annotation(
    Line(points = {{42, 106}, {56, 106}, {56, 104}, {58, 104}}, color = {0, 0, 127}));
  connect(product2.y, product3.u2) annotation(
    Line(points = {{82, 110}, {90, 110}, {90, 124}, {98, 124}, {98, 124}}, color = {0, 0, 127}));
  connect(product3.y, product1.u1) annotation(
    Line(points = {{122, 130}, {128, 130}, {128, 96}, {138, 96}, {138, 96}}, color = {0, 0, 127}));
  connect(product4.y, QCurrentController.u) annotation(
    Line(points = {{162, -70}, {168, -70}, {168, -70}, {170, -70}}, color = {0, 0, 127}));
  connect(qppriority.Iqcmd, product4.u1) annotation(
    Line(points = {{142, -4}, {144, -4}, {144, -52}, {130, -52}, {130, -64}, {138, -64}, {138, -64}}, color = {0, 0, 127}));
  connect(product4.u2, product3.y) annotation(
    Line(points = {{138, -76}, {112, -76}, {112, 100}, {128, 100}, {128, 130}, {122, 130}, {122, 130}}, color = {0, 0, 127}));
  connect(active_power_reference.y, division1.u1) annotation(
    Line(points = {{-38, 40}, {60, 40}, {60, 36}, {78, 36}}, color = {0, 0, 127}));
  connect(freq_ref.y, product3.u1) annotation(
    Line(points = {{62, 150}, {78, 150}, {78, 136}, {98, 136}}, color = {0, 0, 127}));
  annotation(
    Diagram(coordinateSystem(extent = {{-200, -200}, {200, 200}}, initialScale = 0.05), graphics = {Rectangle(origin = {-5, -71}, lineColor = {0, 170, 0}, pattern = LinePattern.Dash, lineThickness = 1, extent = {{-85, 21}, {79, -57}}), Text(origin = {-120, -41}, lineColor = {0, 170, 0}, lineThickness = 1, extent = {{30, -1}, {76, -9}}, textString = "Volt/Var Control", horizontalAlignment = TextAlignment.Left)}),
    Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, initialScale = 0.05), graphics = {Rectangle(origin = {-25.3732, 26.25}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-74.6268, 75.75}, {125.373, -126.25}}), Text(origin = {12.5671, 95.7}, extent = {{-114.568, -35.7}, {87.4328, 6.3}}, textString = "PVD1")}),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
    __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
end PVD1;
