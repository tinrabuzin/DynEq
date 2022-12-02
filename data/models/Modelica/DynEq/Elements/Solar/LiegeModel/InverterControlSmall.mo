within DynEq.Elements.Solar.LiegeModel;

model InverterControlSmall
  parameter SI.PerUnit Imax = 1.1 "Maximum allowable total converter current";
  parameter Boolean PqFlag "Priority on current limit flag: 1=P prio.; 0 = Q prio.";
  parameter SI.Time Tg = 0.02 "Inverter current regulator time constat";
  parameter Real Vmin = 0.2;
  parameter Real c = 0.3;
  parameter Real d = 0.2;
  
  Modelica.Blocks.Interfaces.RealInput Vt annotation(
    Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.Limiter numerical_limit( uMax = Modelica.Constants.inf, uMin = 0.01) annotation(
    Placement(visible = true, transformation(origin = {-50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Division division annotation(
    Placement(visible = true, transformation(origin = {-10, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.Solar.WECC.QPPriority qppriority(Imax = Imax, PqFlag = PqFlag) annotation(
    Placement(visible = true, transformation(origin = {30, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Division division1 annotation(
    Placement(visible = true, transformation(origin = {-10, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Continuous.FirstOrder PCurrentController(T = Tg, initType = Modelica.Blocks.Types.Init.InitialOutput, k = 1, y_start = Pref / u_0) annotation(
    Placement(visible = true, transformation(origin = {82, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Continuous.FirstOrder QCurrentController(T = Tg, initType = Modelica.Blocks.Types.Init.InitialOutput, k = -1, y_start = -Qref / u_0) annotation(
    Placement(visible = true, transformation(origin = {82, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant active_power_reference(k = Pref) annotation(
    Placement(visible = true, transformation(origin = {-50, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant reactive_power_reference(k = Qref) annotation(
    Placement(visible = true, transformation(origin = {-50, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput Ip annotation(
    Placement(visible = true, transformation(origin = {110, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput Iq annotation(
    Placement(visible = true, transformation(origin = {110, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  parameter SI.PerUnit Qref;
  parameter SI.PerUnit Pref;
  parameter SI.PerUnit u_0;
  Modelica.Blocks.Math.Product product1 annotation(
    Placement(visible = true, transformation(origin = {50, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Product product4 annotation(
    Placement(visible = true, transformation(origin = {50, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.Solar.LiegeModel.TrippingSmall trippingSmall( Vmin = Vmin, c = c, d = d)  annotation(
    Placement(visible = true, transformation(origin = {-10, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(numerical_limit.u, Vt) annotation(
    Line(points = {{-62, 0}, {-100, 0}}, color = {0, 0, 127}));
  connect(division.y, qppriority.Iq) annotation(
    Line(points = {{1, -30}, {6, -30}, {6, -5}, {20, -5}}, color = {0, 0, 127}));
  connect(division1.y, qppriority.Ip) annotation(
    Line(points = {{1, 30}, {6, 30}, {6, 5}, {20, 5}}, color = {0, 0, 127}));
  connect(QCurrentController.y, Iq) annotation(
    Line(points = {{93, -60}, {110, -60}}, color = {0, 0, 127}));
  connect(PCurrentController.y, Ip) annotation(
    Line(points = {{93, 60}, {110, 60}}, color = {0, 0, 127}));
  connect(product1.y, PCurrentController.u) annotation(
    Line(points = {{61, 60}, {70, 60}}, color = {0, 0, 127}));
  connect(product4.y, QCurrentController.u) annotation(
    Line(points = {{61, -60}, {70, -60}}, color = {0, 0, 127}));
  connect(trippingSmall.trip, product1.u1) annotation(
    Line(points = {{1, 70}, {23.5, 70}, {23.5, 66}, {38, 66}}, color = {0, 0, 127}));
  connect(reactive_power_reference.y, division.u1) annotation(
    Line(points = {{-39, -50}, {-33.5, -50}, {-33.5, -24}, {-22, -24}}, color = {0, 0, 127}));
  connect(numerical_limit.y, division1.u2) annotation(
    Line(points = {{-38, 0}, {-32, 0}, {-32, 24}, {-22, 24}}, color = {0, 0, 127}));
  connect(division.u2, numerical_limit.y) annotation(
    Line(points = {{-22, -36}, {-32, -36}, {-32, 0}, {-38, 0}}, color = {0, 0, 127}));
  connect(active_power_reference.y, division1.u1) annotation(
    Line(points = {{-39, 50}, {-32, 50}, {-32, 36}, {-22, 36}}, color = {0, 0, 127}));
  connect(trippingSmall.v, Vt) annotation(
    Line(points = {{-20, 70}, {-80, 70}, {-80, 0}, {-100, 0}}, color = {0, 0, 127}));
  connect(product4.u2, trippingSmall.trip) annotation(
    Line(points = {{38, -66}, {14, -66}, {14, 70}, {1, 70}}, color = {0, 0, 127}));
  connect(qppriority.Ipcmd, product1.u2) annotation(
    Line(points = {{42, 4}, {50, 4}, {50, 40}, {30, 40}, {30, 54}, {38, 54}}, color = {0, 0, 127}));
  connect(qppriority.Iqcmd, product4.u1) annotation(
    Line(points = {{42, -4}, {50, -4}, {50, -32}, {30, -32}, {30, -54}, {38, -54}}, color = {0, 0, 127}));
  annotation(
    Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, initialScale = 0.1)),
    Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, initialScale = 0.05), graphics = {Rectangle(origin = {-25.3732, 25}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-74.6271, 75}, {125.373, -125}}), Text(origin = {13.4328, 94}, extent = {{-113.433, -34}, {86.5672, 6}}, textString = "PVSmall")}),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
    __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
end InverterControlSmall;
