within DynEq.Equivalents;

model M6
  extends Modelica.Icons.Example;
  // General Parameters
  parameter Real P_0 = 17.593473230782624e3 ;
  parameter Real Q_0 = 6.44994877688119e3;
  parameter Real v_0 = 1.0003343;
  parameter Real angle_0 = 1;
  parameter Real omega_0 = 0.9997503;
parameter Real fP1 = 0.9294315989160804;
parameter Real fP2 = 5.4851040477831905;
parameter Real fP3 = 3.594410123807072;
parameter Real fP4 = 3.272961387290894;
parameter Real fP5 = 2.0880318042527035;
parameter Real fP6 = 10.841119265101831;
parameter Real vP1 = 2.85992291101969e-09;
parameter Real vP2 = 35.350424103251555;
parameter Real vP3 = 1.7060552497789125;
parameter Real vP4 = 21.289236392865593;
parameter Real vP5 = 85.37454172784466;
parameter Real vP6 = 2.6378594582453543;
parameter Real fQ1 = 14.377025592183436;
parameter Real fQ2 = 24.38909008196476;
parameter Real fQ3 = 12.045096329010525;
parameter Real fQ4 = 756.7456193006351;
parameter Real fQ5 = 255.42075679667417;
parameter Real fQ6 = 33.62408166790542;
parameter Real vQ1 = 1.4861548236676412;
parameter Real vQ2 = 5.3186908455349435;
parameter Real vQ3 = 1.26502429923883;
parameter Real vQ4 = 0.044426096242901836;
parameter Real vQ5 = 11.85201719310349;
parameter Real vQ6 = 2.36029766411587;

//parameter Real deadZone.uMin=0.33514144596496576;
  //parameter Real gain2.k=0.00065755930718596;
  //parameter Real deadZone1.uMin=-0.014252939026995409;
  //parameter Real gain5.k=-0.007892370572483634;
  // Components
  Modelica.Blocks.Interfaces.RealInput u1 annotation(
    Placement(visible = true, transformation(origin = {-200, -78}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(extent = {{-120, 40}, {-80, 80}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput u2 annotation(
    Placement(visible = true, transformation(origin = {-200, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(extent = {{-120, -80}, {-80, -40}}, rotation = 0)));
  Modelica.Blocks.Math.Add add(k2 = -1) annotation(
    Placement(visible = true, transformation(origin = {-150, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const(k = omega_0) annotation(
    Placement(visible = true, transformation(origin = {-190, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const1(k = P_0 / 1e3) annotation(
    Placement(visible = true, transformation(origin = {150, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add add2(k2 = -1) annotation(
    Placement(visible = true, transformation(origin = {-150, -108}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant constant1(k = v_0) annotation(
    Placement(visible = true, transformation(origin = {-190, -128}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add3 intern_P annotation(
    Placement(visible = true, transformation(origin = {110, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add active_power annotation(
    Placement(visible = true, transformation(origin = {190, 62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain gain3(k = 1 / 0.01) annotation(
    Placement(visible = true, transformation(origin = {-120, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain gain4(k = 1 / 0.1) annotation(
    Placement(visible = true, transformation(origin = {-120, -108}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Continuous.TransferFunction f_P(a = {fP6, fP5, fP4}, b = {fP3, fP2, fP1}, initType = Modelica.Blocks.Types.Init.InitialOutput) annotation(
    Placement(visible = true, transformation(origin = {-50, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Continuous.TransferFunction v_Q(a = {vQ6, vQ5, vQ4}, b = {vQ2, vQ1}, initType = Modelica.Blocks.Types.Init.InitialOutput) annotation(
    Placement(visible = true, transformation(origin = {-50, -108}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Continuous.TransferFunction f_Q(a = {fQ6, fQ5, fQ4}, b = {fQ3, fQ2, fQ1}, initType = Modelica.Blocks.Types.Init.InitialOutput) annotation(
    Placement(visible = true, transformation(origin = {-50, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add3 intern_Q annotation(
    Placement(visible = true, transformation(origin = {110, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add reactive_power annotation(
    Placement(visible = true, transformation(origin = {190, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const3(k = Q_0 / 1e3) annotation(
    Placement(visible = true, transformation(origin = {150, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Continuous.TransferFunction v_P(a = { vP6, vP5, vP4}, b = {vP3, vP2, vP1}, initType = Modelica.Blocks.Types.Init.InitialOutput) annotation(
    Placement(visible = true, transformation(origin = {-50, -68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.Solar.LiegeModel.Minimum min1(y_start = 0)  annotation(
    Placement(visible = true, transformation(origin = {-12, 102}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.DeadZone trip_ch_P(deadZoneAtInit = true, uMax=Modelica.Constants.inf, uMin = -1.3418481497282997e-08) annotation(
    Placement(visible = true, transformation(origin = {-50, 102}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain trip_gain_P(k = 0) annotation(
    Placement(visible = true, transformation(origin = {28, 102}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.DeadZone trip_ch_Q(deadZoneAtInit = true, uMax=Modelica.Constants.inf, uMin = -0.28312589198398463) annotation(
    Placement(visible = true, transformation(origin = {-38, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain trip_gain_Q(k = 0) annotation(
    Placement(visible = true, transformation(origin = {40, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.Solar.LiegeModel.Minimum minimum(y_start = 0) annotation(
    Placement(visible = true, transformation(origin = {0, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain gain_vP(k = -3.2596757192924066) annotation(
    Placement(visible = true, transformation(origin = {-4, -68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain gain_vQ(k = -3.2596757192924066) annotation(
    Placement(visible = true, transformation(origin = {-4, -108}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add add1 annotation(
    Placement(visible = true, transformation(origin = {50, -108}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.ExponentialFunction static_load_Q(exponent = 2, u_0 = v_0) annotation(
    Placement(visible = true, transformation(origin = {-96, -148}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain gain_exponential_Q(k = 0.013872105498098692) annotation(
    Placement(visible = true, transformation(origin = {-10, -148}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add add4(k2 = -1)  annotation(
    Placement(visible = true, transformation(origin = {-58, -154}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const2(k = 1) annotation(
    Placement(visible = true, transformation(origin = {-96, -182}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add add3(k2 = -1) annotation(
    Placement(visible = true, transformation(origin = {-14, 168}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant constant2(k = 1) annotation(
    Placement(visible = true, transformation(origin = {-52, 140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain gain_exponential_P(k = 0.013872105498098692) annotation(
    Placement(visible = true, transformation(origin = {34, 174}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Elements.ExponentialFunction static_load_P(exponent = 2, u_0 = v_0) annotation(
    Placement(visible = true, transformation(origin = {-52, 174}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add add5 annotation(
    Placement(visible = true, transformation(origin = {66, 126}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(u2, add.u1) annotation(
    Line(points = {{-200, 72}, {-178, 72}, {-178, 78}, {-162, 78}}, color = {0, 0, 127}));
  connect(const.y, add.u2) annotation(
    Line(points = {{-179, 42}, {-175, 42}, {-175, 66}, {-163, 66}}, color = {0, 0, 127}));
  connect(constant1.y, add2.u2) annotation(
    Line(points = {{-179, -128}, {-175, -128}, {-175, -114}, {-162, -114}}, color = {0, 0, 127}));
  connect(u1, add2.u1) annotation(
    Line(points = {{-200, -78}, {-174, -78}, {-174, -102}, {-162, -102}}, color = {0, 0, 127}));
  connect(add.y, gain3.u) annotation(
    Line(points = {{-139, 72}, {-132, 72}}, color = {0, 0, 127}));
  connect(min1.y, trip_gain_P.u) annotation(
    Line(points = {{-1, 102}, {16, 102}}, color = {0, 0, 127}));
  connect(minimum.y, trip_gain_Q.u) annotation(
    Line(points = {{11, -8}, {28, -8}}, color = {0, 0, 127}));
  connect(add2.y, gain4.u) annotation(
    Line(points = {{-139, -108}, {-133, -108}}, color = {0, 0, 127}));
  connect(f_P.u, gain3.y) annotation(
    Line(points = {{-62, 72}, {-108, 72}}, color = {0, 0, 127}));
  connect(f_Q.u, gain3.y) annotation(
    Line(points = {{-62, 42}, {-80, 42}, {-80, 72}, {-108, 72}}, color = {0, 0, 127}));
  connect(const1.y, active_power.u2) annotation(
    Line(points = {{161, 32}, {167, 32}, {167, 56}, {177, 56}}, color = {0, 0, 127}));
  connect(f_Q.y, intern_Q.u1) annotation(
    Line(points = {{-39, 42}, {91, 42}, {91, 0}, {97, 0}}, color = {0, 0, 127}));
  connect(f_P.y, intern_P.u2) annotation(
    Line(points = {{-39, 72}, {97, 72}}, color = {0, 0, 127}));
  connect(const3.y, reactive_power.u2) annotation(
    Line(points = {{161, -48}, {167, -48}, {167, -24}, {177, -24}}, color = {0, 0, 127}));
  connect(trip_gain_Q.y, intern_Q.u2) annotation(
    Line(points = {{51, -8}, {98, -8}}, color = {0, 0, 127}));
  connect(trip_ch_Q.y, minimum.u) annotation(
    Line(points = {{-27, -8}, {-13, -8}}, color = {0, 0, 127}));
  connect(trip_ch_P.u, gain4.y) annotation(
    Line(points = {{-62, 102}, {-104, 102}, {-104, -108}, {-108, -108}}, color = {0, 0, 127}));
  connect(trip_ch_P.y, min1.u) annotation(
    Line(points = {{-39, 102}, {-24, 102}}, color = {0, 0, 127}));
  connect(trip_ch_Q.u, gain4.y) annotation(
    Line(points = {{-50, -8}, {-104, -8}, {-104, -108}, {-108, -108}}, color = {0, 0, 127}));
  connect(v_P.y, gain_vP.u) annotation(
    Line(points = {{-39, -68}, {-17, -68}}, color = {0, 0, 127}));
  connect(v_Q.y, gain_vQ.u) annotation(
    Line(points = {{-39, -108}, {-17, -108}}, color = {0, 0, 127}));
  connect(gain_vP.y, intern_P.u3) annotation(
    Line(points = {{7, -68}, {83, -68}, {83, 64}, {97, 64}}, color = {0, 0, 127}));
  connect(gain_vQ.y, add1.u1) annotation(
    Line(points = {{7, -108}, {37, -108}, {37, -102}}, color = {0, 0, 127}));
  connect(add1.y, intern_Q.u3) annotation(
    Line(points = {{61, -108}, {91, -108}, {91, -16}, {97, -16}}, color = {0, 0, 127}));
  connect(static_load_Q.u, u1) annotation(
    Line(points = {{-108, -148}, {-174, -148}, {-174, -78}, {-200, -78}}, color = {0, 0, 127}));
  connect(gain_exponential_Q.y, add1.u2) annotation(
    Line(points = {{1, -148}, {19, -148}, {19, -114}, {37, -114}}, color = {0, 0, 127}));
  connect(static_load_Q.y, add4.u1) annotation(
    Line(points = {{-85, -148}, {-71, -148}}, color = {0, 0, 127}));
  connect(add4.y, gain_exponential_Q.u) annotation(
    Line(points = {{-47, -154}, {-35, -154}, {-35, -148}, {-23, -148}}, color = {0, 0, 127}));
  connect(const2.y, add4.u2) annotation(
    Line(points = {{-85, -182}, {-81, -182}, {-81, -160}, {-71, -160}}, color = {0, 0, 127}));
  connect(add3.y, gain_exponential_P.u) annotation(
    Line(points = {{-3, 168}, {9, 168}, {9, 174}, {21, 174}}, color = {0, 0, 127}));
  connect(constant2.y, add3.u2) annotation(
    Line(points = {{-41, 140}, {-37, 140}, {-37, 162}, {-27, 162}}, color = {0, 0, 127}));
  connect(static_load_P.u, u1) annotation(
    Line(points = {{-64, 174}, {-200, 174}, {-200, -78}}, color = {0, 0, 127}));
  connect(static_load_P.y, add3.u1) annotation(
    Line(points = {{-41, 174}, {-27, 174}}, color = {0, 0, 127}));
  connect(add5.y, intern_P.u1) annotation(
    Line(points = {{77, 126}, {87, 126}, {87, 80}, {97, 80}}, color = {0, 0, 127}));
  connect(add5.u2, trip_gain_P.y) annotation(
    Line(points = {{54, 120}, {44, 120}, {44, 102}, {40, 102}}, color = {0, 0, 127}));
  connect(add5.u1, gain_exponential_P.y) annotation(
    Line(points = {{54, 132}, {46, 132}, {46, 162}, {62, 162}, {62, 174}, {46, 174}}, color = {0, 0, 127}));
  connect(intern_Q.y, reactive_power.u1) annotation(
    Line(points = {{122, -8}, {166, -8}, {166, -12}, {178, -12}}, color = {0, 0, 127}));
  connect(intern_P.y, active_power.u1) annotation(
    Line(points = {{122, 72}, {168, 72}, {168, 68}, {178, 68}}, color = {0, 0, 127}));
  connect(gain4.y, v_Q.u) annotation(
    Line(points = {{-108, -108}, {-62, -108}}, color = {0, 0, 127}));
  connect(v_P.u, gain4.y) annotation(
    Line(points = {{-62, -68}, {-104, -68}, {-104, -108}, {-108, -108}}, color = {0, 0, 127}));
  annotation(
    __OpenModelica_simulationFlags(csvInput = "C:/Users/trabu/dev/DynEq/data/experiments/exp3/inputs_om.csv", lv = "LOG_STATS", s = "rungekutta"),
    experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-6, Interval = 0.01),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian --replaceHomotopy=actual",
    Diagram(coordinateSystem(extent = {{-200, -200}, {200, 200}})));
end M6;