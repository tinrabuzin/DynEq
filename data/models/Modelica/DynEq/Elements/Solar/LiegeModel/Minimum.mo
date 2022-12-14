within DynEq.Elements.Solar.LiegeModel;

block Minimum
  extends Modelica.Blocks.Interfaces.SISO;
  parameter Real y_start=Modelica.Constants.inf "Initialization of minimum";
  replaceable model Last = LastLib.Last_omc constrainedby LastLib.LastBase(y_start=y_start) "Last FMU"
    annotation(choices(
     choice(redeclare model Last = LastLib.Last_dymola "Last for Dymola (Windows)"),
     choice(redeclare model Last = LastLib.Last_omc "Last for OpenModelica (C runtime)"),
     choice(redeclare model Last = LastLib.Last_simx "Last for SimulationX")));
  Last last(y_start=y_start, logLevel=0, debugLogging=false) "Last FMU" annotation(Placement(transformation(extent={{-30,10},{-10,30}})));
  Modelica.Blocks.Math.Min min annotation(Placement(transformation(extent={{10,-10},{30,10}})));
  equation
    connect(u, min.u2) annotation(Line(points={{-120,0},{0,0},{0,-6},{8,-6}}, color={0,0,127}));
    connect(min.y, y) annotation(Line(points={{31,0},{110,0}}, color={0,0,127}));
    connect(last.y, min.u1) annotation(Line(points={{-9,20},{0,20},{0,6},{8,6}}, color={0,0,127}));
    connect(min.y, last.u) annotation(Line(points={{31,0},{40,0},{40,40},{-40,40},{-40,20},{-32,20}}, color={0,0,127}));
  annotation(defaultComponentName="min", preferredView="diagram");
end Minimum;
