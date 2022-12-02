within DynEq.Elements;
model ElmVac
  extends DynEq.Elements.BaseClasses.OnePort;

  SI.Angle phiu;
  Modelica.Blocks.Interfaces.RealInput fPu annotation (
    Placement(visible = true, transformation(origin = {-100, 50}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, 50}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput vPu annotation (
    Placement(visible = true, transformation(origin = {-100, -50}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, -50}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Real P1;
  Real Q1;
initial equation
  phiu = angle_0;
equation
  P1 = -S.re/1e3;
  Q1 = -S.im/1e3;
  der(phiu) = 2*fBase*C.pi*(fPu-1);
  p.v = VBase*Complex(vPu*cos(phiu), vPu*sin(phiu));
annotation (
    Icon(graphics={  Ellipse(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}, endAngle = 360), Line(origin = {-5.53254, -0.428994}, points = {{-40, 0}, {-20, 20}, {20, -20}, {40, 0}, {40, 0}})}));
end ElmVac;
