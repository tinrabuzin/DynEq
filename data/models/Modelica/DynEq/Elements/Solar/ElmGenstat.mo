within DynEq.Elements.Solar;

model ElmGenstat
  extends DynEq.Elements.BaseClasses.OnePort(p(i(re(nominal=M_b/VBase),im(nominal=M_b/VBase))));
  parameter Types.ApparentPower M_b "Base power of the PV generation" annotation (Dialog(group="Power flow data"));
  Modelica.Blocks.Interfaces.RealInput id_ref annotation(
    Placement(visible = true, transformation(origin = {-108, 70}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput iq_ref annotation(
    Placement(visible = true, transformation(origin = {-108, 30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput i annotation(
    Placement(visible = true, transformation(origin = {110, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput omega annotation(
    Placement(visible = true, transformation(origin = {0, 100}, extent = {{-20, -20}, {20, 20}}, rotation = -90), iconTransformation(origin = {0, 98}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Types.Angle phi;
  Modelica.Blocks.Interfaces.RealOutput v annotation(
    Placement(visible = true, transformation(origin = {110, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
initial equation
phi = angle_0;
equation
  der(phi) = 2*Modelica.Constants.pi*50*(omega-1);
  p.i.re/(M_b/VBase) = -(id_ref * cos(phi) - iq_ref * sin(phi));
  p.i.im/(M_b/VBase) = -(id_ref * sin(phi) + iq_ref * cos(phi));
  i = CM.'abs'(p.i)/(M_b/VBase);
  v = CM.'abs'(p.v)/VBase;
  annotation(
    Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {100, 100}}), Text(origin = {-69.7, 50.01}, fillPattern = FillPattern.Solid, extent = {{-12.3, -8.18}, {12.3, 8.18}}, textString = "id_ref", fontName = "Arial"), Text(origin = {-63.7, -50.18}, fillPattern = FillPattern.Solid, extent = {{-12.3, -8.18}, {12.3, 8.18}}, textString = "iq_ref", fontName = "Arial"), Text(origin = {24.07, -92.08}, fillPattern = FillPattern.Solid, extent = {{-124.07, -27.92}, {76.07, -8.08}}, textString = "ElmStatgen", fontName = "Arial"), Ellipse(origin = {50, -48}, fillColor = {215, 215, 215}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {-2, -2}}, endAngle = 360), Line(origin = {-0.000710326, 0.400639}, points = {{0, -25}, {0, 25}}, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 10)}, coordinateSystem(initialScale = 0.1)),
    Diagram);
end ElmGenstat;
