within DynEq.Elements;

model ExponentialLoad
  extends DynEq.Elements.BaseClasses.OnePort;
  parameter Real expP = 2;
  parameter Real expQ = 2;
equation
  S = Complex(P_0 * ( CM.'abs'(p.v) / v_0) ^ expP, Q_0 * ( CM.'abs'(p.v) / v_0) ^ expQ);
  annotation(
    Icon(graphics = {Polygon(origin = {0, -20},fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, points = {{0, 60}, {60, 60}, {0, -60}, {-60, 60}, {0, 60}})}, coordinateSystem(initialScale = 0.1)));
end ExponentialLoad;
