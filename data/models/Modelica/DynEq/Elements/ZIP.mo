within DynEq.Elements;

model ZIP
  extends DynEq.Elements.BaseClasses.OnePort;
  parameter Real cPi annotation(Dialog(group = "ZIP Model Parameters"));
  parameter Real cPp annotation(Dialog(group = "ZIP Model Parameters"));
  parameter Real cQi annotation(Dialog(group = "ZIP Model Parameters"));
  parameter Real cQp annotation(Dialog(group = "ZIP Model Parameters"));
  parameter SI.PerUnit ratio;
equation
  ratio = CM.'abs'(p.v) / v_0;
  S.re = P_0 * ((1 - cPi - cPp) * ratio ^ 2 + cPi * ratio + cPp);
  S.im = Q_0 * ((1 - cQi - cQp) * ratio ^ 2 + cQi * ratio + cQp);
  annotation(
    Icon(graphics = {Polygon(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, points = {{0, 40}, {60, 40}, {0, -80}, {-60, 40}, {0, 40}})}, coordinateSystem(initialScale = 0.1)));
end ZIP;
