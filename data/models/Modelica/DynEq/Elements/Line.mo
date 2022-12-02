within DynEq.Elements;
model Line "Model for a transmission Line based on the pi-equivalent circuit"
  extends DynEq.Elements.BaseClasses.TwoPort;
  parameter SI.PerUnit R "Resistance (pu)"
    annotation (Dialog(group="Line parameters"));
  parameter SI.PerUnit X "Reactance (pu)"
    annotation (Dialog(group="Line parameters"));
  parameter Complex Z = Complex(0, X);
equation
  (termA.v - termB.v) = Complex(-termA.i.im*X, termA.i.re*X);
  termB.i = -termA.i;
  annotation (Icon(coordinateSystem(initialScale = 0.1),
        graphics={Rectangle(lineColor = {0, 85, 0},fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-80, 20}, {80, -20}}), Line(origin = {-90, 0}, points = {{10, 0}, {-10, 0}}), Line(origin = {-90, 0}, points = {{10, 0}, {-10, 0}}, color = {0, 85, 0}), Line(origin = {90.39, -0.16}, points = {{10, 0}, {-10, 0}}, color = {0, 85, 0})}));
end Line;
