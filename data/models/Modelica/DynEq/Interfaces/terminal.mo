within DynEq.Interfaces;

connector terminal "Connector for electrical blocks treating voltage and current as complex variables"
  Types.ComplexVoltage v(re(nominal=20e3, start=20e3), im(nominal=20e3, start=0)) "Voltage phasor";
  flow Types.ComplexCurrent i(re(nominal=1e6/20e3),im(nominal=1e6/20e3)) "Current phasor";
  annotation(
    Icon(graphics = {Rectangle(lineColor = {0, 85, 0}, fillColor = {0, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}})}),
    Diagram(graphics = {Text(lineColor = {0, 0, 255}, extent = {{-100, 160}, {100, 120}}, textString = "%name"), Rectangle(lineColor = {0, 85, 0}, fillColor = {0, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}})}));
end terminal;
