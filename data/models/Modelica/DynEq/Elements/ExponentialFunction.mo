within DynEq.Elements;

model ExponentialFunction
  extends Modelica.Blocks.Interfaces.SISO;
  parameter Real u_0;
  parameter Real exponent;
equation
y = (u/u_0)^exponent;
end ExponentialFunction;
