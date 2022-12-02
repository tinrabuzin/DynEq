within DynEq.Elements.BaseClasses;

partial model TwoPort
  outer DynEq.Essentials.SystemBase SysData;
  parameter Types.Voltage VBaseA annotation(Evaluate=true);
  parameter Types.Voltage VBaseB annotation(Evaluate=true);
  parameter Types.ApparentPower SNomA=SysData.SBase annotation(Evaluate=true);
  parameter Types.ApparentPower SNomB=SysData.SBase annotation(Evaluate=true);
  
  parameter Types.Voltage VStartB;
  parameter Types.Angle angleStartB;
  parameter Types.ActivePower PStartB;
  parameter Types.ReactivePower QStartB;
  parameter Types.ComplexCurrent IStartB = CM.conj(Complex(PStartB, QStartB)/Complex(VStartB*cos(angleStartB), VStartB*sin(angleStartB)));

  parameter Types.Voltage VStartA;
  parameter Types.Angle angleStartA;
  parameter Types.ActivePower PStartA;
  parameter Types.ReactivePower QStartA;
  parameter Types.ComplexCurrent IStartA = CM.conj(Complex(PStartA, QStartA)/Complex(VStartA*cos(angleStartA), VStartA*sin(angleStartA)));


  Types.ComplexPower S_AB=termA.v*CM.conj(termA.i);
  Types.ComplexPower S_BA=termB.v*CM.conj(termB.i);
  DynEq.Interfaces.terminal termA(
  v(re(start=VStartA*cos(angleStartA), nominal=VBaseA), im(start=VStartA*sin(angleStartA), nominal=VBaseA)),
  i(re(start=IStartA.re,nominal=SNomA/VBaseA), im(start = IStartA.im,nominal=SNomA/VBaseA))
  ) annotation(
    Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-106, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  DynEq.Interfaces.terminal termB(
  v(re(start=VStartB*cos(angleStartB), nominal=VBaseB), im(start=VStartB*sin(angleStartB), nominal=VBaseB)),
  i(re(start=IStartB.re, nominal=SNomB/VBaseB), im(start = IStartB.im, nominal=SNomB/VBaseB))
  ) annotation(
    Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {106, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
end TwoPort;
