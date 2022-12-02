within DynEq.Elements.BaseClasses;

partial model OnePort
  outer DynEq.Essentials.SystemBase SysData;
  parameter Types.ApparentPower SBase = SysData.SBase "System base power"  annotation (Evaluate = true);
  parameter SI.Frequency fBase=SysData.fBase "System frequency" annotation (Evaluate = true, Dialog(group="Power flow data"));  
  parameter Types.Voltage VBase "Base voltage of the element" annotation (Evaluate = true, Dialog(group="Power flow data"));
  
  parameter Types.ActivePower P_0 "Initial active power" annotation (Dialog(group="Power flow data"));
  parameter Types.ReactivePower Q_0 "Initial reactive power" annotation (Dialog(group="Power flow data"));
  parameter Types.Voltage v_0 "Initial voltage magnitude (pu)" annotation (Dialog(group="Power flow data"));
  parameter Types.Angle angle_0 "Initial voltage angle" annotation (Dialog(group="Power flow data"));
  
  parameter Types.ComplexVoltage U_0 = Complex(v_0*cos(angle_0), v_0*sin(angle_0));
  parameter Types.ComplexPower S_0 = Complex(P_0, Q_0);
  parameter Types.ComplexCurrent I_0 = CM.conj(S_0/U_0);
  DynEq.Interfaces.terminal p(
  v(re(nominal = VBase, start=U_0.re), im(nominal = VBase, start=U_0.im)),
  i(re(nominal = SBase/VBase, start=I_0.re), im(nominal = SBase/VBase, start=I_0.im))) annotation(
    Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

  Types.ActivePower P = S.re;
  Types.ReactivePower Q = S.im;
  Types.ComplexPower S=p.v*CM.conj(p.i);
end OnePort;
