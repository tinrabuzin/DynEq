within DynEq.Essentials;
record SystemBase "System Base Definition"
  parameter SI.ApparentPower SBase(displayUnit="MVA")=100e6 "System base" annotation(Evaluate = true);
  parameter SI.Frequency fBase=50 "System frequency" annotation(Evaluate = true);
  annotation (
    Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-120,-100},{120,100}}), graphics={Rectangle(extent={{-120,100},{120,-100}},
          lineColor={28,108,200}),Text(
          extent={{-100,40},{100,0}},
          lineColor={28,108,200},
          horizontalAlignment=TextAlignment.Left,
          textString="System Base: %S_b"),Text(
          extent={{-100,-20},{100,-60}},
          lineColor={28,108,200},
          horizontalAlignment=TextAlignment.Left,
          textString="Frequency: %fn"),Text(
          extent={{-100,100},{100,60}},
          lineColor={28,108,200},
          horizontalAlignment=TextAlignment.Center,
          textString="System Data")}),
    defaultComponentName="SysData",
    defaultComponentPrefixes="inner",
    missingInnerMessage="
No 'System Data' component is defined. A default component will be used, and generate a system base of 100 MVA, and a frequency of 50 Hz",
    Diagram(coordinateSystem(
        extent={{-120,-100},{120,100}},
        preserveAspectRatio=false)),
    defaultComponentPrefixes="inner");

end SystemBase;
