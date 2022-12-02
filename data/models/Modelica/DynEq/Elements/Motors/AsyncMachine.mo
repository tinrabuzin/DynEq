within DynEq.Elements.Motors;

model AsyncMachine
  import SI = Modelica.SIunits;
  // Initialization of the machine as in
  //J. Powell, G. Radman, 'Initialization for Dynamic Simulation of Stressed Power Systems Considering Induction Motor Components of Loads'
  // Power flow results
  parameter Real P_0 "Initial active power";
  parameter Real v_0 "Initial voltage";
  parameter Real angle_0 "Initial angle";
  parameter Real omega_0 "Initial frequency";
  parameter Real S_b;
  // Parameters of the induction machine
  parameter Real Mlf = 0.8 annotation(
    Dialog(group = "Machine Parameters"));
  parameter Real Rs = 0.028 annotation(
    Dialog(group = "Machine Parameters"));
  parameter Real Xs = 0.01 annotation(
    Dialog(group = "Machine Parameters"));
  parameter Real Xm = 4.236 annotation(
    Dialog(group = "Machine Parameters"));
  parameter Real Rr = 0.00489 annotation(
    Dialog(group = "Machine Parameters"));
  parameter Real Xr = 0.1822 annotation(
    Dialog(group = "Machine Parameters"));
  parameter Real H = 1.67 annotation(
    Dialog(group = "Machine Parameters"));
  parameter Real A = 1 annotation(
    Dialog(group = "Machine Parameters"));
  parameter Real B = 0 annotation(
    Dialog(group = "Machine Parameters"));
  parameter Real C = 1 - A - B;
  // Machine speeds
  Real w;
  // w.r.t. the frame rotating with win
  Real wr;
  // w.r.t the frame rotating with 50Hz, i..e actual mechanical speed
  // d-q values of the machine model
  Real Edp;
  Real Eqp;
  Real Id;
  Real Iq;
  Real Vd;
  Real Vq;
  // Motor active and reactive powers
  Real P;
  Real Q;
  Real Tm;
  Real Te;
  // Angle w.r.t. the reference angle in the 50Hz rotating frame
  Real delta;
  // Power connector
  DynEq.Elements.PwPin p annotation(
    Placement(visible = true, transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput win annotation(
    Placement(visible = true, transformation(origin = {-100, 70}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
protected
  parameter Real M_b = P_0 * S_b / Mlf / v_0;
  // Transformed machine parameters
  parameter Real Tp = (Xr + Xm) / Rr;
  parameter Real X = Xs + Xm;
  parameter Real Xp = Xs + Xm * Xr / (Xm + Xr);
  // Initialization of ,echanical and electrical torques
  parameter Real Te0 = (Edp0 * Id0 + Eqp0 * Iq0) / (w0 * omega_0);
  parameter Real T0 = Te0 / (A * (w0 * omega_0) ^ 2 + B * (w0 * omega_0) + C) "Initial mechanical torque";
  // Initialization of d-q electrical values
  parameter Real Vd0 = v_0 * cos(angle_0);
  parameter Real Vq0 = v_0 * sin(angle_0);
  parameter Real Id0 = (P_0 * S_b / M_b * Vd0 + Qmotor * Vq0) / (Vd0 ^ 2 + Vq0 ^ 2);
  parameter Real Iq0 = (P_0 * S_b / M_b * Vq0 - Qmotor * Vd0) / (Vd0 ^ 2 + Vq0 ^ 2);
  parameter Real Edp0 = Vd0 - Rs * Id0 + Xp * Iq0;
  parameter Real Eqp0 = Vq0 - Rs * Iq0 - Xp * Id0;
  // Computing the initial slip, speed and reactive power from the steady-state equivalent circuit
  parameter Real K1 = Xs * Xm + Xs * Xr + Xm * Xr;
  parameter Real K2 = Xs + Xm;
  parameter Real K3 = Xm + Xr;
  parameter Real a = P_0 * S_b / M_b * K1 ^ 2 + P_0 * S_b / M_b * Rs ^ 2 * K3 ^ 2 - Rs * K3 ^ 2 * v_0 ^ 2;
  parameter Real b = 2 * P_0 * S_b / M_b * Rs * Rr * (K2 * K3 - K1) - Rr * v_0 ^ 2 * (K2 * K3 - K1);
  parameter Real c = Rr ^ 2 * (P_0 * S_b / M_b * Rs ^ 2 + P_0 * S_b / M_b * K2 ^ 2 - Rs * v_0 ^ 2);
  parameter Real s0 = 1 / (2 * a) * ((-b) - sqrt(b ^ 2 - 4 * a * c)) "Initial slip";
  parameter Real w0 = 1 - s0 "Initial speed";
  parameter Real Qmotor = v_0 ^ 2 * (Rr ^ 2 * K2 + s0 ^ 2 * K3 * K1) / ((Rs * Rr - s0 * K1) ^ 2 + (Rr * K2 + Rs * s0 * K3) ^ 2);
initial equation
  w = w0;
  wr = w0 * omega_0;
  Edp = Edp0;
  Eqp = Eqp0;
  delta = angle_0;
equation
// Integration of frequency differences to determine delta for d-q transformation
  der(delta) = 2 * Modelica.Constants.pi * 50 * (win - 1);
// Transformation from 50 Hz rotating frame to the frame rotating with der(phiu)
// i.e. with the frequency of the reference machine
  [Vd; Vq] = [cos(delta), sin(delta); -sin(delta), cos(delta)] * [p.vr; p.vi];
  [Id; Iq] = S_b / M_b * [cos(delta), sin(delta); -sin(delta), cos(delta)] * [p.ir; p.ii];
// Mechanical algebraic equations
  wr = w * win;
  Tm = (A * wr ^ 2 + B * wr + C) * T0;
  Te = (Edp * Id + Eqp * Iq) / wr;
// Induction machine - differential equations
  der(wr) = -1 / (2 * H) * (Tm - Te);
  der(Edp) = (-1 / Tp * (Edp + (X - Xp) * Iq)) - (w - 1) * Eqp;
  der(Eqp) = (-1 / Tp * (Eqp - (X - Xp) * Id)) + (w - 1) * Edp;
// Induction machine - algebraic equations
  Vd - Edp = Id * Rs - Iq * Xp;
  Vq - Eqp = Iq * Rs + Id * Xp;
  P = p.vr * p.ir + p.vi * p.ii;
  Q = (-p.vr * p.ii) + p.vi * p.ir;
  annotation(
    __OpenModelica_simulationFlags(lv = "LOG_RES_INIT,LOG_STATS", outputFormat = "mat", s = "dassl"),
    Icon(graphics = {Ellipse(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}, endAngle = 360), Text(origin = {-42, 66}, extent = {{-40, -6}, {122, -126}}, textString = "AM")}));
end AsyncMachine;