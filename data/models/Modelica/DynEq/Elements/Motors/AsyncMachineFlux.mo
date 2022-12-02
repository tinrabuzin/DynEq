within DynEq.Elements.Motors;

model AsyncMachineFlux
  import SI = Modelica.SIunits;
  import Modelica.Constants.pi;
  // Initialization of the machine as in
  //J. Powell, G. Radman, 'Initialization for Dynamic Simulation of Stressed Power Systems Considering Induction Motor Components of Loads'
  // Power flow results
  parameter Real P_0 "Initial active power";
  parameter Real Q_0 "Initial reactive power";
  parameter Real v_0 "Initial voltage";
  parameter Real angle_0 "Initial angle";
  parameter Real omega_0 "Initial frequency";
  parameter Real S_b;
  parameter Real fn = 50;
  // Parameters of the induction machine
  parameter Real Mlf = 0.8 annotation(
    Dialog(group = "Machine Parameters"));
  parameter Real Rs = 0.201 annotation(
    Dialog(group = "Machine Parameters"));
  parameter Real Xs = 0.1791 annotation(
    Dialog(group = "Machine Parameters"));
  parameter Real Xm = 2.0141 annotation(
    Dialog(group = "Machine Parameters"));
  parameter Real Rr = 0.0424 annotation(
    Dialog(group = "Machine Parameters"));
  parameter Real Xr = 0.1507 annotation(
    Dialog(group = "Machine Parameters"));
  parameter Real H = 1.4529 annotation(
    Dialog(group = "Machine Parameters"));
  parameter Real A = 0.6267 annotation(
    Dialog(group = "Machine Parameters"));
  parameter Real B = 0.4471 annotation(
    Dialog(group = "Machine Parameters"));
  parameter Real C = 1 - A - B;
  // Machine speeds
  Real wr;
  // w.r.t the frame rotating with 50Hz, i..e actual mechanical speed
  // d-q values of the machine model
  Real psi_rd;
  Real psi_rq;
  Real psi_sd;
  Real psi_sq;
  Real i_sd;
  Real i_sq;
  Real u_sd;
  Real u_sq;
  // Real and imaginary quantities
  Real Irsh;
  Real Iish;
  Real Irm;
  Real Iim;
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
  parameter Real T0 = Te0 / (A * (w0 * omega_0) ^ 2 + B * (w0 * omega_0) + C) "Initial mechanical torque";
  parameter Real Te0 = (-Id0 * psi_sq0) + Iq0 * psi_sd0;
  parameter Real x_ss = Xs + Xm;
  parameter Real x_sr = Xm;
  parameter Real x_rr = Xr + Xm;
  parameter Real r_r = Rr;
  parameter Real r_s = Rs;
  parameter Real xhat = x_ss - x_sr / x_rr * x_sr;
  //
  // Initialization of d-q electrical values
  parameter Real Vd0 = v_0 * cos(angle_0);
  parameter Real Vq0 = v_0 * sin(angle_0);
  parameter Real Id0 = (P_0 * S_b / M_b * Vd0 + Qmotor * Vq0) / (Vd0 ^ 2 + Vq0 ^ 2);
  parameter Real Iq0 = (P_0 * S_b / M_b * Vq0 - Qmotor * Vd0) / (Vd0 ^ 2 + Vq0 ^ 2);
  parameter Real a1 = r_r / x_rr * x_sr;
  parameter Real b1 = r_r / x_rr;
  parameter Real c1 = (2 * pi * fn * omega_0 - 2 * pi * fn * w0 * omega_0) / 2 / pi / fn;
  parameter Real psi_rd0 = (a1 * b1 * Id0 + a1 * c1 * Iq0) / (b1 ^ 2 + c1 ^ 2);
  parameter Real psi_rq0 = (a1 * b1 * Iq0 - a1 * c1 * Id0) / (b1 ^ 2 + c1 ^ 2);
  parameter Real psi_sq0 = (r_s * Id0 - Vd0) / omega_0;
  parameter Real psi_sd0 = (Vq0 - r_s * Iq0) / omega_0;
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
  parameter Real Qcomp = Q_0 - Qmotor * M_b / S_b;
  parameter Real Bsh = -Qcomp / v_0 ^ 2;
initial equation
  psi_rd = (psi_sd - xhat * Id0) * x_rr / x_sr;
  psi_rq = (psi_sq - xhat * Iq0) * x_rr / x_sr;
  wr = w0 * omega_0;
  delta = angle_0;
//psi_rd = psi_rd0;
//psi_rq = psi_rq0;
equation
  der(psi_rd) / (2 * pi * fn) = a1 * i_sd - b1 * psi_rd + (win - wr) * psi_rq;
  der(psi_rq) / (2 * pi * fn) = a1 * i_sq - b1 * psi_rq - (win - wr) * psi_rd;
  psi_sd = xhat * i_sd + x_sr / x_rr * psi_rd;
  psi_sq = xhat * i_sq + x_sr / x_rr * psi_rq;
  u_sd = r_s * i_sd - win / 1 * psi_sq;
  u_sq = r_s * i_sq + win / 1 * psi_sd;
// Mechanical equations
  Tm = (A * wr ^ 2 + B * wr + C) * T0;
  Te = (-i_sd * psi_sq) + i_sq * psi_sd;
  der(wr) = -1 / (2 * H) * (Tm - Te);
// Integration of frequency differences to determine delta for d-q transformation
  der(delta) = 2 * Modelica.Constants.pi * 50 * (win - 1);
// Transformation from 50 Hz rotating frame to the frame rotating with der(phiu)
// i.e. with the frequency of the reference machine
  [u_sd; u_sq] = [cos(delta), sin(delta); -sin(delta), cos(delta)] * [p.vr; p.vi];
  [i_sd; i_sq] = S_b / M_b * [cos(delta), sin(delta); -sin(delta), cos(delta)] * [Irm; Iim];
// Shunt compensation
  p.vr * Irsh + p.vi * Iish = 0;
  (p.vr ^ 2 + p.vi ^ 2) * (-Bsh) = (-p.vr * Iish) + p.vi * Irsh;
  Irm + Irsh = p.ir;
  Iim + Iish = p.ii;
  P = p.vr * p.ir + p.vi * p.ii;
  Q = (-p.vr * p.ii) + p.vi * p.ir;
  annotation(
    __OpenModelica_simulationFlags(lv = "LOG_RES_INIT,LOG_STATS", outputFormat = "mat", s = "dassl"),
    Icon(graphics = {Ellipse(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}, endAngle = 360), Text(origin = {-42, 66}, extent = {{-40, -6}, {122, -126}}, textString = "AM")}));
end AsyncMachineFlux;