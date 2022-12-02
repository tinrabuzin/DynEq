within DynEq.Elements;

model CompositeLoad
  extends DynEq.Elements.BaseClasses.OnePort(
  p(i(re(nominal=M_b/VBase), im(nominal=M_b/VBase)))
  );
  parameter SI.PerUnit omega_0 "Initial frequency" annotation (Dialog(group="Power flow data"));
  // General Parameters
  parameter Real Kpm = 0.8 annotation(
    Dialog(group = "General Parameters"));
  parameter Real Mlf = 0.8 annotation(
    Dialog(group = "General Parameters"));
  // Parmeters of the exp model
  parameter Real expP = 2 annotation(
    Dialog(group = "Exponential Load Model Parameters"));
  parameter Real expQ = 2 annotation(
    Dialog(group = "Exponential Load Model Parameters"));
  // Parameters of the induction machine
  parameter SI.PerUnit Rs = 0.028 annotation(
    Dialog(group = "Machine Parameters"));
  parameter SI.PerUnit Xs = 0.01 annotation(
    Dialog(group = "Machine Parameters"));
  parameter SI.PerUnit Xm = 4.236 annotation(
    Dialog(group = "Machine Parameters"));
  parameter SI.PerUnit Rr = 0.00489 annotation(
    Dialog(group = "Machine Parameters"));
  parameter SI.PerUnit Xr = 0.1822 annotation(
    Dialog(group = "Machine Parameters"));
  parameter Real H = 1.67 annotation(
    Dialog(group = "Machine Parameters"));
  parameter Real expT = 2 annotation(
    Dialog(group = "Machine Parameters"));

  final parameter Real omega_base = 2 * C.pi * fBase;
  // Active and reactive powers
  final parameter Types.ActivePower Pm_0 = Kpm * P_0;
  final parameter Types.ApparentPower M_b = Pm_0 / Mlf;
  final parameter Types.ActivePower Pexp_0 = P_0 - Pm_0;
  final parameter Real x_ss = Xs + Xm;
  final parameter Real x_sr = Xm;
  final parameter Real x_rr = Xr + Xm;
  final parameter Real r_r = Rr;
  final parameter Real r_s = Rs;
  final parameter Real xhat = x_ss - x_sr / x_rr * x_sr;

// Start parameters
  parameter Types.ComplexCurrent I_m0 = Complex(i_sd0, i_sq0)*M_b/VBase;
  parameter Types.ComplexCurrent I_e0 = I_0 - I_m0;
  final parameter SI.PerUnit u_sd0 = v_0 * cos(angle_0) / VBase;
  final parameter SI.PerUnit u_sq0 = v_0 * sin(angle_0) / VBase;
  final parameter SI.PerUnit psi_rd0(fixed=false);
  final parameter SI.PerUnit psi_rq0(fixed=false);
  final parameter SI.PerUnit psi_sq0(fixed=false);
  final parameter SI.PerUnit psi_sd0(fixed=false);
  final parameter SI.PerUnit i_sd0Start = I_0.re/(M_b/VBase);
  final parameter SI.PerUnit i_sq0Start = I_0.im/(M_b/VBase);
  final parameter SI.PerUnit i_sd0(start=i_sd0Start, fixed=false);
  final parameter SI.PerUnit i_sq0(start=i_sq0Start, fixed=false);
  final parameter SI.PerUnit w_0(start=omega_0, fixed=false);
  final parameter Types.ReactivePower Qexp_0 = Q_0;
  final parameter Types.ReactivePower Qm0(fixed=false);

  // Model variables
  // Angle w.r.t. the reference angle in the 50Hz rotating frame
  SI.PerUnit delta;
  // Synchronous frequency
  Modelica.Blocks.Interfaces.RealInput win annotation(
    Placement(visible = true, transformation(origin = {-100, 70}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, 0},  extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Types.ComplexCurrent imach(re(start=I_m0.re), im(start=I_m0.im));
  Types.ComplexCurrent iexp(re(start=I_e0.re), im(start=I_e0.im));
  SI.PerUnit wr(start=w_0);
  SI.PerUnit psi_rd(start=psi_rd0);
  SI.PerUnit psi_rq(start=psi_rq0);
  SI.PerUnit psi_sd(start=psi_sd0);
  SI.PerUnit psi_sq(start=psi_sq0);
  SI.PerUnit i_sd(start=i_sd0);
  SI.PerUnit i_sq(start=i_sq0);
  SI.PerUnit u_sd(start=u_sd0);
  SI.PerUnit u_sq(start=u_sq0);
  SI.PerUnit T_m(start=-i_sd0 * psi_sq0 + i_sq0 * psi_sd0);
  SI.PerUnit T_e(start=-i_sd0 * psi_sq0 + i_sq0 * psi_sd0);
  final parameter SI.PerUnit T_0(fixed=false, start=-i_sd0 * psi_sq0 + i_sq0 * psi_sd0);
initial equation
// Start values for the machine model
  0 = r_r / x_rr * x_sr * i_sd0 - r_r / x_rr * psi_rd0 + (omega_0 - w_0) * psi_rq0;
  0 = r_r / x_rr * x_sr * i_sq0 - r_r / x_rr * psi_rq0 - (omega_0 - w_0) * psi_rd0;
  psi_sd0 = xhat * i_sd0 + x_sr / x_rr * psi_rd0;
  psi_sq0 = xhat * i_sq0 + x_sr / x_rr * psi_rq0;
  u_sd0 = r_s * i_sd0 - omega_0  * psi_sq0;
  u_sq0 = r_s * i_sq0 + omega_0  * psi_sd0;
  Pm_0  = (u_sd0 * i_sd0 + u_sq0 * i_sq0) * M_b;
// Initialization
  Qm0 = (u_sq0 * i_sd0 - u_sd0 * i_sq0) * M_b;
  P = P_0;
  delta = 0;
  der(wr) = 0;
  der(psi_rd) = 0;
  der(psi_rq) = 0;
equation
// Integration of frequency differences to determine delta for d-q transformation
//  assert(w_0 < omega_0, "Motor failed to initialize to a stable point");
  der(delta) = omega_base * (win - 1);
// d-q transformations
  [u_sd; u_sq] * VBase = [cos(delta), sin(delta); -sin(delta), cos(delta)] * [p.v.re; p.v.im];
  [i_sd; i_sq] * (M_b/VBase) = [cos(delta), sin(delta); -sin(delta), cos(delta)] * [imach.re; imach.im];
// Mechanical equations
  T_m = T_0 * wr ^ expT;
  T_e = -i_sd * psi_sq + i_sq * psi_sd;
  der(wr) = 1 / (2 * H) * (T_e - T_m);
// Electrical equations
  u_sd = r_s * i_sd - win * psi_sq;
  u_sq = r_s * i_sq + win * psi_sd;
  psi_sd = xhat * i_sd + x_sr / x_rr * psi_rd;
  psi_sq = xhat * i_sq + x_sr / x_rr * psi_rq;
  der(psi_rd) / omega_base = r_r / x_rr * x_sr * i_sd - r_r / x_rr * psi_rd + (win - wr) * psi_rq;
  der(psi_rq) / omega_base = r_r / x_rr * x_sr * i_sq - r_r / x_rr * psi_rq - (win - wr) * psi_rd;
// Exponential load model
  p.v * CM.conj(iexp) = Complex(Pexp_0 * (CM.'abs'(p.v) / v_0) ^ expP, -Qm0 + Qexp_0 * (CM.'abs'(p.v) / v_0) ^ expQ);
  p.i = imach + iexp;
  annotation(
    Icon(graphics = {Ellipse(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}, endAngle = 360), Text(origin = {-40.989, -17}, extent = {{-39.011, -3}, {118.989, -63}}, textString = "COMPOSITE")}, coordinateSystem(initialScale = 0.1)));
end CompositeLoad;
