within Flex;
model FlexibleThinBeam "Flexible thin beam model"
  import      Modelica.Units.SI;
  import Cv = Modelica.Units.Conversions;
  import Modelica.Mechanics.MultiBody.Frames;

  parameter SI.Density rho=7800 "Material density" annotation(Evaluate = true);
  parameter SI.Length L "Beam length" annotation(Evaluate = true);

  parameter Boolean clampedFree=true
    "Clamped-Free model if true, else simply-supported model"
    annotation(Dialog(group="Boundary Conditions"));
  parameter Boolean circularSection=true
    "CircularSection if true, else Rectangular"
    annotation(Dialog(group="3D Graphics"));
  parameter Modelica.Mechanics.MultiBody.Types.Color beamColor={128,128,128}
    "Beam color "
    annotation(Dialog(group="3D Graphics"));
  parameter SI.Area A "Cross sectional area" annotation(Evaluate = true);
  parameter SI.ModulusOfElasticity E "Young modulus of the beam material" annotation(Evaluate = true) ;
  parameter SI.SecondMomentOfArea J "Second moment of area of cross section around z-axis"  annotation (Evaluate = true);
  parameter Integer N(min=1) = 5 "Number of elements";
  parameter SI.PerUnit csi1 = 0 "Desired damping coefficient at frequency omega1"
    annotation(Dialog(group="Damping"));
  parameter SI.PerUnit csi2 = 0 "Desired damping coefficient at frequency omega2"
      annotation(Dialog(group="Damping"));
  parameter SI.AngularFrequency omega1 = 10 "Angular frequency 1"
    annotation(Dialog(group="Damping"));
  parameter SI.AngularFrequency omega2 = 100 "Angular frequency 2"
    annotation(Dialog(group="Damping"));
  final parameter SI.DampingCoefficient alpha(fixed = false)
    "Rayleigh structural damping proportional to mass";
  final parameter SI.Time beta(fixed = false)
    "Rayleigh structural damping proportional to stiffness";
  final parameter SI.Length h=L/N "Length of one element";
  final parameter SI.Mass m=rho*L*A "Mass of the beam";
protected
  constant Real pi=Modelica.Constants.pi;
  final parameter Real SbarEl[3, 6]=
    m/N/12*[6, 0, 0, 6, 0,  0;
            0, 6, h, 0, 6, -h;
            0, 0, 0, 0, 0,  0];
  final parameter Real Sbar11[6, 6]=
    m/N*[1/3, 0, 0, 1/6, 0, 0;
           0, 0, 0,   0, 0, 0;
           0, 0, 0,   0, 0, 0;
         1/6, 0, 0, 1/3, 0, 0;
           0, 0, 0,   0, 0, 0;
           0, 0, 0,   0, 0, 0];
  final parameter Real Sbar12[6, 6]=transpose(Sbar21);
  final parameter Real Sbar13[6, 6]=zeros(6, 6);
  final parameter Real Sbar21[6, 6]=
    m/N*[   0, 0, 0,      0, 0, 0;
         7/20, 0, 0,   3/20, 0, 0;
       1/20*h, 0, 0, 1/30*h, 0, 0;
            0, 0, 0,      0, 0, 0;
         3/20, 0, 0,   7/20, 0, 0;
      -1/30*h, 0, 0,-1/20*h, 0, 0];
  final parameter Real Sbar22[6, 6]=
    m/N*[0,         0,          0, 0,         0,         0;
         0,     13/35,   11/210*h, 0,      9/70, -13/420*h;
         0,  11/210*h,  1/105*h^2, 0,  13/420*h,-1/140*h^2;
         0,         0,          0, 0,         0,         0;
         0,      9/70,   13/420*h, 0,     13/35, -11/210*h;
         0, -13/420*h, -1/140*h^2, 0, -11/210*h, 1/105*h^2];
  final parameter Real Sbar23[6,6]=zeros(6, 6);
  final parameter Real Sbar31[6,6]=zeros(6, 6);
  final parameter Real Sbar32[6,6]=zeros(6, 6);
  final parameter Real Sbar33[6,6]=zeros(6, 6);
  final parameter Real Ibar11[1,6] =  h*m/N*[1/6, 0, 0, 1/3, 0, 0];
  final parameter Real Ibar11adj[1,6]=h*m/N*[1/2, 0, 0, 1/2, 0, 0];
  final parameter Real Ibar12[1,6] =  h*m/N*[0, 3/20, 1/30*h, 0, 7/20, -1/20*h];
  final parameter Real Ibar12adj[1,6]=h*m/N*[0,  1/2, 1/12*h, 0,  1/2, -1/12*h];
  final parameter Real Ibar13[1,6]=zeros(1, 6);
  final parameter Real Ibar21[1,6]=zeros(1, 6);
  final parameter Real Ibar22[1,6]=zeros(1, 6);
  final parameter Real Ibar23[1,6]=zeros(1, 6);
  final parameter Real Ibar31[1,6]=zeros(1, 6);
  final parameter Real Ibar32[1,6]=zeros(1, 6);
  final parameter Real Ibar33[1,6]=zeros(1, 6);
  final parameter Real KffEl[6, 6]=
    E*[A/h,  0,         0,       -A/h,  0,         0;
       0,    12*J/h^3,  6*J/h^2,  0,   -12*J/h^3,  6*J/h^2;
       0,    6*J/h^2,   4*J/h,    0,   -6*J/h^2,   2*J/h;
      -A/h,  0,         0,        A/h,  0,         0;
       0,   -12*J/h^3, -6*J/h^2,  0,    12*J/h^3, -6*J/h^2;
       0,    6*J/h^2,   2*J/h,    0,   -6*J/h^2,   4*J/h];
  final parameter Real S0[3, 6]=
    [1, 0, 0, 0, 0, 0;
     0, 1, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0];
  final parameter Real dS0[3, 6]=
    [0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0];
  final parameter Real S1[3, 6]=
    [0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 0, 0];
  final parameter Real S2[3, 6]=
    [0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 0, 0];
  final parameter Real dS1[3, 6]=
    [0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 1];
  type vec3D = SI.Position[3];

public
  Real qf[3*N](start = zeros(3*N)) "Elastic coordinates"
    annotation (Dialog(tab="Initialization",showStartAttribute=true));
  Real dqf[3*N](start = zeros(3*N)) "Elastic velocities"
    annotation (Dialog(tab="Initialization",showStartAttribute=true));
  Real ddqf[3*N](start = zeros(3*N)) "Elastic accelerations"
    annotation (Dialog(tab="Initialization",showStartAttribute=true));
  SI.Acceleration g_0[3] "Gravity acceleration resolved in world frame";

  /* 3D graphics variables */
  vec3D r0shape[N] "Position of left boundary of finite elements resolved in frame_a";
  vec3D rrelshape[N] "Left boundary position - Right boundary position vector resolved in frame_a";
  Real Lshape[N] "Length of each finite volume";

protected
  Real B[N,6,3*N];
  Real Stbar[3];
  Real StbarCross[3, 3];
  Real Ithth_bar[3, 3];
  Real Ithth_bar11;
  Real Ithth_bar22;
  Real Ithth_bar33;
  Real Ithth_bar12;
  Real Ithth_bar_der[3, 3];
  Real Ithth_bar11_der;
  Real Ithth_bar22_der;
  Real Ithth_bar33_der;
  Real Ithth_bar12_der;
  Real QvR[3, 1];
  Real QvAlpha[3, 1];
  Real wCross[3, 3];
  Real wCross2[3, 3];
  Real Sbar[3, 3*N];
  Real Ithf_bar[3, 3*N];
  Real mff[3*N, 3*N];
  Real Kff[3*N, 3*N];
  Real F[6, 6];
  Real H[6, 6];
  Real Qef[3*N, 1];
  Real Qvf[3*N, 1];


  Real ra[3];
  Real va[3];
  Real v_0a[3];
  Real aa[3];
  Real wa[3];
  Real za[3];
  Real fa[3];
  Real ta[3];

  Real rb[3];
  Real fb[3];
  Real tb[3];
  Real fb_a[3];
  Real tb_a[3];

public
  Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a "Coordinate system fixed to the component with one cut-force and cut-torque" annotation (Placement(
        transformation(extent={{-106,-8},{-86,12}}, rotation=0),
        iconTransformation(extent={{-108,-10},{-86,12}})));
  Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b "Coordinate system fixed to the component with one cut-force and cut-torque" annotation (Placement(
        transformation(extent={{88,-8},{108,12}}, rotation=0),
        iconTransformation(extent={{86,-10},{108,12}})));

protected
  Modelica.Mechanics.MultiBody.Frames.Orientation R_rel;
  outer Modelica.Mechanics.MultiBody.World world;
  Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape Segment[N](
    r_shape=r0shape,
    length=Lshape,
    each width=if circularSection then 2*sqrt(A/pi) else A/(0.4*sqrt(A)),
    each height=if circularSection then 2*sqrt(A/pi) else (0.4*sqrt(A)),
    lengthDirection=rrelshape,
    each widthDirection={0,0,1},
    each shapeType=if circularSection then "cylinder" else "box",
    each color=beamColor,
    each extra=0.0,
    each r=frame_a.r_0,
    each R=frame_a.R);

initial equation
   // initial equations to determine the Rayleigh damping coefficients alpha and beta
   csi1 = 1/2*(alpha/omega1 + beta*omega1);
   csi2 = 1/2*(alpha/omega2 + beta*omega2);
equation
  Connections.branch(frame_a.R, frame_b.R);
  assert(cardinality(frame_a) > 0 or cardinality(frame_b) > 0, "Neither connector frame_a nor frame_b of FlexBeamFem object is connected");

  //connectivity matrices
  if clampedFree then
    B[1,:,:] = [zeros(3, 3*N); identity(3), zeros(3, 3*(N - 1))];
    B[N,:,:] = [zeros(6, 3*(N - 2)), identity(6)];
  else
    B[1,:,:] = [zeros(2, 3*N);
                [zeros(1, 3*(N - 1)), [0, 1, 0]];
                identity(3), zeros(3, 3*(N - 1))];
    B[N,:,:] = [zeros(3, 3*(N - 2)), identity(3), zeros(3, 3);
                zeros(3, 3*(N - 2)), zeros(3, 3), [1, 0, 0; 0, 0, 0; 0, 0, 1]];
  end if;

  for i in 2:N - 1 loop
    B[i,:,:] = [zeros(6, 3*(i - 2)), identity(6), zeros(6, 3*(N - i))];
  end for;

  //integral shape function matrices
  Sbar = sum(SbarEl*B[i,:,:] for i in 1:N);

  //rotation-deformation coupling
  Ithf_bar =
    [zeros(2, 3*N);
     sum((Ibar12 + (i - 1)*Ibar12adj +
      transpose(matrix(qf))*transpose(B[i,:,:])*(Sbar12 - Sbar21))*B[i,:,:]
     for i in 1:N)];

  //mass structural matrix
  mff = sum(transpose(B[i,:,:])*(Sbar11 + Sbar22)*B[i,:,:] for i in 1:N);

  //Inertia matrix and derivative
  Ithth_bar11 =1e-10+
   scalar(transpose(matrix(qf))*sum((transpose(B[i,:,:])*Sbar22*B[i,:,:]) for i in 1:N)*qf);

  Ithth_bar22 = m*L^2/3 +
    scalar(2*sum((Ibar11 + (i - 1)*Ibar11adj)*B[i,:,:] for i in 1:N)*qf +
    transpose(matrix(qf))*sum(transpose(B[i, :, :])*Sbar11*B[i,:,:] for i in 1:N)*qf);

  Ithth_bar33 = Ithth_bar11 + Ithth_bar22;

  Ithth_bar12 = scalar(-sum((Ibar12 + (i - 1)*Ibar12adj)*B[i, :, :] for i in 1:
    N)*qf - transpose(matrix(qf))*sum(transpose(B[i, :, :])*Sbar21*B[i, :, :]
    for i in 1:N)*qf);

  Ithth_bar = [Ithth_bar11, Ithth_bar12, 0;
               Ithth_bar12, Ithth_bar22, 0;
               0,           0,    Ithth_bar33];

  Ithth_bar11_der =
    2*scalar(transpose(matrix(qf))*sum((transpose(B[i,:,:])*Sbar22*B[i,:,:]) for i in 1:N)*dqf);

  Ithth_bar22_der =
    scalar(2*sum((Ibar11 + (i - 1)*Ibar11adj)*B[i,:,:] for i in 1:N)*dqf +
    2*transpose(matrix(qf))*sum(transpose(B[i,:,:])*Sbar11*B[i,:,:] for i in 1:N)*dqf);

  Ithth_bar12_der =
    scalar(-sum((Ibar12 + (i - 1)*Ibar12adj)*B[i,:,:] for i in 1:N)*dqf -
           2*transpose(matrix(qf))*sum(transpose(B[i, :, :])*Sbar21*B[i,:,:] for i in 1:N)*dqf);

  Ithth_bar33_der = Ithth_bar11_der + Ithth_bar22_der;

  Ithth_bar_der = [Ithth_bar11_der, Ithth_bar12_der,    0;
                   Ithth_bar12_der, Ithth_bar22_der,    0;
                   0,               0,           Ithth_bar33_der];

  //stiffness matrix
  Kff = sum(transpose(B[i,:,:])*(KffEl)*B[i,:,:] for i in 1:N);

  /* Flange A quantities definitions */
  g_0 =Frames.resolve2(frame_a.R, world.gravityAcceleration(frame_a.r_0 +
    Frames.resolve1(frame_a.R, {L/2,0,0})));

  ra =frame_a.r_0;
  va =Frames.resolve2(frame_a.R, der(frame_a.r_0));
  wa =Frames.angularVelocity2(frame_a.R);
  v_0a = der(ra);
  aa =Frames.resolve2(frame_a.R, der(v_0a));
  za = der(wa);
  fa =frame_a.f;
  ta =frame_a.t;

  /* Time-variant quantities */
  dqf = der(qf);
  ddqf = der(dqf);

  QvR = matrix(-cross(wa, cross(wa, Stbar)) - 2*cross(wa, (Sbar*dqf)));
  QvAlpha = matrix(-cross(wa, Ithth_bar*wa) - der(Ithth_bar)*wa - cross(wa, (Ithf_bar*dqf)));
  Stbar = {m*L/2,0,0} + Sbar*qf;
  StbarCross = skew(Stbar);

  wCross = skew(wa);
  wCross2 = wCross*wCross;

  F = Sbar11*wCross2[1, 1] + Sbar12*wCross2[2, 1] +
      Sbar21*wCross2[1, 2] + Sbar22*wCross2[2, 2];

  H = -(Sbar12*wCross[2, 1] + Sbar21*wCross[1, 2]);

  Qvf = -sum((transpose(B[i, :, :])*(transpose([(Ibar11 + (i - 1)*Ibar11adj)*
    wCross2[1, 1] + (Ibar12 + (i - 1)*Ibar12adj)*wCross2[2, 1]]))) + matrix((
    transpose(B[i,:,:])*F*B[i,:,:])*qf) + matrix(2*(transpose(B[i,:,:])*H
    *B[i,:,:])*dqf) for i in 1:N);

  transpose(Qef) = (transpose(matrix(fb_a))*S1*B[N,:,:] +
                    transpose(matrix(tb_a))*dS1*B[N,:,:] +
                    transpose(matrix(ta))*dS0*B[1,:,:]);

  /* Dynamics equations */
  [m*identity(3), transpose(StbarCross), Sbar]*[aa - g_0; za; ddqf] =
    QvR + matrix(fa + fb_a);

  [StbarCross, Ithth_bar, Ithf_bar]*[aa - g_0; za; ddqf] =
    QvAlpha + matrix(ta + tb_a + cross(({L,0,0} + S1*B[N,:,:]*qf), fb_a));

  [transpose(Sbar), transpose(Ithf_bar), mff]*[aa - g_0; za; ddqf] =
    Qvf + Qef - matrix(Kff*qf) - matrix((alpha*mff + beta*Kff)*dqf);

  /* Flange B quantities definitions */
  rb =frame_b.r_0;
  fb =frame_b.f;
  tb =frame_b.t;
  fb_a = Frames.resolve1(R_rel, fb);
  tb_a = Frames.resolve1(R_rel, tb);
  rb =ra + Frames.resolve1(frame_a.R, ({L,0,0} + S1*B[N, :, :]*qf));
  frame_b.R = Frames.absoluteRotation(frame_a.R, R_rel);
  R_rel = Frames.planarRotation({0,0,1}, qf[3*N], dqf[3*N]);

  /* 3D Visual Representation */
  r0shape[1, :] = {0,0,0};
  rrelshape[1, :] = r0shape[2, :];
  Lshape[1] = sqrt(rrelshape[1, :]*rrelshape[1, :]);

  r0shape[N, :] = {L*(N - 1)/N,0,0} + S1*B[N - 1, :, :]*qf;
  rrelshape[N, :] = {L,0,0} + S1*B[N, :, :]*qf - r0shape[N, :];
  Lshape[N] = sqrt(rrelshape[N, :]*rrelshape[N, :]);

  for i in 2:N - 1 loop
    r0shape[i, :] = {L*(i - 1)/N,0,0} + S1*B[i - 1, :, :]*qf;
    rrelshape[i, :] = (r0shape[i + 1, :] - r0shape[i, :]);
    Lshape[i] = sqrt(rrelshape[i, :]*rrelshape[i, :]);
  end for;
  annotation (
    Icon(graphics={
        Text(
          extent={{-54,-30},{50,-52}},
          lineColor={0,0,255},
          textString=
               "%name"),
        Rectangle(
          extent={{-94,20},{88,-20}},
          lineColor={135,135,135},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{88,20},{94,-20}},
          lineColor={135,135,135},
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid)}),
    uses(Modelica(version="3.2.2")),
    Documentation(info="<html>
<p>Flexible beam model compatible with the MultiBody library. The beam model is based on
<a href=\"https://en.wikipedia.org/wiki/Euler%E2%80%93Bernoulli_beam_theory\">Euler-Bernoulli theory</a>, 
accounting for planar deformation on the x-y plane, including compression and bending, and neglecting the effect
of shear stress in the formulation of the dynamic equations. The beam is infinitely rigid with respect to 
out-of-plane deformation and torsional deformation. It is thus suitable to model thin beams subject
to forces and torques that mainly cause x-y plane deformation.</p>
<p>The required data are the material properties, the beam length and cros section, and the second moment of area
<code>J</code> of the beam cross section around the z-axis. For a rectangular beam with thickness <code>h</code> (in the y-axis
direction) and width <code>h</code> (in the z-axis direction), <code>J = b*h^3/12</code>. Here is a 
<a href=\"https://en.wikipedia.org/wiki/List_of_second_moments_of_area\">list of moments of area for other shapes</a>.</p>
<p>The Euler-Bernoulli PDEs are discretized into <code>N</code> elements by the Finite Element Method, 
assuming small deformations in the x-y plane. The beam position is thus determined by the position and rotation of the reference
frame <code>frame_a</code>, to which a small deformation described by the elastic coordinates <code>qf</code>, by means of suitable
shape functions, is superimposed. More specifically, the left boundary of the first element is rigidly attached to <code>frame_a</code>,
and then the right boundary of each i-th element has a displacement in the x-direction proportional to <code>q[1+(i-1)*3]</code>,
a displacement in the y-direction proportional to <code>q[(2+(i-1)*3]</code>, and a rotation around the z-axis proportional
to <code>q[(3+(i-1)*3]</code>.</p>
<p>The position of the left boundary of each element, resolved in the <code>frame_a</code> reference, is computed in the vector 
<code>r0shape</code>, which is also used for visualization.</p>
<p>The shape of the basis functions that describe the beam deformation can be chosen between the clamped-free and the simply supported
configuration. In general, it is possible to connect the two frame connectors at the beam boundary to any other rigid or flexible body connector;
the dynamic behaviour due to flexibility will be better approximated by increasing the number of elements <code>N</code>. However,
an appropriate choice of the basis shape function allows to obtain a better approximation of the exact motion with a lower number of
elements.</p>
<p>The application of the FEM method to the Euler-Bernoulli PDEs allows to compute the mass and stiffness matrices of the model, which has
no inherent damping. Some damping may be introduced by other components connected to the beam, though this is usually not enough to
inject sufficient damping on high-frequency oscillation modes. The presence of high-frequency poorly damped or undamped modes can have a
catastrophic effect on simulation time when variable step-size algorithms with error control are used, because they may require a very high
number of simulation steps to follow the fast undamped oscillations with the required accuracy. This unwanted phenomenon can be avoided
by means of two strategies.</p>
<p>The first strategy is to employ integration algorithms which are A-stable, or at least such that the stability boundary remains parallel
to the imaginary axis for the widest possible frequency range. This ensures that poorly damped modes are not further undamped by the
integration algoritm, ensuring that their oscillations die out after a short enough transients. Implicit Runge-Kutta algorithms such as Radau II
are thus recommended, while BDF algorithms such as DASSL or IDA may not work well, because of their a marked tendency to destabilize
poorly damped modes at frequencies between one and 15 times the inverse of the time step.</p>
<p>The second strategy is to introduce some structural damping in the model.
Unfortunately, this can only be done empirically, by fitting the damping coefficient of some vibration modes of the system. The flexible
beam model allows to do so using Rayleigh damping, whereby the damping matrix is a linear combination of the mass and stiffness matrices
with coefficients <code>alpha</code> and <code>beta</code>. The coefficients are computed to obtain the desired damping coefficients 
<code>csi1</code> and <code>csi2</code> at the angular frequencies <code>omega1</code> and <code>omega2</code>; 
damping will be higher outside this interval, and slighly lower within this interval.
The recommended way to use this feature is to first run a simulation with the default zero damping coefficients, 
identify the angular frequency of the two most important vibration modes, and then set reasonable damping coefficients at their angular
frequency values. Note that a poor choice of these damping parameters could lead to over-damped, non-physical behaviour of the beam, so it is essential to choose them
carefully.</p>
<p>The flexible beam model is built on the assumption that the elastic deformations are small, allowing to describe them by linear equations.
However, it is possible to represent large deformations correctly by connecting several beam models in series, as demonstrated in the
<a href=\"modelica://Flex.Examples.LargeDeformation\">LargeDeformation</a> example case.</p>
<h3>References</h3>
<p>
[1] F. Schiavo, G. Ferretti, L. Vigan&ograve;, Object-Oriented Modelling and Simulation of 
Flexible Multibody Thin Beams in Modelica with the Finite Element Method, 
Proc. 4th International Modelica Conference, Hamburg, 2005, pp. 25-34,
<a href=\"https://modelica.org/events/Conference2005/online_proceedings/Session1/Session1a2.pdf\">online</a>.<br>
[2] F. Schiavo, L. Vigan&ograve;, G. Ferretti, Object-oriented modelling of flexible beams, 
Multibody Systems Dynamics, vol.5, n. 3, 2006, pp. 263-286, <a href=\"https://doi.org/10.1007/s11044-006-9012-8\">
DOI:10.1007/s11044-006-9012-8</a>.
</p>
</html>"));
end FlexibleThinBeam;
