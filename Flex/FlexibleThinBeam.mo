within Flex;
model FlexibleThinBeam "Flexible thin beam model"
  import SI = Modelica.SIunits;
  import Cv = Modelica.SIunits.Conversions;
  import Modelica.Mechanics.MultiBody.Frames;

  parameter SI.Density rho=7800 "Material density";
  parameter SI.Length L "Beam length";

  parameter SI.Position r_0_start[3]={0,0,0}
    "Initial values of frame_a.r_0 (vector from origin of world frame to origin of frame_a resolved in world frame)"
    annotation (Dialog(tab="Initialization"));
  parameter Real qf_start[3*N]=zeros(3*N)
    "Initial beam deformation"
    annotation(Dialog(tab="Initialization"));
  parameter Real dqf_start[3*N]=zeros(3*N)
    "Initial velocity of beam deformation"
    annotation(Dialog(tab="Initialization"));
  parameter Real ddqf_start[3*N]=zeros(3*N)
    "Initial acceleration of beam deformation"
    annotation(Dialog(tab="Initialization"));
  parameter Boolean ClampedFree=true
    "Clamped-Free model if true, else simply-supported model"
    annotation(Dialog(tab="Boundary Conditions"));
  parameter Boolean CircularSection=true
    "CircularSection if true, else Rectangular"
    annotation(Dialog(tab="3D Graphics"));
  parameter Modelica.Mechanics.MultiBody.Types.Color ColorBeam={128,128,128}
    "Beam color "
    annotation(Dialog(tab="3D Graphics"));
  parameter SI.Area A "Cross sectional area";
  parameter SI.ModulusOfElasticity E "Young modulus of the beam material";
  parameter SI.SecondMomentOfArea J "Second moment of area of cross section around z-axis";
  parameter Real alpha(unit="1/s")=0
    "Rayleigh structural damping proportional to mass [sec^-1]";
  parameter Real beta(unit="s")=0
    "Rayleigh structural damping proportional to stiffness [sec]";
  parameter Integer N(min=1) = 5 "Number of elements";
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

public
  Real qf[3*N](start=qf_start) "Elastic coordinates";
  Real dqf[3*N](start=dqf_start) "Elastic velocities";
  Real ddqf[3*N](start=ddqf_start) "Elastic accelerations";
  SI.Acceleration g_0[3] "Gravity acceleration resolved in world frame";

protected
  Real B[N, 6, 3*N];
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

  /* 3D graphics variables */
  type vec3D = Real[3];
  vec3D r0shape[N];
  vec3D rrelshape[N];
  Real Lshape[N];

public
  Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a annotation (Placement(
        transformation(extent={{-106,-8},{-86,12}}, rotation=0),
        iconTransformation(extent={{-108,-10},{-86,12}})));
  Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b annotation (Placement(
        transformation(extent={{88,-8},{108,12}}, rotation=0),
        iconTransformation(extent={{86,-10},{108,12}})));

protected
  Modelica.Mechanics.MultiBody.Frames.Orientation R_rel;
  outer Modelica.Mechanics.MultiBody.World world;
  Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape Segment[N](
    r_shape=r0shape,
    length=Lshape,
    each width=if CircularSection then 2*sqrt(A/pi) else A/(0.4*sqrt(A)),
    each height=if CircularSection then 2*sqrt(A/pi) else (0.4*sqrt(A)),
    lengthDirection=rrelshape,
    each widthDirection={0,0,1},
    each shapeType=if CircularSection then "cylinder" else "box",
    each color=ColorBeam,
    each extra=0.0,
    each r=frame_a.r_0,
    each R=frame_a.R);

equation
  Connections.branch(frame_a.R, frame_b.R);
  assert(cardinality(frame_a) > 0 or cardinality(frame_b) > 0,
    "Neither connector frame_a nor frame_b of FlexBeamFem object is connected");

  //connectivity matrices
  if ClampedFree then
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
