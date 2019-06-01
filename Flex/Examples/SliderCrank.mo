within Flex.Examples;
model SliderCrank
  Modelica.Blocks.Sources.RealExpression rodMidpointDeformation(y=-rod.r0shape[
        4, 2]) annotation (Placement(transformation(extent={{-14,54},{14,66}})));
  Modelica.Blocks.Sources.RealExpression sliderPosition(y=slider.frame_a.r_0[1])
    annotation (Placement(transformation(extent={{-16,68},{12,82}})));
  inner Modelica.Mechanics.MultiBody.World world(g=0, enableAnimation=true,
    animateWorld=false,
    animateGravity=false)
    annotation (Placement(transformation(extent={{-84,-32},{-64,-12}}, rotation=
           0)));
  Modelica.Mechanics.MultiBody.Joints.Revolute actuatedRevolute(
    useAxisFlange=true,
    stateSelect=StateSelect.avoid,
    phi(fixed=true),
    w(fixed=true))
                 annotation (Placement(transformation(
        origin={-46,4},
        extent={{-10,-10},{10,10}},
        rotation=90)));
  Modelica.Mechanics.MultiBody.Joints.Prismatic prismatic2(
    stateSelect=StateSelect.prefer,
    s(start=0),
    animation=false,
    n={-1,0,0}) annotation (Placement(transformation(extent={{-10,-32},{10,-12}},
          rotation=0)));
  Modelica.Mechanics.Rotational.Sources.Torque torqueGenerator(useSupport=false)
    annotation (Placement(transformation(extent={{-84,-6},{-64,14}}, rotation=0)));
  Modelica.Mechanics.MultiBody.Joints.Revolute revolute1(stateSelect=
        StateSelect.avoid)
    annotation (Placement(transformation(extent={{-6,20},{14,40}}, rotation=0)));
  Flex.FlexibleThinBeam rod(
    rho=2770,
    A=7.854*1e-5,
    J=4.909*1e-10,
    ClampedFree=false,
    E=5e7,
    N=8,
    L=0.304) annotation (Placement(transformation(extent={{34,20},{54,40}},
          rotation=0)));
  Modelica.Mechanics.MultiBody.Joints.RevolutePlanarLoopConstraint
                                               revolute2
                                      annotation (Placement(transformation(
        origin={46,-22},
        extent={{10,-10},{-10,10}},
        rotation=180)));
  Flex.FlexibleThinBeam crank(
    rho=2770,
    L=0.152,
    A=7.854*1e-5,
    E=1e9,
    J=4.909*1e-10,
    N=3) annotation (Placement(transformation(extent={{-38,20},{-18,40}},
          rotation=0)));
  Modelica.Mechanics.MultiBody.Parts.FixedRotation fixedRotation1(
    n={0,0,1},
    animation=false,
    r={0,0,0},
    angle=180) annotation (Placement(transformation(
        origin={78,6},
        extent={{-10,-10},{10,10}},
        rotation=90)));
  Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation1(r={0.152 + 0.304,0,0},
      animation=false) annotation (Placement(transformation(extent={{-38,-32},{
            -18,-12}}, rotation=0)));
  Modelica.Mechanics.MultiBody.Parts.BodyCylinder slider(
    density=0,
    r={0.15,0,0},
    diameter=0.15) annotation (Placement(transformation(extent={{30,-58},{50,
            -38}}, rotation=0)));
  Modelica.Mechanics.MultiBody.Visualizers.FixedShape cylinder(
    shapeType="pipecylinder",
    width=0.2,
    height=0.2,
    extra=0.15/0.2,
    r_shape={0.152 + 0.15,0,0},
    length=0.304 + 0.15 - 0.15,
    color={255,255,255},
    animation=false) annotation (Placement(transformation(extent={{-28,-66},{-8,
            -46}}, rotation=0)));
  Modelica.Blocks.Sources.RealExpression torque(y=if time < 0.7 then 0.01*(1 -
        exp(-time/0.167)) else 0)
    annotation (Placement(transformation(extent={{-122,-6},{-102,14}})));
equation
  connect(actuatedRevolute.frame_a, world.frame_b) annotation (Line(
      points={{-46,-6},{-46,-22},{-64,-22}},
      color={0,0,0},
      thickness=0.5));
  connect(crank.FrameB, revolute1.frame_a) annotation (Line(
      points={{-18.3,30.1},{-14,30.1},{-14,30},{-6,30}},
      color={0,0,0},
      thickness=0.5));
  connect(crank.FrameA, actuatedRevolute.frame_b) annotation (Line(
      points={{-37.7,30.1},{-46,30.1},{-46,14}},
      color={0,0,0},
      thickness=0.5));
  connect(revolute1.frame_b, rod.FrameA) annotation (Line(
      points={{14,30},{24.3,30},{24.3,30.1},{34.3,30.1}},
      color={0,0,0},
      thickness=0.5));
  connect(fixedTranslation1.frame_a, world.frame_b) annotation (Line(
      points={{-38,-22},{-64,-22}},
      color={0,0,0},
      thickness=0.5));
  connect(fixedTranslation1.frame_b,prismatic2. frame_a) annotation (Line(
      points={{-18,-22},{-10,-22}},
      color={0,0,0},
      thickness=0.5));
  connect(revolute2.frame_a,prismatic2. frame_b) annotation (Line(
      points={{36,-22},{10,-22}},
      color={0,0,0},
      thickness=0.5));
  connect(revolute2.frame_b,fixedRotation1. frame_a) annotation (Line(
      points={{56,-22},{78,-22},{78,-4}},
      color={0,0,0},
      thickness=0.5));
  connect(fixedRotation1.frame_b, rod.FrameB) annotation (Line(
      points={{78,16},{78,30.1},{53.7,30.1}},
      color={0,0,0},
      thickness=0.5));
  connect(slider.frame_a,prismatic2. frame_b) annotation (Line(
      points={{30,-48},{22,-48},{22,-22},{10,-22}},
      color={0,0,0},
      thickness=0.5));
  connect(cylinder.frame_a, world.frame_b) annotation (Line(
      points={{-28,-56},{-46,-56},{-46,-22},{-64,-22}},
      color={0,0,0},
      thickness=0.5));
  connect(torqueGenerator.flange, actuatedRevolute.axis)
    annotation (Line(points={{-64,4},{-56,4}}, color={0,0,0}));
  connect(torque.y, torqueGenerator.tau)
    annotation (Line(points={{-101,4},{-86,4}}, color={0,0,127}));
  annotation (Diagram(coordinateSystem(extent={{-140,-100},{100,100}})),
    Icon(coordinateSystem(extent={{-100,-100},{100,100}})),
    experiment(StopTime=1.6, __Dymola_NumberOfIntervals=1000),
    Documentation(info="<html>
<p>This model replicates the test case presented in [1]. The connecting rod has the same cross-section
and density of the crank, but it is assumed to be 20 times less stiff, while the slider block is
assumed to be rigid and massless. The crank is modelled by 3 finite elements, the rod by 8 finite
elements.</p>

<p>In the initial position, both the connecting rod and crank are assumed to be
horizontal. The mechanism is assumed to be driven by a moment applied at the crankshaft, that increases
with a negative exponential time law.</p>

<p><code>sliderPosition.y</code> and <code>rodMidpointDeformation.y</code> reproduce the results shown in
Fig. 7 and 9 of [1], that were obtained with a different finite element multibody code, thus confirming
the correctness of the flexible beam model.</p>

<p>References</p>
[1] J.L. Escalona, H. A. Hussien and A. A. Shabana, Application of the absolute nodal coordinate formulation
to multibody systems dynamics, Journal of Sound and Vibration, 214(5), 1998, pp. 833-851,
<a href=\"https://doi.org/10.1006/jsvi.1998.1563\">DOI:10.1006/jsvi.1998.1563</a>.<br>
[2] F. Schiavo, G. Ferretti, L. Vigan&ograve;, Object-Oriented Modelling and Simulation of 
Flexible Multibody Thin Beams in Modelica with the Finite Element Method, 
Proc. 4th International Modelica Conference, Hamburg, 2005, pp. 25-34,
<a href=\"https://modelica.org/events/Conference2005/online_proceedings/Session1/Session1a2.pdf\">online</a>.<br>
[3] F. Schiavo, L. Vigan&ograve;, G. Ferretti, Object-oriented modelling of flexible beams, 
Multibody Systems Dynamics, vol.5, n. 3, 2006, pp. 263-286, <a href=\"https://doi.org/10.1007/s11044-006-9012-8\">
DOI:10.1007/s11044-006-9012-8</a>.
</html>"));
end SliderCrank;
