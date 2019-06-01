within Flex.Examples;
model CartFlexibleBeam
 parameter Modelica.SIunits.Length L=2;
  parameter Real lato=0.01;
  parameter Modelica.SIunits.Area A=lato^2;
  parameter Real E=2.1e11;
  parameter Real J=lato^4/12;
  parameter Integer N=4 "number of elements";
  inner Modelica.Mechanics.MultiBody.World world(
    enableAnimation=true,
    animateWorld=false,
    animateGravity=false)     annotation (Placement(transformation(extent={{-64,-84},
            {-48,-68}},      rotation=0)));
  FlexibleThinBeam Beam2(
    A=A,
    E=E,
    J=J,
    N=N,
    qf_start=zeros(3*N),
    ColorBeam={128,128,128},
    L=L/3,
    Alpha=0.05,
    Beta=5e-4) annotation (Placement(transformation(
        origin={16,54},
        extent={{-10,-10},{10,10}},
        rotation=90)));
  Modelica.Mechanics.MultiBody.Parts.Body Body1(r_CM={0,0,0}, m=1)
                                                annotation (Placement(transformation(
        origin={16,108},
        extent={{-10,-10},{10,10}},
        rotation=90)));
  FlexibleThinBeam Beam3(
    A=A,
    E=E,
    J=J,
    N=N,
    qf_start=zeros(3*N),
    L=L/3,
    Alpha=0.05,
    Beta=5e-4) annotation (Placement(transformation(
        origin={16,82},
        extent={{-10,-10},{10,10}},
        rotation=90)));
  Modelica.Mechanics.MultiBody.Joints.Prismatic ActuatedPrismatic(
    useAxisFlange=true,
    s_offset=0,
    animation=false) annotation (Placement(transformation(extent={{-32,-86},{-12,
            -66}}, rotation=0)));
  Modelica.Mechanics.MultiBody.Parts.BodyBox cart(
    color={255,0,0},
    r={0.8,0,0},
    width=0.1,
    height=0.4) annotation (Placement(transformation(extent={{6,-86},{26,-66}},
          rotation=0)));
  Modelica.Mechanics.MultiBody.Parts.FixedTranslation FixedTranslation1(animation=false, r={-0.4,
        0.05,0})
    annotation (Placement(transformation(extent={{36,-38},{16,-18}}, rotation=0)));
  Modelica.Mechanics.MultiBody.Parts.FixedRotation FixedRotation1(
    animation=false,
    n={0,0,1},
    angle=90) annotation (Placement(transformation(
        origin={16,-4},
        extent={{-10,-10},{10,10}},
        rotation=90)));
  Modelica.Mechanics.Translational.Sources.Position actuator annotation (
      Placement(transformation(extent={{-42,-62},{-22,-42}}, rotation=0)));
  Modelica.Blocks.Sources.Sine sineSignal(offset=1, freqHz=0.5) annotation (
      Placement(transformation(extent={{-88,-60},{-72,-44}}, rotation=0)));
  FlexibleThinBeam Beam1(
    A=A,
    E=E,
    J=J,
    N=N,
    qf_start=zeros(3*N),
    ColorBeam={128,128,128},
    L=L/3,
    Alpha=0.05,
    Beta=5e-4) annotation (Placement(transformation(
        origin={16,26},
        extent={{-10,-10},{10,10}},
        rotation=90)));
  Modelica.Mechanics.MultiBody.Visualizers.FixedShape wheel1(
    shapeType="cylinder",
    lengthDirection={0,0,-1},
    color={255,255,0},
    animation=true,
    r_shape={-0.1/2,-0.05 - 0.1/2,0.2},
    length=0.05,
    width=0.1,
    height=0.1) annotation (Placement(transformation(extent={{56,-124},{76,-104}},
          rotation=0)));
  Modelica.Mechanics.MultiBody.Visualizers.FixedShape wheel2(
    shapeType="cylinder",
    color={255,255,0},
    animation=true,
    length=0.05,
    width=0.1,
    height=0.1,
    r_shape={-0.1/2,-0.05 - 0.1/2,-0.2},
    lengthDirection={0,0,1}) annotation (Placement(transformation(extent={{56,-68},
            {76,-48}}, rotation=0)));
  Modelica.Mechanics.MultiBody.Visualizers.FixedShape wheel3(
    shapeType="cylinder",
    color={255,255,0},
    animation=true,
    length=0.05,
    width=0.1,
    height=0.1,
    r_shape={0.1/2,-0.05 - 0.1/2,-0.2},
    lengthDirection={0,0,1}) annotation (Placement(transformation(extent={{-18,-28},
            {-38,-8}},  rotation=0)));
  Modelica.Mechanics.MultiBody.Visualizers.FixedShape wheel4(
    shapeType="cylinder",
    color={255,255,0},
    animation=true,
    length=0.05,
    width=0.1,
    height=0.1,
    lengthDirection={0,0,-1},
    r_shape={0.1/2,-0.05 - 0.1/2,0.2}) annotation (Placement(transformation(
          extent={{-16,-124},{-36,-104}},rotation=0)));
equation
  connect(world.frame_b, ActuatedPrismatic.frame_a) annotation (Line(
      points={{-48,-76},{-32,-76}},
      color={0,0,0},
      thickness=0.5));
  connect(ActuatedPrismatic.frame_b, cart.frame_a) annotation (Line(
      points={{-12,-76},{6,-76}},
      color={0,0,0},
      thickness=0.5));
  connect(FixedTranslation1.frame_a, cart.frame_b) annotation (Line(
      points={{36,-28},{42,-28},{42,-52},{26,-52},{26,-76}},
      color={0,0,0},
      thickness=0.5));
  connect(FixedTranslation1.frame_b, FixedRotation1.frame_a) annotation (Line(
      points={{16,-28},{16,-14}},
      color={0,0,0},
      thickness=0.5));
  connect(actuator.flange, ActuatedPrismatic.axis)
    annotation (Line(points={{-22,-52},{-14,-52},{-14,-70}}, color={0,191,0}));
  connect(Beam3.frame_b, Body1.frame_a) annotation (Line(
      points={{15.9,91.7},{15.9,94.2},{16,94.2},{16,98}},
      color={0,0,0},
      thickness=0.5));
  connect(Beam2.frame_b, Beam3.frame_a) annotation (Line(
      points={{15.9,63.7},{15.9,72.3}},
      color={0,0,0},
      thickness=0.5));
  connect(Beam1.frame_b, Beam2.frame_a) annotation (Line(
      points={{15.9,35.7},{15.9,44.3}},
      color={0,0,0},
      thickness=0.5));
  connect(FixedRotation1.frame_b, Beam1.frame_a) annotation (Line(
      points={{16,6},{16,16.3},{15.9,16.3}},
      color={0,0,0},
      thickness=0.5));
  connect(wheel1.frame_a, cart.frame_b) annotation (Line(
      points={{56,-114},{26,-114},{26,-76}},
      color={0,0,0},
      thickness=0.5));
  connect(wheel2.frame_a, cart.frame_b) annotation (Line(
      points={{56,-58},{42,-58},{42,-76},{26,-76}},
      color={0,0,0},
      thickness=0.5));
  connect(wheel3.frame_a, cart.frame_a) annotation (Line(
      points={{-18,-18},{6,-18},{6,-76}},
      color={0,0,0},
      thickness=0.5));
  connect(wheel4.frame_a, cart.frame_a) annotation (Line(
      points={{-16,-114},{6,-114},{6,-76}},
      color={0,0,0},
      thickness=0.5));
  connect(sineSignal.y, actuator.s_ref)
    annotation (Line(points={{-71.2,-52},{-44,-52}}, color={0,0,127}));
  annotation (uses(Modelica(version="3.2.2")),
      Diagram(coordinateSystem(extent={{-100,-140},{100,120}})),
    experiment(StopTime=10),
    Icon(coordinateSystem(extent={{-100,-100},{100,100}})));
end CartFlexibleBeam;
