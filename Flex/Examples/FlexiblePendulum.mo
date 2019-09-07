within Flex.Examples;
model FlexiblePendulum
  inner Modelica.Mechanics.MultiBody.World world(animateGravity=false, animateWorld=false)
                              annotation (Placement(transformation(extent={{-64,
            -4},{-44,16}}, rotation=0)));
  Modelica.Mechanics.MultiBody.Joints.Revolute revolute1(phi(fixed=true), w(
        fixed=true))                                     annotation (Placement(transformation(
          extent={{-28,-4},{-8,16}}, rotation=0)));
  FlexibleThinBeam flexiblePendulum(
    rho=5540,
    circularSection=false,
    A=0.0018,
    E=1e9,
    J=1.215e-8,
    csi1=0.01,
    csi2=50,
    omega1=0.01,
    omega2=100,
    L=0.4,
    N=10,
    dqf(fixed=true),
    ddqf(fixed=true))
    annotation (Placement(transformation(extent={{0,-4},{40,16}}, rotation=0)));
  Modelica.Mechanics.MultiBody.Joints.Revolute revolute2(phi(fixed=true), w(
        fixed=true))                                     annotation (Placement(transformation(
          extent={{-28,-30},{-8,-10}},
                                     rotation=0)));
  Modelica.Mechanics.MultiBody.Parts.BodyBox rigidPendulum(
    r={0.4,0,0},
    width=0.02,
    height=0.005)
    annotation (Placement(transformation(extent={{10,-30},{30,-10}})));
equation
  connect(revolute1.frame_a, world.frame_b) annotation (Line(
      points={{-28,6},{-44,6}},
      color={0,0,0},
      thickness=0.5));
  connect(revolute1.frame_b, flexiblePendulum.frame_a) annotation (Line(
      points={{-8,6},{-1.7,6},{-1.7,6.1},{0.6,6.1}},
      color={0,0,0},
      thickness=0.5));
  connect(world.frame_b, revolute2.frame_a) annotation (Line(
      points={{-44,6},{-36,6},{-36,-20},{-28,-20}},
      color={95,95,95},
      thickness=0.5));
  connect(revolute2.frame_b, rigidPendulum.frame_a) annotation (Line(
      points={{-8,-20},{10,-20}},
      color={95,95,95},
      thickness=0.5));
  annotation (
    experiment(StopTime=5), uses(Modelica(version="3.2.2")));
end FlexiblePendulum;
