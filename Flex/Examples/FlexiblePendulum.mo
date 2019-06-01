within Flex.Examples;
model FlexiblePendulum
  inner Modelica.Mechanics.MultiBody.World world(animateGravity=false, animateWorld=false)
                              annotation (Placement(transformation(extent={{-64,
            -4},{-44,16}}, rotation=0)));
  Modelica.Mechanics.MultiBody.Joints.Revolute Revolute2 annotation (Placement(transformation(
          extent={{-28,-4},{-8,16}}, rotation=0)));
  FlexibleThinBeam FEMPendulum(
    rho=5540,
    A=0.0018,
    E=1e9,
    J=1.215e-8,
    CircularSection=false,
    L=0.4,
    N=10,
    Alpha=0.02,
    Beta=6e-5,
    ClampedFree=true)
    annotation (Placement(transformation(extent={{0,-4},{40,16}}, rotation=0)));
equation
  connect(Revolute2.frame_a, world.frame_b) annotation (Line(
      points={{-28,6},{-44,6}},
      color={0,0,0},
      thickness=0.5));
  connect(Revolute2.frame_b, FEMPendulum.FrameA) annotation (Line(
      points={{-8,6},{-1.7,6},{-1.7,6.1},{0.6,6.1}},
      color={0,0,0},
      thickness=0.5));
  annotation (
    experiment(StopTime=5), uses(Modelica(version="3.2.2")));
end FlexiblePendulum;
