within Flex.Test;

model FreeVibration
  inner Modelica.Mechanics.MultiBody.World world(gravityType = Modelica.Mechanics.MultiBody.Types.GravityTypes.NoGravity)  annotation(
    Placement(transformation(extent = {{-10, -10}, {10, 10}})));

 FlexibleThinBeam beam(
    L=2, 
    A=0.01^2,
    E=7.2e10,
    J=0.01^4/12,
    N=3,
    rho=2700,
    clampedFree=true,
    circularSection=false, qf(each fixed = true), dqf(each fixed = true)) annotation(
    Placement(transformation(origin = {36, 0}, extent = {{-10, -10}, {10, 10}})));equation
  connect(world.frame_b, beam.frame_a) annotation(
    Line(points = {{10, 0}, {26, 0}}));

end FreeVibration;
