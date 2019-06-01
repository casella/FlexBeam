within Flex.Examples;
model SliderCrank2
  extends SliderCrank(torque(y=if time < 0.7 then 0.01*(1-exp(-time/0.167)) else
                                                                                0));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<p>Same as <a href=\"modelica://Flex.Examples.SliderCrank\">SliderCrank</a>, except that the applied torque is 
set to zero for time > 0.7.</p>
<p><code>sliderPosition.y</code> and <code>rodMidpointDeformation.y</code> reproduce the results shown in
Fig. 8 and 10 of [1], that were obtained with a different finite element multibody code, thus confirming
the correctness of the flexible beam model.</p>
<h3>References</h3>
<p>
[1] J.L. Escalona, H. A. Hussien and A. A. Shabana, Application of the absolute nodal coordinate formulation
to multibody systems dynamics, Journal of Sound and Vibration, 214(5), 1998, pp. 833-851,
<a href=\"https://doi.org/10.1006/jsvi.1998.1563\">DOI:10.1006/jsvi.1998.1563</a>.</p>
</html>"),
    experiment(StopTime=1.6, __Dymola_NumberOfIntervals=2000));
end SliderCrank2;
