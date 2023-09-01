model test
  extends OpenModelica;
  Modelica.Blocks.Sources.Ramp ramp(duration = 15, height = 1, offset = 0, startTime = 5)  annotation(
    Placement(visible = true, transformation(origin = {-10, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput momentumModelica annotation(
    Placement(visible = true, transformation(origin = {50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput power annotation(
    Placement(visible = true, transformation(origin = {-28, -42}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-28, -42}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Continuous.Integrator integrator annotation(
    Placement(visible = true, transformation(origin = {38, -46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput powerIntegral annotation(
    Placement(visible = true, transformation(origin = {82, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {82, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Ramp ramp1(duration = 15, height = 100, offset = 850, startTime = 5) annotation(
    Placement(visible = true, transformation(origin = {-42, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput temperatureModelica annotation(
    Placement(visible = true, transformation(origin = {38, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput HModelica annotation(
    Placement(visible = true, transformation(origin = {40, 74}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Ramp ramp2(duration = 5, height = -9990, offset = 10000, startTime = 15) annotation(
    Placement(visible = true, transformation(origin = {-40, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput hxPower annotation(
    Placement(visible = true, transformation(origin = {-68, -74}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-28, -42}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Continuous.Integrator integrator1 annotation(
    Placement(visible = true, transformation(origin = {-2, -78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput hxPowerIntegral annotation(
    Placement(visible = true, transformation(origin = {42, -82}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {82, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(ramp.y, momentumModelica) annotation(
    Line(points = {{2, -2}, {50, -2}, {50, 0}}, color = {0, 0, 127}));
  connect(power, integrator.u) annotation(
    Line(points = {{-28, -42}, {3, -42}, {3, -46}, {26, -46}}, color = {0, 0, 127}));
  connect(integrator.y, powerIntegral) annotation(
    Line(points = {{50, -46}, {82, -46}, {82, -50}}, color = {0, 0, 127}));
  connect(ramp1.y, temperatureModelica) annotation(
    Line(points = {{-31, 30}, {6.5, 30}, {6.5, 32}, {38, 32}}, color = {0, 0, 127}));
  connect(ramp2.y, HModelica) annotation(
    Line(points = {{-29, 72}, {8.5, 72}, {8.5, 74}, {40, 74}}, color = {0, 0, 127}));
  connect(hxPower, integrator1.u) annotation(
    Line(points = {{-68, -74}, {-37, -74}, {-37, -78}, {-14, -78}}, color = {0, 0, 127}));
  connect(integrator1.y, hxPowerIntegral) annotation(
    Line(points = {{9, -78}, {41, -78}, {41, -82}}, color = {0, 0, 127}));
  annotation(
    uses(Modelica(version = "3.2.3")));
end test;