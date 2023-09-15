classdef testC4 < matlab.unittest.TestCase
  % TESTC4 test c_4_grad routine

  properties
  end

  methods (Static)
    function R=defineSatellitePosition(scSeparation)
      %
      R1 = scSeparation * [1 0 0];
      R2 = scSeparation * [cos(1/3*2*pi) sin(1/3*2*pi) 0];
      R3 = scSeparation * [cos(2/3*2*pi) sin(2/3*2*pi) 0];
      R4 = scSeparation * [0 0 sqrt(norm(R2 - R1)^2-1)];
      zMassCenter=[0 0 R4(3)/4];
      R.C1 = R1 - zMassCenter;
      R.C2 = R2 - zMassCenter;
      R.C3 = R3 - zMassCenter;
      R.C4 = R4 - zMassCenter;
    end
  end
  % 	% construction of field
  % 	Bconst	=@(x) [0 0 1]; %
  methods (Test)
    function testCurvature(testCase)
      %% SUBTEST: position of tetrahedron, base plane in xy
      % Cylindrical field with symmetry axis along z and center at rc
      % such field gives curvature equal to 1/R where R is distance to
      % the cylinder symmetry axis
      Bcircle	=@(x,rc) [(x(:,2)-rc(2))./sqrt((x(:,1)-rc(1)).^2 +(x(:,2)-rc(2)).^2) ...
        -(x(:,1)-rc(1))./sqrt((x(:,1)-rc(1)).^2 +(x(:,2)-rc(2)).^2) 0];
      scSeparation = 1;
      R = testC4.defineSatellitePosition(scSeparation);
      symmetryAxisLocation = [10 10 0]; % cylinder symmetry axis location, B field is clockwise
      distanceSymmetryAxisToSpacecraft = norm(symmetryAxisLocation);
      B.C1=Bcircle(R.C1,symmetryAxisLocation);
      B.C2=Bcircle(R.C2,symmetryAxisLocation);
      B.C3=Bcircle(R.C3,symmetryAxisLocation);
      B.C4=Bcircle(R.C4,symmetryAxisLocation);
      curvature = c_4_grad(R,B,'curvature');
      curvatureTheoretical = 1/distanceSymmetryAxisToSpacecraft^2*symmetryAxisLocation;
      testCase.verifyEqual(curvature,curvatureTheoretical,...
        'AbsTol',sqrt(scSeparation / distanceSymmetryAxisToSpacecraft));
    end
    function testConstantCurrentInZ(testCase)
      jz  = 10;
      Bjz = @(x) [0 jz*x(1) 0];
      scSeparation = 1;
      R   = testC4.defineSatellitePosition(scSeparation);
      B.C1 = Bjz(R.C1);
      B.C2 = Bjz(R.C2);
      B.C3 = Bjz(R.C3);
      B.C4 = Bjz(R.C4);
      actCurl = c_4_grad(R,B,'curl');
      expCurl = [0 0 jz];
      testCase.verifyEqual(actCurl,expCurl,'AbsTol',2*eps);
    end
  end
  %
  % 	% Test drift gradient
  % 	% Fitzpatrick book page 30
  % 	% V=mv^2/2/e B x gradB / B^3
  %

end

