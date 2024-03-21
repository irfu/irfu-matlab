classdef test_geocentric_coordinate_transformation < matlab.unittest.TestCase
  %TEST_GEOCENTRIC_COORDINATE_TRANSFORMATION
  properties (ClassSetupParameter)
    % Change rand(2,1) to rand(X,1) where X is the number of vectors to run
    % tests on for each permutation of coordinate systems.
    % Default 2 means each possible permutation is run twice.
    % Use random dates a maximum of 365*50 days, ie 50 years, before today.
    t = struct('time', irf_time(now - 365*50*rand(2, 1), 'datenum>tt'));
  end
  properties (TestParameter)
    coord0 = {'gei', 'geo', 'gse', 'gsm', 'mag'};
    coord1 = {'gei', 'geo', 'gse', 'gsm', 'mag'};
    coord2 = {'gei', 'geo', 'gse', 'gsm', 'mag'};
    coord3 = {'gei', 'geo', 'gse', 'gsm', 'mag'};
  end

  methods (Test, ParameterCombination='exhaustive')

    function test_cyclic_transformations(testCase, coord0, coord1, coord2, coord3)
      % Test transforming from "coord0">"coord1">"coord2">"coord3" and
      % then finally back to "coord0" and verify result is the same.
      if ( strcmp(coord0, coord1) || strcmp(coord1, coord2) || ...
          strcmp(coord2, coord3) || strcmp(coord3, coord0) )
        % No transformations to speak of in these, skip it.
        irf.log('debug', ['Skipping: ', coord0, '>', coord1, '>', ...
          coord2, '>', coord3, '>', coord0]);
      else
        irf.log('debug', ['Each test consists of ', num2str(length(testCase.t.time)),' separate times.']);
        vecStart = [testCase.t.time, rand(length(testCase.t.time),3)*rand(1)*1e4];
        irf.log('debug', ['First row of vecStart = ', num2str(vecStart(1,:))]);
        % First coord. transformation
        vec = irf.geocentric_coordinate_transformation(vecStart, [coord0 '>' coord1]);
        irf.log('debug', ['First row of Vec1 = ' num2str(vec(1,:))]);
        % Second
        vec = irf.geocentric_coordinate_transformation(vec, [coord1 '>' coord2]);
        irf.log('debug', ['First row of Vec2 = ', num2str(vec(1,:))]);
        % Third
        vec = irf.geocentric_coordinate_transformation(vec, [coord2 '>' coord3]);
        irf.log('debug', ['First row of Vec3 = ', num2str(vec(1,:))]);
        % And back again
        vec = irf.geocentric_coordinate_transformation(vec, [coord3 '>' coord0]);
        irf.log('debug', ['First row of Vec = ', num2str(vec(1,:))]);
        % Same result?
        testCase.verifyEqual(vec, vecStart, 'AbsTol', 1e-10);
      end
    end
  end
end
