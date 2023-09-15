classdef TestTimeArray < matlab.unittest.TestCase
  % Test cases for GenericTimeArray and derived classes

  % ----------------------------------------------------------------------------
  % SPDX-License-Identifier: Beerware
  % "THE BEER-WARE LICENSE" (Revision 42):
  % <yuri@irfu.se> wrote this file.  As long as you retain this notice you
  % can do whatever you want with this stuff. If we meet some day, and you think
  % this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
  % ----------------------------------------------------------------------------

  properties (TestParameter)
    class = {'EpochUnix','EpochTT','EpochCdf','EpochCdf16'};
    utc = {'1970-01-01T00:00:00.000000Z',...
      '1998-07-17T12:00:00.123456789Z',...
      '2002-03-04T11:59:59.999999999Z',...
      '2015-08-01T12:00:00'};
    ttns = {int64(-946727959814622001), ...
      int64(-46051136692543211), ...
      int64(68515264183999999), int64(491702468184000000)};
    unixEpoch = {0, 900676800.1234568, 1015243200, 1438430400};
  end

  methods (Test, ParameterCombination='sequential')
    function testEpochUnixConstructorFromUtc(testCase,utc,unixEpoch)
      t = EpochUnix(utc);
      testCase.verifyEqual(t.epoch, unixEpoch)
    end
    function testEpochUnixConstructorFromTtns(testCase,ttns,unixEpoch)
      t = EpochUnix(ttns);
      testCase.verifyEqual(t.epoch, unixEpoch)
    end
    function testEpochTTConstructorFromUtc(testCase,utc,ttns)
      t = EpochTT(utc);
      testCase.verifyEqual(t.epoch, ttns)
    end
  end
end