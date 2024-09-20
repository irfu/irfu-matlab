%
% Hard-coded constants generated by code.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef sconst
% PROPOSAL: Merge class back into bicas.const.



  %###########
  %###########
  % CONSTANTS
  %###########
  %###########
  properties(Constant)



    C = bicas.sconst.init_const();



  end    % properties(Constant)



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static, Access=public)



    % Defines all constants for ASID, SSID, SDID, and Routing.
    %
    % IMPLEMENTATION NOTE: Defining one constant struct, which contains
    % multiple constants as fields. Makes it possible to access constants
    % through a variable copy of this constant rather than using the long
    % qualifiers.
    function C = init_const()
      % PROPOSAL: Naming convention that is consistent for k<-->s dictionaries.
      %   PROPOSAL: "O" = object = ASID/SSID/SDID/Routing

      C = struct();
      C.S_ASID_DICT    = configureDictionary('string', 'bicas.proc.L1L2.AntennaSignalId');
      C.S_SSID_DICT    = configureDictionary('string', 'bicas.proc.L1L2.SignalSourceId');
      C.S_SDID_DICT    = configureDictionary('string', 'bicas.proc.L1L2.SignalDestinationId');
      C.S_ROUTING_DICT = configureDictionary('string', 'bicas.proc.L1L2.Routing');

      C.SSID_S_K_DICT  = configureDictionary('string', 'uint8');
      C.SDID_S_K_DICT  = configureDictionary('string', 'uint8');
      C.SSID_K_DICT    = configureDictionary('bicas.proc.L1L2.SignalSourceId',      'uint8');
      C.SDID_K_DICT    = configureDictionary('bicas.proc.L1L2.SignalDestinationId', 'uint8');

      % Global list of k values for ALL classes. For avoiding collisions.
      kGlobalArray = [];

      function add_ASR(fn, k, asidCategory, asidAntennas)
        Asid = bicas.proc.L1L2.AntennaSignalId(asidCategory, asidAntennas);
        Ssid = bicas.proc.L1L2.SignalSourceId(Asid);

        C.S_ASID_DICT(fn) = Asid;
        add_SSID(fn, k, Ssid)
        add_SDID(fn, k, bicas.proc.L1L2.SignalDestinationId(Asid))
        C.S_ROUTING_DICT(fn) = bicas.proc.L1L2.Routing(Ssid);
      end

      function assert_new_k(k)
        assert(isnumeric(k))
        assert(~ismember(k, kGlobalArray))
        kGlobalArray(end+1) = k;
      end

      function add_SSID(fn, k, Ssid)
        k2 = k+100;
        assert_new_k(k2)

        C.S_SSID_DICT(fn)   = Ssid;
        C.SSID_S_K_DICT(fn) = k2;
        %C.K_SSID_DICT = add_key_value(C.K_SSID_DICT, k+100, Ssid);
        C.SSID_K_DICT(Ssid) = uint8(k2);
      end

      function add_SDID(fn, k, Sdid)
        k2 = k+200;
        assert_new_k(k2)

        C.S_SDID_DICT(fn)   = Sdid;
        C.SDID_S_K_DICT(fn) = k2;
        %C.K_SDID_DICT = add_key_value(C.K_SDID_DICT, k+200, Sdid);
        C.SDID_K_DICT(Sdid) = uint8(k2);
      end

      function main()
        % =====================================
        % Add every possible unique ASID object
        % =====================================
        add_ASR("DC_V1",  1, 'DC_SINGLE', [1   ]);
        add_ASR("DC_V2",  2, 'DC_SINGLE', [2   ]);
        add_ASR("DC_V3",  3, 'DC_SINGLE', [3   ]);

        add_ASR("DC_V12", 4, 'DC_DIFF',   [1, 2]);
        add_ASR("DC_V13", 5, 'DC_DIFF',   [1, 3]);
        add_ASR("DC_V23", 6, 'DC_DIFF',   [2, 3]);

        add_ASR("AC_V12", 7, 'AC_DIFF',   [1, 2]);
        add_ASR("AC_V13", 8, 'AC_DIFF',   [1, 3]);
        add_ASR("AC_V23", 9, 'AC_DIFF',   [2, 3]);

        add_SSID("REF25V",  10, bicas.proc.L1L2.SignalSourceId('2.5V_REF'));
        add_SSID("GND",     11, bicas.proc.L1L2.SignalSourceId('GND'));
        add_SSID("UNKNOWN", 12, bicas.proc.L1L2.SignalSourceId('UNKNOWN'));

        add_SDID("NOWHERE", 13, bicas.proc.L1L2.SignalDestinationId('NOWHERE'));

        C.S_ROUTING_DICT("REF25V_TO_DC_V1")    = bicas.proc.L1L2.Routing(C.S_SSID_DICT("REF25V"),  C.S_SDID_DICT("DC_V1"));
        C.S_ROUTING_DICT("REF25V_TO_DC_V2")    = bicas.proc.L1L2.Routing(C.S_SSID_DICT("REF25V"),  C.S_SDID_DICT("DC_V2"));
        C.S_ROUTING_DICT("REF25V_TO_DC_V3")    = bicas.proc.L1L2.Routing(C.S_SSID_DICT("REF25V"),  C.S_SDID_DICT("DC_V3"));
        C.S_ROUTING_DICT("GND_TO_DC_V1")       = bicas.proc.L1L2.Routing(C.S_SSID_DICT("GND"),     C.S_SDID_DICT("DC_V1"));
        C.S_ROUTING_DICT("GND_TO_DC_V2")       = bicas.proc.L1L2.Routing(C.S_SSID_DICT("GND"),     C.S_SDID_DICT("DC_V2"));
        C.S_ROUTING_DICT("GND_TO_DC_V3")       = bicas.proc.L1L2.Routing(C.S_SSID_DICT("GND"),     C.S_SDID_DICT("DC_V3"));
        C.S_ROUTING_DICT("UNKNOWN_TO_NOWHERE") = bicas.proc.L1L2.Routing(C.S_SSID_DICT("UNKNOWN"), C.S_SDID_DICT("NOWHERE"));

        assert(~ismember(255, kGlobalArray))
        assert(~ismember(0,   kGlobalArray))
      end

      main()
    end



  end    % methods(Static, Access=public)



end
