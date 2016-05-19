{
    "identification": {
        "project":      "ROC-SGSE",
        "name":         "BIASRCS (temporary name)",
        "identifier":   "ROC-SGSE-BIASRCS",
        "description":  "BIAS calibration software (temporary description)"
    },
    "release":      {
        "version":      "0.0.1",
        "date":         "2016-05-19",
        "author":       "Erik P G Johansson",
        "contact":      "erik.johansson@irfu.se",
        "institute":    "IRF-U",
        "modification": "None (Initial release)"
    },
    "environment":  {
        "executable":   "bin/biasrcs"
    },
    "modes":        [
    {
        "name":         "testmode1",
        "purpose":      "Mode 1 (temporary purpose description)",
        "inputs":       {
            "input_SCI":    {
                "identifier":   "ROC-SGSE_L2R_RPW-LFR-SBM1-CWF",
                "version":      "01"
            },
            "input_HK":     {
                "identifier":   "ROC-SGSE_HK_RPW-BIA",
                "version":      "01"
            }
        },
        "outputs":      {
            "output_SCI":   {
                "identifier":   "ROC-SGSE_L2S_RPW-BIA-xxxxx",
                "name":         "xxxxx (temporary name)",
                "description":  "Contains xxxxx (temporary description)",
                "level":        "L2S",
                "release":      {
                    "author":       "Erik P G Johansson",
                    "date":         "2016-05-19",
                    "version":      "01",
                    "contact":      "erik.johansson@irfu.se",
                    "institute":    "IRF-U",
                    "modification": "None (initial release)"
                }
            }
        }
    },
    {
        "name":         "testmode2",
        "purpose":      "Mode 2 (temporary purpose description)",
        "inputs":       {
            "input_SCI":    {
                "identifier":   "ROC-SGSE_L2R_RPW-TDS-SURV-RSWF",
                "version":      "01"
            }
        },
        "outputs":      {
            "output_SCI":   {
                "identifier":   "ROC-SGSE_L2S_RPW-BIA-xxxxx",
                "name":         "xxxxx (temporary name)",
                "description":  "Contains xxxxx (temporary description)",
                "level":        "L2S",
                "release":      {
                    "author":       "Erik P G Johansson",
                    "date":         "2016-05-19",
                    "version":      "01",
                    "contact":      "erik.johansson@irfu.se",
                    "institute":    "IRF-U",
                    "modification": "None (initial release)"
                }
            }
        }
    }
    ]
}
