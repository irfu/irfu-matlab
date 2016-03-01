{
    "identification":   {
        "project":  "ROC-SGSE",
        "name": "BIASPSW (temporary name)",
        "identifier":   "ROC-SGSE-BIASPSW",
        "description":  "BIAS pipeline software (temporary description)"
    },  
    "release":  {   
        "version":  "1.0.0",
        "date": "2016-03-01",   
        "author":   "Erik Johansson",
        "contact":  "erik.johansson@irfu.se", 
        "institute":    "IRF-U",    
        "modification": "None (Initial release)"
    },  
    "environment":  {   
        "activation":   "scripts/setup.sh",
        "deactivation": "scripts/unset.sh"   
        "executable":   "bin/biassw",
    }
    "modes":    [
    {
        "name":     "mode_1",
        "purpose":  "Mode 1 (temporary purpose description)",
        "inputs":   {
            "input_1": {   
                "identifier":   "ROC-SGSE_LL01_RPW-BIA",
                "version":      "01"
            }   
        },  
        "outputs":  {
            "output_2":   {   
                "identifier":   "ROC-SGSE_L2S_RPW-BIA-SWEEP",
                "name":         "ROC SGSE RPW BIAS L2S Sweep data",
                "description":  "Contains RPW BIAS L2S Sweep data (temporary description)",
                "level":    "L2S",
                "release":  {   
                    "author":   "Erik P G Johansson",
                    "date":     "2016-03-01",
                    "version":  "01",
                    "contact":  "erik.johansson@irfu.se",
                    "institute":    "IRF-U",
                    "modification": "None (initial release)"
                }
            }
        }
    }
    