import json
import pickle
import urllib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon

samples = {
    "Nonhomologous end-joining": ["XRCC6",
                                  "XRCC5",
                                  "PRKDC",
                                  "LIG4",
                                  "XRCC4",
                                  "DCLRE1C",
                                  "NHEJ1"],

    "Microhomology end-joining": ["MRE11",
                                  "RAD50",
                                  "NBN",
                                  "RBBP8",
                                  "ERCC4",
                                  "ERCC1",
                                  "LIG1",
                                  "POLL",
                                  "POLB",
                                  "PARP1",
                                  "LIG3",
                                  "XRCC1"
                                  ],
    "Homologous recombination": ["RAD51",
                                 "RAD51B",
                                 "RAD51D",
                                 "DMC1",
                                 "XRCC2",
                                 "XRCC3",
                                 "RAD52",
                                 "RAD54L",
                                 "RAD54B",
                                 "BRCA1",
                                 "SHFM1",
                                 "RAD50",
                                 "MRE11",
                                 "NBN",
                                 "RBBP8",
                                 "MUS81",
                                 "EME1",
                                 "EME2",
                                 "GIYD1",
                                 "GIYD2",
                                 "GEN1",
                                 ],
    "Base excision repair": ["UNG",
                             "SMUG1",
                             "MBD4",
                             "TDG",
                             "OGG1",
                             "MUTYH",
                             "NTHL1",
                             "MPG",
                             "NEIL1",
                             "NEIL2",
                             "NEIL3",
                             "APEX1",
                             "APEX2",
                             "LIG3",
                             "XRCC1",
                             "PNKP",
                             "APLF",
                             "PARP1",
                             "PARP2",
                             "PARP3",
                             "MGMT",
                             "ALKBH2",
                             "ALKBH3",
                             ],
    "Repair of DNA-topoisomerase crosslinks": ["TDP1",
                                               "TDP2"
                                               ],
    "Mismatch excision repair": ["MSH2",
                                 "MSH3",
                                 "MSH6",
                                 "MLH1",
                                 "PMS2",
                                 "MSH4",
                                 "MSH5",
                                 "MLH3",
                                 "PMS1",
                                 "PMS2L3"
                                 ],
    "Nucleotide excision repair": ["RAD23B",
                                   "CETN2",
                                   "RAD23A",
                                   "XPA",
                                   "DDB1",
                                   "DDB2",
                                   "RPA1",
                                   "RPA2",
                                   "RPA3",
                                   "TFIIH",
                                   "ERCC3",
                                   "ERCC2",
                                   "GTF2H1",
                                   "GTF2H2",
                                   "GTF2H3",
                                   "GTF2H4",
                                   "GTF2H5",
                                   "CDK7",
                                   "CCNH",
                                   "MNAT1",
                                   "ERCC5",
                                   "ERCC1",
                                   "ERCC4",
                                   "LIG1",
                                   "ERCC8",
                                   "ERCC6",
                                   "UVSSA",
                                   "XAB2",
                                   "MMS19",
                                   ],
    "Fanconi anemia": ["FANCA",
                       "FANCB",
                       "FANCC",
                       "BRCA2",
                       "FANCD2",
                       "FANCE",
                       "FANCF",
                       "FANCG",
                       "FANCI",
                       "BRIP1",
                       "FANCL",
                       "FANCM",
                       "PALB2",
                       "RAD51C",
                       "BTBD12",
                       "FAAP20",
                       "FAAP24",
                       ],
    "Modulation of nucleotide pools": ["NUDT1",
                                       "DUT",
                                       "RRM2B",
                                       ],
    "DNA polymerases": ["POLB",
                        "POLG",
                        "POLD1",
                        "POLE",
                        "PCNA",
                        "REV3L",
                        "MAD2L2",
                        "REV1L",
                        "POLH",
                        "POLI",
                        "POLQ",
                        "POLK",
                        "POLL",
                        "POLM",
                        "POLN",
                        ],
    "Editing and processing nucleases": ["FEN1",
                                         "FAN1",
                                         "TREX1",
                                         "TREX2",
                                         "EXO1",
                                         "APTX",
                                         "SPO11",
                                         "ENDOV",
                                         ],
    "Ubiquitination and modification": ["UBE2A",
                                        "UBE2B",
                                        "RAD18",
                                        "SHPRH",
                                        "HLTF",
                                        "RNF168",
                                        "SPRTN",
                                        "RNF8",
                                        "RNF4",
                                        "UBE2V2",
                                        "UBE2N",
                                        ],
    "Chromatin Structure and Modification": ["H2AFX",
                                             "CHAF1A",
                                             "SETMAR",
                                             ],
    "Other conserved DNA damage response genes": ["ATR",
                                                  "ATRIP",
                                                  "MDC1",
                                                  "RAD1",
                                                  "RAD9A",
                                                  "HUS1",
                                                  "RAD17",
                                                  "CHEK1",
                                                  "CHEK2",
                                                  "TP53",
                                                  "TP53BP1",
                                                  "RIF1",
                                                  "TOPBP1",
                                                  "CLK2",
                                                  "PER1",
                                                  ]
}


def diff(li1, li2):
    return (list(set(li1) - set(li2)))


def save_obj(obj, name):
    with open('obj/' + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


def getStartandEnd(symb):
    try:
        symbolResponse = json.load(urllib.request.urlopen(
            f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symb}?content-type=application/json;expand=1"))

        chromosome = int(symbolResponse['seq_region_name'])
        if symbolResponse['assembly_name'] == 'GRCh37':
            start = plotBuilder.lift(symbolResponse['start'], chromosome)
            end = plotBuilder.lift(symbolResponse['end'], chromosome)
        else:
            start = symbolResponse['start']
            end = symbolResponse['end']
        return chromosome, start, end
    except Exception as e:
        print(e)
        print(
            f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symb}?content-type=application/json;expand=1")
        print(f"symb got fricked")
        return None
