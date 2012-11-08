"""This file contains lists of organisms"""

test_orgs = ["Eschericia_coli",
             "Pseudomonas_aeruginosa",
             "Haemophilus_influenzae"]

orgs = ["Acinetobacter_ADP1",
        "Pseudomonas_fluorescens", 
        "Acinetobacter_baumannii", 
        "Pseudomonas_fulva", 
        "Azotobacter_vinelandii", 
        "Pseudomonas_mendocina", 
        "Cellvibrio__gilvus", 
        "Pseudomonas_putida", 
        "Cellvibrio_japonicus", 
        "Pseudomonas_stutzeri", 
        "Pseudomonas_syringae", 
        "Moraxella_catarrhalis", 
        "Psychrobacter_arcticus_273", 
        "Pseudomonas_aeruginosa",  
        "Psychrobacter_cryohalolentis", 
        "Pseudomonas_brassicacearum", 
        "Psychrobacter_PRwf", 
        "Pseudomonas_entomophila"]

gammas = ["Haemophilus_influenzae", 
          "Escherichia_coli", 
          "Salmonella_typhimurium", 
          "Shewanella_oneidensis", 
          "Vibrio_cholerae"
          ]

firmicutes = ["Bacillus_subtilis", 
              "Lactobacillus_plantarum", 
              "Lactococcus_lactis", 
              "Bacillus_halodurans", 
              "Enterococcus_faecalis", 
              ]
#
actinos = ["Arthrobacter", 
           "Mycobacterium_smegmatis", 
           "Streptomyces_coelicolor", 
           "Bifidobacterium_breve", 
           "Corynebacterium_diphtheriae"]
misc = ["Helicobacter_pylori_HPAG1",
        ]
enteros = ["Enterococcus_faecalis",
           "Bifidobacterium_breve",
           "Bacteroides_thetaiotaomicron",
           "Fusobacterium_nucleatum",
           "Lactobacillus_gasseri",
           "Proteus_mirabilis",
           "Escherichia_coli",
           "Eubacterium_rectale"]

acinetos = ["Acinetobacter_baumannii",
               "Acinetobacter_ADP1",
               "Acinetobacter_calcoaceticus",
               "Acinetobacter_oleivorans"]

pseudos = [org for org in orgs if "Pseudo" in org]
psychros = [org for org in orgs if "Psychro" in org]
all_orgs = list(set(orgs + gammas + firmicutes + actinos + enteros))
new_orgs = gammas + firmicutes
pvp_orgs = pseudos + psychros
last_orgs = pvp_orgs + new_orgs + actinos
ecoli = "Escherichia_coli"
group_names = ["pseudos","psychros","pvp_orgs","gammas",
               "firmicutes","enteros","actinos","last_orgs"]

def org_group(org):
    # lookup = {"pseudos":pseudos,
    #           "psychros":psychros,
    #           "gammas":gammas,
    #           "firmicutes":firmicutes,
    #           "actinos":actinos,
    #           "enteros":enteros,
    #           "all_orgs":all_orgs
    #           }
    for group in ["pseudos",
              "psychros",
              "gammas",
              "firmicutes",
              "actinos",
              "enteros",
              "all_orgs"]:
        if org in lookup[group]:
            return group
    assert False, "couldn't find %s in lookup table" % org
    
del(org) # namespace leak
