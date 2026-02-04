# iCHO3K

<p align="center">
  <a href="https://doi.org/10.1101/2025.04.10.647063">
    <img alt="bioRxiv DOI" src="https://img.shields.io/badge/bioRxiv-10.1101%2F2025.04.10.647063-B31B1B">
  </a>
  <img alt="COBRApy" src="https://img.shields.io/badge/COBRApy-compatible-blue">
  <img alt="MATLAB" src="https://img.shields.io/badge/MATLAB-COBRA%20Toolbox-orange">
  <img alt="Model" src="https://img.shields.io/badge/iCHO3K-model%20repository-brightgreen">
  <img alt="Formats" src="https://img.shields.io/badge/formats-SBML%20%7C%20JSON%20%7C%20MATLAB-informational">
</p>

> A community-driven, genome-scale metabolic reconstruction and analysis toolkit for Chinese Hamster Ovary (CHO) cells.

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="iCHO3K/banner/icho3k_figure.png">
    <source media="(prefers-color-scheme: light)" srcset="iCHO3K/banner/icho3k_figure.png">
    <img src="iCHO3K/banner/icho3k_figure.png" alt="iCHO3K metabolic network overview" width="700">
  </picture>
</p>

**Quick links:** [Installation](#installation--setup) · [Create model](#create-icho3k-model) · [Pre-process](#pre-process-icho3k-for-context-specific-model-generation) · [Citing](#citing)

**Highlights**
- **Scope:** 11,004 reactions • 7,377 metabolites • 3,597 genes • 3,489 mapped protein structures  
- **Use cases:** context-specific modeling (WT vs ZeLa), pFBA/ecFBA, flux sampling, subsystem coverage, structure-guided hypotheses (e.g., putative PEP allosteric inhibition of PFK)  
- **Artifacts:** curated datasets, notebooks (Python & MATLAB), figures, and utilities to reproduce key analyses

**Model provided in this repository**
> **Note:** The model files distributed in `iCHO3K/Model/` correspond to an **unblocked** version of the full iCHO3K network. In this version, blocked reactions have been removed and boundary reactions have been configured to support context-specific model generation with tools such as **mCADRE**.  
> The unblocked iCHO3K model contains **8,370 reactions, 4,695 metabolites, and 2,929 genes**.  
> To generate the **full** iCHO3K model, follow the steps in [Create iCHO3K Model](#create-icho3k-model).

If you use this repository or the iCHO3K model, please see **Citing** below.
```

---

|  |  |
|---|---|
| **A community-driven, genome-scale metabolic reconstruction and analysis toolkit for Chinese Hamster Ovary (CHO) cells.**<br><br>**Quick links:** [Installation](#installation--setup) · [Create model](#create-icho3k-model) · [Pre-process](#pre-process-icho3k-for-context-specific-model-generation) · [Citing](#citing) | <picture><source media="(prefers-color-scheme: dark)" srcset="iCHO3K/banner/icho3k_figure.png"><source media="(prefers-color-scheme: light)" srcset="iCHO3K/banner/icho3k_figure.png"><img src="iCHO3K/banner/icho3k_figure.png" alt="iCHO3K metabolic network overview" width="520"></picture> |

---

## Table of contents

- [Repository layout](#repository-layout)
- [Installation & setup](#installation--setup)
- [Create iCHO3K Model](#create-icho3k-model)
- [Pre-process iCHO3K for context-specific model generation](#pre-process-icho3k-for-context-specific-model-generation)
- [Reproducing key analyses](#reproducing-key-analyses)
- [Data & provenance](#data--provenance)
- [Model formats & I/O](#model-formats--io)
- [Solvers & performance tips](#solvers--performance-tips)
- [Contributing](#contributing)
- [Citing](#citing)
- [Maintainers & contact](#maintainers--contact)
- [Acknowledgments](#acknowledgments)
- [FAQ / Troubleshooting](#faq--troubleshooting)

---

## Repository layout

```
Analyses/
├── conf_score_distribution.png            # Confidence Score distribution across all reactions from iCHO3K
├── data_preprocessing/                    # ZeLa vs WT growth rate and spent media data analysis
├── ecFBA and Flux Sum/                    # ZeLa vs WT ecFBA and Flux Sum datasets
├── growth_rate_pred/                      # pFBA simulations from ZeLa and WT context-specific models
├── recons_comparisons/                    # Plot comparisions between iCHO3K and previous CHO reconstructions
├── Relevant_mets/                         # Analysis of subsystems related to metabolites relevant to biomass sysnthesis
├── sampled_results/
├── subsystem_overview/                    #Subsystem/System classification sunburst plot
└── tSNE/                                  #tSNE embedding analysis

Data/                         
├── Context_specific_models/              # Context-specific ZeLa and WT reduced and ec models (MAT)
│   └── unblocked_ecModel_generic/        # Generic iCHO3K ec model
│
├── GPR_Mapping/                          # Supplementary data for GPR Mapping from Recon3D to iCHO3K
├── Gene_Essentiality/                    # Set of experimentally validated CHO essential genes
├── Metabolites/                          # Supplementary data for metabolites information
├── Orthologs/                            # Ortholog mapping information from Human to Chinese Hamster
├── Reconciliation/                       # Source reconstructions & derived models and datasets
│   ├── datasets/
│   └── models/
│
├── Uptake_Secretion_Rates/               # Pre-processed uptake and secretion rates from ZeLa fed-batch data
└── ZeLa Data/                            # ZeLa 14 fed-batch raw transcriptomics and spent media data

iCHO3K/
├── Dataset/                              # iCHO3K source dataset for model generation
└── Model/                                # iCHO3K generic **unblocked** model variants

Matlab/
├── ecFBA/                               # ecFBA scripts
├── Model_Extraction/                    # Model Extraction with mCADRE scripts
├── Standep/                             # ubiquityScore calculation with Standep
├── main_Model_Extraction.m              # Main code for mCADRE model extraction
└── main_standep.m                       # Main code for Standep data preprocessing

Python/
├── Network Reconstruction/              # Notebooks related to the reconciliation of previous reconstructions and building of iCHO3K
│   ├── Genes.ipynb                      # Retrieval of Gene information from databases
│   ├── Metabolites.ipynb                # Integration of metabolite information, de-duplication and analysis
│   ├── Reactions.ipynb                  # Reconciliation of previous CHO and Recon3D reconstructions, de-duplication, subsytem re-organization 
│   └── retrieveTurnoverNumber.ipynb     # Fetch turnover numbers and molecular weights from the BRENDA
│                 
├── Supplementary Notebooks/             # Supplementary Notebooks with extra information of previous reconstructions
├── Comparison..Reconstructions.ipynb    # Comparison of iCHO3K with previous CHO reconstructions
├── Computational_Tests.ipynb            # pFBA Simulations & tSNE embeddings
├── Model_builder.ipynb                  # iCHO3K Model Buider
├── Calculate_specific_rates.ipynb       # Preprocess of spent media data into GEM fluxes
└── ZeLa_fluxomics.ipynb                 # ZeLa fluxomics data

```

> **Large files:** Some assets may use **Git LFS**. If you see pointer files, run:
> ```bash
> git lfs install
> git lfs pull --include="Data/**"
> ```
> **Optional:** clone without data, then fetch only Data/:
> ```bash
> GIT_LFS_SKIP_SMUDGE=1 git clone <repo-url>
> cd iCHO3K && git lfs install
> git lfs pull --include="Data/**" --exclude=""
> ```


---

## Installation & setup

### Option 1 — Recommended (via environment.yml)
Create the exact conda environment shipped with the repo:

```bash
conda env create -f environment.yml
conda activate icho3k
```

### Option 2 — Pip-only (subset)
If you prefer pip, use the lightweight set (no native solvers or libSBML):
```bash
python -m venv .venv && source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

### Optional: graph utilities & network export
```bash
# conda-forge (recommended)
conda install -c conda-forge graphviz pygraphviz ndex2
# If you use pip for pygraphviz, ensure system Graphviz is installed first.
```

> If you plan to use commercial solvers (e.g., Gurobi/CPLEX), install them separately and set the solver in your code (see **Solvers & performance tips**).


### MATLAB (optional)
- MATLAB R2022b+ recommended (earlier versions likely workable).
- Add `Notebooks/Matlab/` to your MATLAB path.

---

## Create iCHO3K Model
```python
import pandas as pd
import cobra
from cobra import Model, Reaction, Metabolite

#Read iCHO3K Full Dataset (11,004 reactions • 7,377 metabolites • 3,597 genes)
FILE_PATH = 'iCHO3K/Dataset/iCHO3K.xlsx'
metabolites      = pd.read_excel(FILE_PATH, sheet_name='Metabolites')
rxns             = pd.read_excel(FILE_PATH, sheet_name='Rxns')
rxns_attributes  = pd.read_excel(FILE_PATH, sheet_name='Attributes')
boundary_rxns    = pd.read_excel(FILE_PATH, sheet_name='BoundaryRxns')
genes_df         = pd.read_excel(FILE_PATH, sheet_name='Genes')

# Create a cobra model and add reactions
model = Model("iCHO3K_vX")
lr = []
for _, row in rxns.iterrows():
    r = Reaction(row['Reaction'])
    lr.append(r)    
model.add_reactions(lr)

# Add reaction-specific information
for i,r in enumerate(tqdm(model.reactions)):
    r.build_reaction_from_string(rxns['Reaction Formula'][i])
    r.name = rxns['Reaction Name'][i]
    r.subsystem = rxns['Subsystem'][i]
    if not (pd.isna(rxns['GPR_iCHO3K'][i]) or rxns['GPR_iCHO3K'][i] == ''):
        r.gene_reaction_rule = str(rxns['GPR_iCHO3K'][i])
    r.lower_bound = float(rxns_attributes['Lower bound'][i])
    r.upper_bound = float(rxns_attributes['Upper bound'][i])
    r.annotation['confidence_score'] = str(rxns['Conf. Score'][i])
    r.annotation['molwt'] = str(rxns_attributes['Mol wt'][i])
    r.annotation['kcat_f'] = str(rxns_attributes['kcat_forward'][i])
    r.annotation['kcat_b'] = str(rxns_attributes['kcat_backward'][i])

# Add Boundary Reactions
dr = []
for _, row in boundary_rxns.iterrows():
    r = Reaction(row['Reaction'])
    dr.append(r)    
model.add_reactions(dr)

boundary_rxns_dict = boundary_rxns.set_index('Reaction').to_dict()
boundary_rxns_dict

for i,r in enumerate(tqdm(model.reactions)):
    if r in dr:
        r.build_reaction_from_string(boundary_rxns_dict['Reaction Formula'][r.id])
        r.name = boundary_rxns_dict['Reaction Name'][r.id]
        r.subsystem = boundary_rxns_dict['Subsystem'][r.id]
        r.lower_bound = float(boundary_rxns_dict['Lower bound'][r.id])
        r.upper_bound = float(boundary_rxns_dict['Upper bound'][r.id])
        r.annotation['confidence_score'] = str(1)

# Add information for each metabolite
metabolites_dict = metabolites.set_index('BiGG ID').to_dict('dict')
for met in model.metabolites:
    try:
        met.name = metabolites_dict['Name'][f'{met}']
        met.formula = metabolites_dict['Formula'][f'{met}']
        met.compartment = metabolites_dict['Compartment'][f'{met}'].split(' - ')[0]
        try:
            met.charge = int(metabolites_dict['Charge'][f'{met}'])
        except (ValueError, TypeError):
            print(f'{met} doesnt have charge')
    except (KeyError):
        print('----------------------------')
        print(f'{met} doesnt exist in the df')
        print('----------------------------')

#OPTIONAL: Add Gene Name information
genes_dict = genes_df.iloc[:,:2].set_index('Gene Entrez ID').to_dict('dict')
for g in model.genes:
    if g.id in list(genes_dict['Gene Symbol'].keys()):
        g.name = genes_dict['Gene Symbol'][f'{g}']

# Set biomass objective function and save the model
name  = 'CHO-Generic'  # OPTIONS: 'CHO-S', 'CHO-Generic', 'CHO-Prod-Generic'
if name == 'CHO-Generic':
    model.objective = 'biomass_cho'
elif name == 'CHO-S':
    model.objective = 'biomass_cho_s'
elif name == 'CHO-Prod-Generic':
    model.objective = 'biomass_cho_prod'
else:
    raise ValueError(f"Unrecognized model name: {name}")

# OPTIONAL: remove infeasible loop-generating reactions:
glucloop_reactions = [
    'GapFill-R01206', 'GAUGE-R00557', 
    'GAUGE-R00558', 'FNOR', 'GGH', 'r0741', 'r1479', 'XOLESTPOOL'
]

try:
    model.remove_reactions(glucloop_reactions, remove_orphans=True)
except KeyError:
    print(f'Reaction {reaction_id} not in model {model.id}')

atp_loop_reactions = [
    'SCP22x','TMNDNCCOAtx','OCCOAtx','r0391','BiGGRxn67','r2247','r2280',
    'r2246','r2279','r2245','r2305','r2317','r2335','HMR_0293','HMR_7741',
    'r0509','r1453','HMR_4343','ACONTm','PDHm','r0426','r0383','r0555',
    'r1393','NICRNS','GAUGE-R00648','GAUGE-R03326','GapFill-R08726','RE2915M',
    'HMR_3288','HMR_1325','HMR_7599','r1431','r1433','RE2439C','r0791',
    'r1450','GAUGE-R00270','GAUGE-R02285','GAUGE-R04283','GAUGE-R06127','GAUGE-R06128',
    'GAUGE-R06238','GAUGE-R00524','RE3477C','AAPSAS','RE3347C','HMR_0960','HMR_0980',
    'RE3476C','r0708','r0777','r0424','r0698','3HDH260p','HMR_3272','ACOAD183n3m',
    'HMR_1996','GapFill-R01463','GapFill-R04807','r1468','r2435','r0655','r0603','r0541',
    'RE0383C','HMR_1329','TYRA','NRPPHRt_2H','GAUGE-R07364','GapFill-R03599','ARD',
    'RE3095C','RE3104C','RE3104R','ACONT','ICDHxm','ICDHy',
    'r0425','r0556','NH4t4r','PROPAT4te','r0085','r0156','r0464','ABUTDm',
    'OIVD1m','OIVD2m','OIVD3m','r2194','r2202','HMR_9617','r2197','r2195',
    '2OXOADOXm','r2328','r0386','r0451','FAS100COA','FAS120COA','FAS140COA',
    'FAS80COA_L','r0604','r0670','r2334','r0193','r0595','r0795','GLYCLm',
    'MACACI','r2193','r0779','r0669','UDCHOLt','r2146','r2139'
]

try:
    model.remove_reactions(atp_loop_reactions, remove_orphans=True)
except KeyError:
    print(f'Reaction {reaction_id} not in model {model.id}')


##### ----- Save iCHO3K full version (11,004 reactions • 7,377 metabolites • 3,597 genes) ----- #####

suffix = name.lower().replace('-', '_')

# XML
model_name_xml = f'../iCHO3K/Model/iCHO3K_{suffix}.xml' 
write_sbml_model(model, model_name_xml)

# JSON, because the sbml doesnt save the subsystems
model_name_json = f'../iCHO3K/Model/iCHO3K_{suffix}.json' 
save_json_model(model, model_name_json)

# MATLAB
model_name_matlab = f'../iCHO3K/Model/iCHO3K_{suffix}.mat' 
save_matlab_model(model, model_name_matlab)

```
---

## Pre-process iCHO3K for context-specific model generation
```python

from cobra.flux_analysis import find_blocked_reactions, flux_variability_analysis

## Step 1: Identify blocked reactions

for rxn in model.boundary:
    if rxn.id.startswith("EX_"):
        rxn.bounds = (-1000,1000)
    if rxn.id.startswith("SK_"):
        rxn.bounds = (-1000,1000)
    if rxn.id.startswith("DM_"):
        rxn.bounds = (0,1000)

model.solver = 'gurobi'
blocked_reactions = find_blocked_reactions(model)

## Step 2: Remove blocked reactions from the model

model_unblocked = model.copy()
blocked_reaction_objects = [model_unblocked.reactions.get_by_id(rxn_id) for rxn_id in blocked_reactions]
model_unblocked.remove_reactions(blocked_reaction_objects, remove_orphans=True)
print(f"Removed {len(blocked_reaction_objects)} blocked reactions from the model.")

## Step 3: Set the bounds for context-specific model generation

for rxn in model.boundary:
        
    # Models that are forced to secrete ethanol are not feasible
    if rxn.id == 'EX_etoh_e':
        rxn.bounds = (-1,1)
        continue
    
    # Keep boundaries open for essential metabolites
    if rxn.id == 'EX_h2o_e':
        rxn.bounds = (-1000,1000)
        continue
    if rxn.id == 'EX_h_e':
        rxn.bounds = (-1000,1000)
        continue
    if rxn.id == 'EX_o2_e':
        rxn.bounds = (-1000,1000)
        continue
    if rxn.id == 'EX_hco3_e':
        rxn.bounds = (-1000,1000)
        continue
    if rxn.id == 'EX_so4_e':
        rxn.bounds = (-1000,1000)
        continue
    if rxn.id == 'EX_pi_e':
        rxn.bounds = (-1000,1000)
        continue

    # Boundaries from Sink reactions on iCHO_v1 (100 times lower)
    if rxn.id == 'SK_Asn_X_Ser_Thr_r':
        rxn.bounds = (-0.001,1000)
        continue
    if rxn.id == 'SK_Tyr_ggn_c':
        rxn.bounds = (-0.001,1000)
        continue
    if rxn.id == 'SK_Ser_Thr_g':
        rxn.bounds = (-0.001,1000)
        continue
    if rxn.id == 'SK_pre_prot_r':
        rxn.bounds = (-0.001,1000)
        continue
    
    # Close uptake rates for the rest of the boundaries
    if rxn.id.startswith("EX_"):
        rxn.bounds = (0,1000) 
    if rxn.id.startswith("SK_"):
        rxn.bounds = (0,1000)
    if rxn.id.startswith("DM_"):
        rxn.bounds = (0,1000)

##### ----- Save iCHO3K unblocked version ----- #####

# XML
model_name_xml = f'../iCHO3K/Model/iCHO3K_{suffix}_unblocked.xml' 
write_sbml_model(model_unblocked, model_name_xml)

# JSON, because the sbml doesnt save the subsystems
model_name_json = f'../iCHO3K/Model/iCHO3K_{suffix}_unblocked.json' 
save_json_model(model_unblocked, model_name_json)

# MATLAB
model_name_matlab = f'../iCHO3K/Model/iCHO3K_{suffix}_unblocked.mat' 
save_matlab_model(model_unblocked, model_name_matlab)

```

## Reproducing key analyses

Many figures in **Analyses/** and the accompanying manuscript (https://doi.org/10.1101/2025.04.10.647063) were generated from notebooks in **Python/**:

- **Reconstruction comparisons and Gene Essentiality Analysis** →  
  *Python/*`Comparison of Metabolic Reconstructions.ipynb` → *Analyses/recons_comparisons/*
- **Subsystem coverage & sunbursts** →  
  *Python/*`Comparison of Computational Tests.ipynb /1. Subsystem Overview and Analysis` → *Analyses/subsystem_overview/*
- **Topology & t-SNE** →  
  *Python/*`Comparison of Computational Tests.ipynb /4. tSNE Comparison of Models` → *Analyses/tSNE/*
- **Growth rate prediction (WT vs ZeLa)** →  
  *Python/*`Comparison of Computational Tests.ipynb /5. pFBA` → *Analyses/growth_rate_pred/*

> Most notebooks begin with a “Paths & Environment” cell—update paths as needed. For strict reproducibility, pin exact package versions via `environment.yml` and use releases/DOI snapshots.

---

## Data & provenance

Curated inputs and derived artifacts are organized under **Data/**. Key elements:

- **WT & ZeLa context-specific models** → `Data/Context_specific_models/`
- **Experimentally validated CHO essential genes** → `Data/Gene_Essentiality/`  
- **GPR from Recon3D to CHO mapping** → `Data/GPR_Mapping/`   
- **Metabolite Curation Files** → `Data/Metabolites/`  
- **Soruce Reconstructions Datasets and Models** → `Data/Metabolites/`  
- **ZeLa vs WT Fed-batch Process Raw Data** → `Data/ZeLa Data/`
- **ZeLa vs WT Fed-batch Process pre-processed Data** → `Data/Uptake_Secretion_Rates/`
  
---

## Model formats & I/O
> The models in `iCHO3K/Model/` are the **unblocked** iCHO3K variants (blocked reactions removed; boundaries configured for context-specific extraction such as mCADRE).  
> To generate the **full** iCHO3K model directly from the Excel dataset, see [Create iCHO3K Model](#create-icho3k-model).
- **Excel**: Final iCHO3K set of reactions, genes and metabolites lives in `iCHO3K/Dataset` for inspection and conversion.
- **COBRApy (JSON / SBML)**
  ```python
  import cobra

  # JSON
  m_json = cobra.io.load_json_model("iCHO3K/Model/<model_name>.json")

  # SBML
  m_sbml = cobra.io.read_sbml_model("iCHO3K/Model/<model_name>.xml")
  ```
---

## Solvers & performance tips

- **LP/QP solvers:** GLPK (free), CPLEX/Gurobi (commercial, academic licenses). Set via COBRApy:
  ```python
  import cobra
  cobra.Configuration().solver = "glpk"  # or "gurobi", "cplex"
  ```
- **Speed:** Prefer commercial solvers for large sampling tasks; reduce model size using context-specific models; cache solutions where possible.
- **Numerics:** Tighten feasibility/optimality tolerances for sensitive analyses.

---

## Contributing

Contributions are welcome!

1. **Issues**: Report bugs, request features, or flag data discrepancies.  
2. **PRs**: Use feature branches; include a clear description, minimal reproducible example or notebook, and updated docs.  
3. **Style**: PEP 8 for Python; strip heavy notebook outputs before committing.  
4. **Data**: Avoid committing large binaries; use **Git LFS** or attach to Releases/Zenodo.

If contributing new datasets or model variants, please include:
- Data dictionary (column descriptions, units)
- Provenance (source links/versions)
- Minimal script/notebook to regenerate derived artifacts

---

## Citing

If you use **iCHO3K** or materials from this repository, please cite the **bioRxiv preprint**:

Di Giusto, P., Choi, D.-H., *et al.* (2025). *A community-consensus reconstruction of Chinese Hamster metabolism enables structural systems biology analyses to decipher metabolic rewiring in lactate-free CHO cells.* **bioRxiv**. https://doi.org/10.1101/2025.04.10.647063 (v1 posted April 17, 2025).

---

## Maintainers & contact

- **Pablo Di Giusto** — pdigiusto@health.ucsd.edu · pablodigiusto91@gmail.com  
  Systems Biology & Cell Engineering Lab (Lewis Lab), University of Georgia

---

## Acknowledgments

We thank the iCHO3K community contributors and collaborators. This work builds upon public resources: iCHO1766, iCHO2101, iCHO2291, iCHO2441, Recon3D, Human-GEM, BiGG, MetaNetX, Rhea, UniProt, BRENDA, and others referenced throughout the notebooks and the manuscript.

---

## FAQ / Troubleshooting

**I see *.gitattributes* LFS pointers instead of files.**  
Run:
```bash
git lfs install
git lfs pull
```

**Solver not found / poor performance.**  
Install an LP solver (GLPK works; Gurobi/CPLEX recommended for speed). Set `cobra.Configuration().solver = "gurobi"` once installed.

**Model won’t optimize (infeasible).**  
- Harmonize exchange bounds across conditions.
- Check blocked reactions / dead-end metabolites.
- For comparative runs (WT vs ZeLa), ensure identical media constraints.

**Notebook paths are wrong.**  
Edit the first “Paths & Environment” cell—most notebooks expose a single place to set root paths.
