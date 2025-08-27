# iCHO3K — CHO Whole-Cell Network Repository

> A community-driven, systems-scale reconstruction and analysis toolkit for Chinese Hamster Ovary (CHO) cells.

**Highlights**
- **Scope:** 11,004 reactions • 7,377 metabolites • 3,597 genes • 3,489 mapped protein structures  
- **Use cases:** context-specific modeling (WT vs ZeLa), pFBA/ecFBA, flux sampling, subsystem coverage, structure-guided hypotheses (e.g., putative PEP allosteric inhibition of PFK)  
- **Artifacts:** curated datasets, notebooks (Python & MATLAB), figures, and utilities to reproduce key analyses

> If you use this repository or the iCHO3K model, please see **Citing** below. A manuscript is in preparation; formal citation details will be added upon release.

---

## Table of contents

- [Repository layout](#repository-layout)
- [Installation & setup](#installation--setup)
- [Quickstart (Python)](#quickstart-python)
- [Quickstart (MATLAB)](#quickstart-matlab)
- [Reproducing key analyses](#reproducing-key-analyses)
- [Data & provenance](#data--provenance)
- [Model formats & I/O](#model-formats--io)
- [Solvers & performance tips](#solvers--performance-tips)
- [Contributing](#contributing)
- [Citing](#citing)
- [License](#license)
- [Maintainers & contact](#maintainers--contact)
- [Acknowledgments](#acknowledgments)
- [FAQ / Troubleshooting](#faq--troubleshooting)

---

## Repository layout

```
Analyses/                     # Results & figures from analyses
├── conf_score_distribution.png
├── flux_enrichment_analysis/
├── growth_rate_pred/
├── recons_comparisons/
├── sampled_results/
├── subsystem_overview/
└── tSNE/

Data/                         # Curated datasets & model assets
├── Context_specific_models/  # Context-specific models (MAT, JSON) + scores
├── GPR_Curation/
├── Gene_Essentiality/
├── Metabolites/
├── Orthologs/
├── Reconciliation/           # source reconstructions & derived models
│   ├── datasets/
│   └── models/
├── Sec_Recon_shared_genes/
├── Subsystem/
├── Uptake_Secretion_Rates/
├── kcat_values/
├── iCHO3K_final/             # finalized iCHO3K (Excel format)
├── sampling_ADSB/
├── ZeLa Data/
└── Supplementary Data.xlsx, mart_GO.txt, ncbi_dataset.tsv

Networks/                     # Static images illustrating mass/flux flow

Notebooks/                    # Jupyter & MATLAB notebooks/scripts
├── Matlab/                   # extraction, flux sampling, context-specific
├── Supplementary Notebooks/
├── Website Data/
├── metabolite_identifiers.py
├── google_sheet.py
├── utils.py
└── *Comparison of Metabolic Reconstructions.ipynb,* *Final CHO Model.ipynb*, etc.

Other top-level files
├── CHO_Genome.xlsx
└── Notebooks.ipynb           # overview notebook linking to others
```

> **Large files:** Some assets may use **Git LFS**. If you see pointer files, run:
> ```bash
> git lfs install && git lfs pull
> ```

---

## Installation & setup

### Recommended: conda environment

```bash
conda create -n icho3k python=3.11 -y
conda activate icho3k
pip install cobra pandas numpy scipy matplotlib optlang networkx jupyterlab escher seaborn
```

Optional (for graph utilities & network export):
```bash
pip install ndex2 pygraphviz
```

> If `environment.yml` or `requirements.txt` is provided in the repo or a release, prefer installing from those for exact reproducibility.

### MATLAB (optional)
- MATLAB R2022b+ recommended (earlier versions likely workable).
- Add `Notebooks/Matlab/` to your MATLAB path.

---

## Quickstart (Python)

### 1) Load a context-specific model and run pFBA
```python
import cobra
from cobra.flux_analysis import pfba

model = cobra.io.load_json_model("Data/Context_specific_models/ZeLa_model.json")
solution = pfba(model)
print(f"Objective ({model.objective.direction}): {solution.objective_value:.4f}")

# Top 10 absolute fluxes
top = sorted(solution.fluxes.items(), key=lambda x: abs(x[1]), reverse=True)[:10]
for rxn, v in top:
    print(f"{rxn:25s} {v:10.3f}")
```

### 2) Compare WT vs ZeLa growth under the same media
```python
import cobra, pandas as pd

wt   = cobra.io.load_json_model("Data/Context_specific_models/WT_model.json")
zela = cobra.io.load_json_model("Data/Context_specific_models/ZeLa_model.json")

# Example: harmonize key exchange bounds
for ex in ["EX_glc__D_e", "EX_gln__L_e", "EX_o2_e"]:
    for m in (wt, zela):
        if ex in m.reactions:
            m.reactions.get_by_id(ex).lower_bound = -10.0

res = []
for name, m in [("WT", wt), ("ZeLa", zela)]:
    sol = m.optimize()
    res.append({"model": name, "mu": sol.objective_value})

print(pd.DataFrame(res))
```

### 3) Flux sampling (optlang-compatible solver required)
```python
from cobra.sampling import sample
model = cobra.io.load_json_model("Data/Context_specific_models/WT_model.json")
samples = sample(model, n=1000)   # DataFrame
samples.to_csv("Analyses/sampled_results/wt_samples.csv", index=False)
```

---

## Quickstart (MATLAB)

```matlab
% Ensure COBRA Toolbox is installed & initialized
initCobraToolbox(false)  % without updates
changeCobraSolver('glpk', 'LP');

% Load a JSON context-specific model
model = readCbModel('Data/Context_specific_models/WT_model.json');

% Optimize & print objective
solution = optimizeCbModel(model);
fprintf('Growth objective: %.4f\n', solution.f);
```

> MATLAB scripts for extraction, flux sampling, and context-specific modeling are under `Notebooks/Matlab/`.

---

## Reproducing key analyses

Many figures in **Analyses/** are generated from notebooks in **Notebooks/**:

- **Reconstruction comparisons** →  
  *Notebooks/*`Comparison of Metabolic Reconstructions.ipynb` → *Analyses/recons_comparisons/*
- **Subsystem coverage & sunbursts** →  
  *Analyses/subsystem_overview/*
- **Flux enrichment & sampling** →  
  *Analyses/flux_enrichment_analysis/*, *Analyses/sampled_results/*
- **Growth rate prediction (WT vs ZeLa)** →  
  *Analyses/growth_rate_pred/*
- **Topology & t-SNE** →  
  *Analyses/tSNE/*

> Most notebooks begin with a “Paths & Environment” cell—update paths as needed. For strict reproducibility, pin exact package versions via `environment.yml` and use releases/DOI snapshots.

---

## Data & provenance

Curated inputs and derived artifacts are organized under **Data/**. Key elements:

- **Source reconstructions** → `Reconciliation/datasets/` (inputs) and `Reconciliation/models/` (intermediate models).  
- **Annotations & mappings** → `Metabolites/`, `Subsystem/`, `Orthologs/`.  
- **Evidence & curation** → `GPR_Curation/`, `Gene_Essentiality/`, `kcat_values/`.  
- **Experimental constraints** → `Uptake_Secretion_Rates/`, `ZeLa Data/`.  
- **Secretory overlap** → `Sec_Recon_shared_genes/`.  
- **Final model** → `iCHO3K_final/` (Excel format; conversion notebooks provided).

During manual curation, compartment and subsystem information were inherited from source reconstructions; discrepancies were resolved using authoritative resources (see notes within notebooks).

---

## Model formats & I/O

- **Excel**: Final iCHO3K lives in `Data/iCHO3K_final/` for inspection and conversion.
- **SBML / JSON**: Preferred for simulation. Use conversion notebooks (e.g., *Notebooks/Final CHO Model.ipynb*) or COBRApy I/O:
  ```python
  import cobra
  m = cobra.io.load_json_model("path/to/model.json")
  cobra.io.save_json_model(m, "out.json")
  cobra.io.write_sbml_model(m, "out.xml")
  ```

> Some scripts expect standardized BiGG-style IDs. See `Notebooks/metabolite_identifiers.py` for mapping helpers.

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

Until a formal publication is available, please cite this repository:

```
Di Giusto P., Richelle A., Lewis N.E., et al. (2025).
iCHO3K — CHO Whole-Cell Network Reconstruction. GitHub repository.
https://github.com/<org-or-user>/<repo>
```

A `CITATION.cff` will be provided upon manuscript/publication release. If you reuse specific datasets or scripts, also cite the upstream resources referenced in the notebooks (e.g., Recon3D, BiGG, MetaNetX, Rhea, UniProt, BRENDA).

---

## License

See the `LICENSE` file in this repository for terms of use.  
If no license is present, usage defaults to **“all rights reserved”** until a license is added.

---

## Maintainers & contact

- **Pablo Di Giusto** — pdigiusto@health.ucsd.edu · pablodigiusto91@gmail.com  
  Systems Biology & Cell Engineering Lab (Lewis Lab), UC San Diego & University of Georgia

---

## Acknowledgments

We thank the iCHO3K community contributors and collaborators (including secRecon curators). This work builds upon public resources: Recon3D, BiGG, MetaNetX, Rhea, UniProt, BRENDA, and others referenced throughout the notebooks.

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
