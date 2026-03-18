# README: Simulation Project Setup with Signac

---

## 1. Environment Setup

### Create a Virtual Environment

First, create an isolated environment. Replace `yourname` with your preferred environment name (e.g., `sim_env`).

```bash
python -m venv yourname

```

### Activate the Environment

Activate the environment to load changes and use the isolated shell:

**Bash (Linux/macOS):**

```bash
source ./yourname/bin/activate

```

### Install Dependencies

Install the required packages using the provided requirements files. Ensure you provide the correct path to your files.

**For Cascade Computing (CC):**

```bash
pip install -r ./PATH_TO/requirements_cc.txt

```

## 2. Path and Project Initialization

Open `setup_analysis.ipynb` and configure your local paths:

1. **Repository Path**: Set `path_git` to your local clone of the `cascade_computing` repository.
2. **Project Name**: Set `p_name` to a unique identifier (e.g., `'MUT16_MUT8'`). This creates a dedicated workspace for your project.

---

## 3. System Requirements & Configuration

To analyze your system, you must specify the proteins and their structural properties.

### Protein Specification

For every protein in your system, create a configuration dictionary including:

* **`prot`**: The name of the protein.
* **`file`**: Path to the full-length PDB file.
* **`orig_na`**: Total amino acids in the original full-length protein.
* **`cut_na`**: Number of amino acids in the simulated fragment.
* **`min`** / **`max`**: The residue indices in the full-length protein defining your simulated fragment.

### Domain Mapping

If your proteins contain specific domains (e.g., `NTERM`), define them using boundary columns. Use `DOM_MIN` and `DOM_MAX` to define the start and end residues relative to the full sequence.

**Example Configuration:**

```python
domains = ['NTERM', 'FULL']

mut16 = {
    'prot': "MUT16",
    'orig_na': 2204,
    'cut_na': 140,
    'min': 633,
    'max': 772,
    'NTERM_MIN': 633,
    'NTERM_MAX': 700,
    'file': f"{path_input}/{filenameA}.pdb"
}

```

### The df_domains Reference Table

The notebook aggregates these dictionaries into a `df_domains` table. This serves as the primary metadata source, allowing the notebook to automatically extract amino acid sequences and map analysis results to specific domains.

---

## 4. Running Signac Operations

### Load the Project

The notebook connects to the database created during setup. It looks for a `signac.rc` file in your directory to identify the project.

```python
import signac
project = signac.get_project()

```

### Filter and Select Jobs

You can filter for specific proteins or concentrations to limit the scope of your operations:

```python
# Select all jobs
for job in project:
    pass

# Filter for a specific protein
for job in project.find_jobs({'prot': 'MUT16'}):
    pass

```

### Execute Operations

Use the `project.run()` method. This is the preferred execution method as it handles environment variables and logging automatically.

**Example: Running a transformation**

```python
project.run(
    names=['transform'], 
    jobs=[job], 
    progress=True, 
    num_passes=1, 
    order="by-job"
)

```

**Parameter Breakdown:**

* **`names`**: List of operation names to execute (e.g., `['transform', 'contacts']`).
* **`jobs`**: A list containing the specific job object(s) to run.
* **`progress`**: Set to `True` for a progress bar.
* **`order="by-job"`**: Executes all operations for one job before moving to the next.

---

## 6. Available Operations

The project is configured with the following operations:

* **`post_processing`**: Sets up the analysis framework.
* **`transform`**: Trajectory transformation (centering, PBC wrapping).
* **`contacts`**`: Computing residue-residue contact maps.
* **`eval_contacts`**: Contact evaluation and creation of the contact record
* **``analysis**`: Downstream analysis e.g. pivot tables from the contact record
* **`visualization`**: Generating plots .

**Note on Modifications:** To modify these functions, edit the source file located at:
`cascade_computing/src/compute/signac/sgnc.py`

---

## 6. Examples & Test

The notebook cascade_computing/tests/debug_test.ipynb allows you test how to use the contact calculation and evaluation on your system. 
We provide a small trajectory and the resulting contact data of the test run here: https://seafile.rlp.net/seafhttp/f/f448889532704bf98a28/?op=view
After downloading, unpack the data folder into the cascade_computing/examples folder. You can follow through the debug_test.ipynb using this data.
