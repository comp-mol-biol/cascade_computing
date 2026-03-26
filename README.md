# README: Cascade computing - contact analysis

---

## 1. Environment Setup & Run

### Create a Virtual Environment

First, create an isolated environment. Replace `yourname` with your preferred environment name (e.g., `sim_env`).

```bash
python -m venv yourname

```

### Activate the Environment

Activate the environment:

**Bash:**

```bash
source ./yourname/bin/activate

```

### Install Dependencies

Install the required packages using the provided requirements files. Ensure you provide the correct path to your files.

**For Cascade Computing (CC):**

```bash
pip install -r ./cascade_computing/requirements_cc.txt

```

### Create IPython Kernel 
Create a kernel that you can use inside Jupyter Lab to use your enviroment

```bash
python -m ipykernel install --user --name=yourname --display-name "cc_kernel"

```

### Run Jupyter lab
Jupyter lab runs in the local browser

```bash
jupyer lab &

```
or you can direct it to a port on your HPC system

```bash
nohup jupyter notebook --no-browser --port=8008 &
```
in the second case you need to log onto your HPC system with the same port: `ssh -L 8008:localhost:8008 yourHPC`.

Check if its running using:
```bash
jupyer list

```


## 2. Explore examples

### MUT16-FFR example

Download the corresponding dataset (MUT16_FFR_atm_analysis_part.zip) from `https://zenodo.org/uploads/19239689` and unpack the folder into the examples directory - then open & run: 
```bash
examples/minimal_example_MUT16_FFR.ipnyb
```

You might have to set the path to your cascade_computing directory.
The examples shows you what the results of the analysis can look like:

-what does the contact record include
-how are the contact frequencies and lifetimes saved
-how to evaluate the corresponding distributions

### Tests
Download the corresponding dataset (data.zip) from `https://zenodo.org/uploads/19239689 and unpack the folder into the examples directory.  There are two test scripts included, that you can use to get insight into the contact calculations.
The first one in
```bash
tests/unit_test.ipnyb
```
shows contact calculations and runs a small example that calculated also specific bond interaction.
The second one in 
```bash
tests/debug_test.ipnyb
```
allows you to verify agains pre-computed data, in case you modify or extend your version of the code. 


## 3. Run analysis on your data

Open 

```bash
`setup_analysis.ipynb`
```

then start by configuring your local paths:

1. **Repository Path**: Set `path_git` to your local clone of the `cascade_computing` repository.
2. **Project Name**: Set `p_name` to a unique identifier (e.g., `'MUT16_MUT8'`). This creates a dedicated workspace for your project.

---

### System Requirements & Configuration

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
