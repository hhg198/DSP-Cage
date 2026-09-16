# DSP-Cage

DSP-Cage is a Python-based computational tool for calculating and visualizing the cavity volume of supramolecular cages using a dynamic probe-splitting strategy.

## Installation and Usage

1. Ensure that Python is installed on your system.

   Python 3.9 or later is recommended.

2. Install the required Python packages.

```bash
# Skip packages that are already installed in your Python environment.
pip install numpy
pip install scipy
pip install pymeshlab
pip install biopython
pip install streamlit
```

Alternatively, install all dependencies with one command:

```bash
pip install numpy scipy pymeshlab biopython streamlit
```

3. Place the supramolecular cage structure file in the `datas` folder.

4. Open `Main_Process.py` and modify the `filename` parameter so that it matches the name of the cage structure file.

```python
filename = "datas/your_cage_file.pdb"
```

5. Run the main program.

```bash
python Main_Process.py
```

## Run the Web Interface

DSP-Cage also provides a Streamlit-based visual interface.

Run the following command in a terminal:

```bash
streamlit run D:\daima\DSP_Cage\Visual_APP_En.py
```

> **Important:** `D:\daima\DSP_Cage\Visual_APP_En.py` is only an example path.  
> Replace it with the actual path of `Visual_APP_En.py` on your own computer.

For example:

```bash
streamlit run "C:\Users\YourName\Documents\DSP_Cage\Visual_APP_En2.py"
```

You can also first enter the DSP-Cage project directory and then run:

```bash
cd D:\daima\DSP_Cage
streamlit run Visual_APP_En.py
```

Again, replace `D:\daima\DSP_Cage` with the actual DSP-Cage project path on your own computer.

After the command is executed successfully, Streamlit will normally open the DSP-Cage interface automatically in your default web browser.

## Notes

- Make sure the input file path and file name are correct.
- Run the command from the root directory of the DSP-Cage project.
- The calculation results and visualization files will be generated according to the output settings defined in `Main_Process.py`.
- If the `streamlit` command is not recognized, install Streamlit with:

```bash
pip install streamlit
```

- On some systems, the interface can also be started with:

```bash
python -m streamlit run Visual_APP_En.py
```

## Project Structure

```text
DSP-Cage/
├── Main_Process.py
├── Cavity_Calculation.py
├── Balloon.py
├── Visual_APP_En.py
├── datas/
│   └── your_cage_file.pdb
└── README.md
```

# Acknowledgement and Third-Party Components

DSP-Cage is an independent framework developed for supramolecular cage cavity
characterization based on the Dynamic Splitting of Probes (DSP) strategy.

The core methodology of DSP-Cage, including:

- Dynamic Splitting of Probes strategy;
- adaptive probe sampling mechanism;
- dynamic cavity exploration;
- cavity reconstruction workflow;
- GUI-oriented visualization and interaction modules;

was independently developed in this work.

---

## Relation to CageCavityCalc (C3)

DSP-Cage incorporates several low-level computational components and scientific
parameters from the open-source CageCavityCalc (C3) framework developed by:

**Vicente Marti-Centelles**

Repository:

https://github.com/VicenteMartiCentelles/CageCavityCalc


The original CageCavityCalc software is distributed under the MIT License.
The original copyright notice and license information have been retained.

Please refer to:

LICENSE_C3
and
THIRD_PARTY_LICENSES.md


Users should also cite the original CageCavityCalc work:
CageCavityCalc
Vicente Marti-Centelles et al.
DOI: 10.1021/acs.jcim.4c00355


---

## License

DSP-Cage is released for academic research and non-commercial use.

Third-party components derived from CageCavityCalc remain under their
original MIT License.

Commercial use of DSP-Cage requires permission from the authors.

