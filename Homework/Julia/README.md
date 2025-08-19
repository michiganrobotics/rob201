# ROB201-JuliaHW
"Calculus for the Modern Engineer" Julia HW assignments from the University of Michigan's Robotics Department. Modified for online learning. 

---

## Setup Guide
This guide will walk you through setting up your environment to run Julia files, specifically Jupyter Notebooks written in Julia, using Anaconda Navigator.

**Note:** All required tools are completely free to download and use.

Tools you'll need:
* Anaconda Navigator - a GUI for managing packages, environments, and launching Jupyter Notebooks.
* Julia - the programming language used for ROB 201.
* IJulia - a Julia package that integrates Julia with Jupyter. 

### 1. Install Anaconda
Download and install the latest version of Anaconda (Python 3.x) from the [official website](https://www.anaconda.com/download/success). Once installed, launch Anaconda Navigator from your Start Menu or Applications folder.

### 2. Install Julia
Download and install the Julia language from the [official website](https://julialang.org/downloads/). After installation, note the install location (e.g., `C:\Users\<user>\AppData\Local\Programs\Julia-1.10.0` on Windows).

### 3. Add Julia to Jupyter (IJulia)
Now you need to connect Julia to Jupyter Notebooks.
- Open Julia (from your Start Menu or terminal).
- In the Julia REPL (the command-line interface), run:
  ```julia
  using Pkg
  Pkg.add("IJulia")
  ```
This will install the IJulia package and automatically register Julia as a kernel in Jupyter.

## Running Julia Files in Jupyter
**Option 1: Using Anaconda Navigator**
1. Open Anaconda Navigator.
2. Launch Jupyter Notebook.
3. In your browser, navigate to the location of `JuliaHW0X.ipynb`.
4. Click the notebook to open it.
5. At the top-right, make sure the kernel is set to `Julia` (not Python). You can change it via `Kernel > Change kernel > Julia`.

**Option 2: Using the Terminal**
Alternatively, you can launch Jupyter from the terminal:
```bash
jupyter notebook
```
This will open a local session in your browser. Navigate to the `.ipynb` file and open it.

## Troubleshooting/FAQ

### Launching from Terminal Fails
Make sure that Anaconda and Julia are added to your system's PATH.

### Missing ID Field/Validation Error
You might encounter a message like: “Missing 'id' field in cell” or “Notebook validation failed.”

Don’t worry. Your code and notebook will still run normally.

This warning appears because Jupyter Notebooks now use a newer file format (version 4.5 and above) that includes an id field in each cell. Older versions of Jupyter (especially the classic interface) may not recognize this field and display a warning.

What can you do?
You have two perfectly safe options:
* Ignore the warning. It doesn’t affect functionality.
* Update Jupyter to a more recent version that fully supports the latest notebook format (4.5+).

In either case, your Julia notebook will continue to run as expected.

### Can I Use Something Other Than Anaconda?
Yes, you can use other tools or workflows besides Anaconda to run Julia notebooks, such as installing Jupyter manually via `pip` or using VS Code with the Jupyter extension.

However, make sure the following requirements are met:
* Python is installed on your system (preferably Python 3.7 or later).
* Jupyter Notebook or JupyterLab is installed (via pip install notebook or pip install jupyterlab).
* You’ve added the Julia kernel to Jupyter by installing the IJulia package in Julia (see step 3 above)

Anaconda simply simplifies this setup by bundling Python, Jupyter, and environment management into one package, but it's not required. If you're comfortable managing dependencies on your own, feel free to use your preferred setup!
