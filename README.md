# is-2dspu

Implementation of the Iterative Search algorithm proposed in the paper [*An open space based heuristic for the 2D strip packing problem with unloading constraints*](https://doi.org/10.1016/j.apm.2019.01.022) to solve the Two-Dimensional Strip Packing Problem with Unloading Constraints.

## Pre-requisites 🛠️

You must have C++17 and GNU g++ 11.4.0

## How to run ⚙️

```bash
    ./run.sh -f instance_path -s seed
```

Available flags:

* **-f <file>**: Specifies the path of the input instance. 📂 [required]
* **-s <string>**: Sets the seed for generating random numbers. 🌱 [required]
* **-t <number>**: Sets the time limit for execution (defaults to 60s). ⏳
* **-d**: Activates debug mode, the output is a file you can use as input in the `debug.py` script. 🐞
* **-h**: Shows this help menu. 📖

## How to debug 🐞

An auxiliary Python tool was developed to help debug the strip packing layout solution. It plots the arranged items on the strip and their maximum width to visualize the generated solutions.

**Visualization Details:**
* Numbers on the items indicate their position in the item vector.
* Items from the same client share the same color.
* The color palette maps client classes: **Higher classes** are closer to **red**, and **lower classes** are closer to **violet**.

```bash
    ./debug.py solution_file_path
```