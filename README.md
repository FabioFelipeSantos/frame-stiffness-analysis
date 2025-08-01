# 2D Frame Analysis using the Stiffness Method

This is a MATLAB project for the structural analysis of 2D (plane) frames using the Direct Stiffness Method. The code is capable of calculating nodal displacements, support reactions, and internal forces, and it visually presents the results through diagrams of the deformed structure and internal forces.

## Features

-   **2D Frame Analysis**: Analyzes structures composed of beams and columns in a 2D plane.
-   **Direct Stiffness Method**: Implements the matrix-based stiffness method to solve for the behavior of the structure.
-   **Multiple Load Types**: The code can handle various types of loads, including:
    -   Nodal forces and moments applied directly to the joints.
    -   Concentrated forces and moments applied at any point along an element.
    -   Uniformly or linearly distributed loads applied over a portion or the entire length of an element.
-   **Comprehensive Results**: The analysis provides a complete set of results:
    -   Nodal displacements and rotations.
    -   Support reactions.
    -   Internal forces (normal force, shear force, and bending moment) for each element.
-   **Graphical Outputs**: Generates plots to visualize the results:
    -   Deformed shape of the structure.
    -   Normal Force Diagram.
    -   Shear Force Diagram.
    -   Bending Moment Diagram.
-   **Detailed Report**: Creates a `Dados Estrutura.txt` file with detailed step-by-step results, including local and global stiffness matrices, rotation matrices, force vectors, and final displacements for debugging and verification.

## How It Works

The main script, `Portico_Plano.m`, orchestrates the entire analysis process:

1.  **Input Data**: The script reads the structural data from a specified `.xlsx` Excel file. This file must contain three sheets:
    -   `Nós` (Nodes): Defines the coordinates and support conditions for each node.
    -   `Elementos (Vigas)` (Elements/Beams): Defines the elements by their connecting nodes and material/geometric properties (E, A, I).
    -   `Forças` (Forces): Details all the loads applied to the structure.
2.  **Stiffness Matrix Assembly**: For each element, the script calculates the local stiffness matrix (`matRigVigaPortico.m`), rotates it to the global coordinate system (`matrizRotacao.m`), and assembles the global stiffness matrix `K` for the entire structure.
3.  **Force Vector Assembly**:
    -   It calculates the fixed-end forces and moments for any loads applied along the elements (`momentosDeExtremoFixoCargaDistribuida.m`, `momentosDeExtremoFixoCargaConcentrada.m`). These are assembled into a global fixed-end force vector `Qf`.
    -   Nodal loads from the input file are assembled into the global force vector `f`.
4.  **Applying Boundary Conditions**: The script identifies and partitions the degrees of freedom based on the support conditions (restraints) defined in the `Nós` sheet.
5.  **Solving the System**: It solves the fundamental system of equations `K * u = f - Qf` to find the unknown nodal displacements and rotations `u`.
6.  **Post-processing**:
    -   Calculates the support reactions using the computed displacements.
    -   Determines the local internal forces for each element.
    -   Generates plots for the deformed shape and internal force diagrams (`EsforcoNormal.m`, `EsforcoCortante.m`, `EsforcoMomento.m`).

## How to Run

1.  **Prerequisites**: Make sure you have MATLAB installed.
2.  **Setup**:
    -   Place all the `.m` files in the same directory.
    -   Prepare your input data in an Excel file (e.g., `MyStructure.xlsx`). You can use the provided `.csv` files as a template for creating the sheets.
3.  **Input File Structure**:
    -   **Sheet 1: `Nós`**
        -   `Coord X`, `Coord Y`: Nodal coordinates.
        -   `Rest X`, `Rest Y`, `Rest MZ`: Boundary conditions. Use `1` for a restrained (fixed) degree of freedom and `0` for a free one.
    -   **Sheet 2: `Elementos (Vigas)`**
        -   `No i`, `No f`: The start and end nodes for the element.
        -   `E`: Modulus of Elasticity.
        -   `A`: Cross-sectional Area.
        -   `I`: Moment of Inertia.
    -   **Sheet 3: `Forças`**
        -   `Tipo`: Load type (`1` for concentrated, `2` for distributed, `3` for concentrated moment).
        -   `Viga`: The element number on which the load is applied (leave blank for nodal loads).
        -   `Coor X/Y Inicial/Final`: Coordinates for the start and end points of the load.
        -   `Valor Inicial/Final`: Magnitude of the load. For uniform loads, the initial and final values are the same.
        -   `Direção`: The angle of the load application in degrees (e.g., `-90` for a downward vertical load).
4.  **Execution**:
    -   Open the `Portico_Plano.m` file in MATLAB.
    -   On line 12, change the `nome` variable to your Excel file's name:
        ```matlab
        nome = 'MyStructure.xlsx';
        ```
    -   Run the script (press F5 or click the "Run" button).

## Expected Results

After running the script, you will get:

-   **Console Output**: The support reactions will be printed in the MATLAB command window.
-   **Detailed Text File**: A file named `Dados Estrutura.txt` will be generated in the same directory, containing all the intermediate matrices and vectors used in the calculation.
-   **Plot Windows**: Four figures will be generated, showing:
    1.  The deformed shape of the structure (in red) superimposed on the original structure (in blue dashed lines).
    2.  The Normal Force Diagram.
    3.  The Shear Force Diagram.
    4.  The Bending Moment Diagram.

## Project Files

-   `Portico_Plano.m`: The main executable script.
-   `matRigVigaPortico.m`: Calculates the 6x6 local stiffness matrix for a frame element.
-   `matrizRotacao.m`: Creates the transformation matrix to rotate element matrices from local to global coordinates.
-   `momentosDeExtremoFixoCarga*.m`: Functions to calculate fixed-end moments for concentrated and distributed loads.
-   `Esforco*.m`: Functions to calculate the internal normal, shear, and moment values along an element for plotting.
-   `FuncoesDeForma*.m`: Shape functions for beam and bar elements used for plotting the deformed shape.
-   `coisas_simbolicas.m`: A script with symbolic math derivations of the shape functions and fixed-end moment formulas.
-   `*.xlsx`: Excel files containing the input data for different example structures.
-   `*.csv`: CSV versions of the sheets in the Excel files.
