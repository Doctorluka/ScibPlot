# ScibPlot
**Codes for:**
- Processing the output of benchmarking of single-cell sequencing integration methods (SCIB).
- Producing nice plot.

**The main functions includes:**
- `scib_calculate_scores` calulates batch correction, bio conservation and overscore following the weighting algorithm proposed by SCIB.
  - The input data should have row names as methods and column names as parameters for SCIB calculations.
- `scib_score2tab` produces and checks necessary parameters for plotting. 
  - You need to manually add some necessary information to the `scib_calculate_scores` output file: 
  - **Methods:** usually the row name of the SCIB result, like "scVI" and "Harmony".
  - **Features:** "HVG" or "FULL", which depends on the corresponding integration method and calculation process.
  - **Scaling:** "scaled" or "unscaled", which depends on the content of the matrix needed for the integration algorithm. For example, scVI requires the input of the raw counts matrix (labeled "unscaled").
  - **Output:** one of "Embedding", "Graph" and "Features", which depends on who the integration method works for.
- `scib_tab2plot` organizes and generates all the elements needed for the drawing based on the objects generated by `scib_score2tab`.
- `scib_NicePlot` produces the final image based on the objects generated by `scib_tab2plot`.
- `scib_OneClick` further encapsulates the three functions (`scib_score2tab`, `scib_tab2plot` and `scib_NicePlot`), so that after running `scib_calculate_scores` and adding the necessary custom data, you can run `scib_OneClick` directly to get the plot.

**Other functions for internal data processing:**
- `add_column_if_missing`
- `check_column_value`

### Code Declaration
1. This project was first released on July 8, 2024, published on GitHub. It is intended for research purposes only and strictly prohibits commercial use.
2. The main code is directly contributed by the HLCA paper. It can be viewed here: [HLCA_reproducibility](https://github.com/LungCellAtlas/HLCA_reproducibility). Importantly, this is not original work, it's just a reorganization and adjustment of existing work.
3. The ATAC part about the source code is removed here, please read the source code if needed.
4. Thanks to the GZDLab team for their support.

Note: When using this code, please comply with related laws and respect the rights of the original author.

### Presentation of the final drawing
![scib_plot](https://github.com/Doctorluka/ScibPlot/blob/main/images/scib_nice_plot.png)

