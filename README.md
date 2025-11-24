
# BeadBuddy

**BeadBuddy** is a flexible, user-friendly MATLAB tool for correcting chromatic errors in multi-channel single molecule imaging datasets. Users simply image a fluorescent calibration bead slide at the beginning of their imaging session and then BeadBuddy handles the rest. 

BeadBuddy analyzes the bead images to model arbitrary chromatic errors in three dimensions as a function of position in the imaging plane. It then applies the corrections to the user's dataset resulting in a more precise, more reproducible downstream analysis!

**Watch the [BeadBuddy video tutorial](https://www.youtube.com/watch?v=72aow-Y2Re4)**
**Test it out on our  [demo dataset](https://doi.org/10.5281/zenodo.17605122)**

---
## BeadBuddy Overview & GUI

*(Insert GUI screenshot here)*

![BeadBuddy GUI](bb_gui.png)
---

## 🧰 Features

- Intuitive, step-by-step graphical user interface
- Automated localization and analysis of fluorescent calibration bead images
- Visualization of chromatic errors present in your imaging system and dataset
- Handles 2D or 3D single-molecule localization datasets


---

## 📦 Installation

### Requirements
- MATLAB (R2020b or newer recommended)
- The following MATLAB Toolboxes:
  - Optimization Toolbox
  - Image Processing Toolbox
  - Statistics and Machine Learning Toolbox
  - Curve Fitting Toolbox

### Setup
1. Clone or download this repository:
   ```bash
   git clone https://github.com/timotheelionnet/BeadBuddy
   ```
2. In MATLAB, use the folder browser to navigate to the BeadBuddy-main folder. This is now your working directory.

---

## 🚀 Getting Started

You can practice with the demo dataset included in the repository. A full [video tutorial](#) is available to guide you through the process.

### Workflow
1. **Thermalize the microscope** (~2.5–3 hours).
2. **Acquire calibration bead dataset**:
   - 3×3 grid of bead Z-stacks (must be imaged in 3D regardless) ~2.5% overlap
   - Use all channels relevant to the dataset you will go on to correct
3. **Acquire and localize your single molecule data**.
   - Your dataset must be .csv or .xlsx table and include either 5 or 4 columns:
     - `channel`
     - `x`
     - `y`
     - `z` *(omit this column for 2D datasets)*
     - `FOV`
4. **Launch BeadBuddy GUI**:
   - In MATLAB, after navigating to BeadBuddy-main, go to the Command Window and execute:
      - `BeadBuddy_MAIN`
   - The GUI will guide you through the rest.

---

## 📊 Output

BeadBuddy generates the following output folders in your project folder:
- results
   - "Your data_BEADCRXN" --> corrected user single molecule dataset in same units as input (pix or nm)
   - NN_analysis_CDF_plots
      - Nearest neighbor analysis of the user's dataset before and after BeadBuddy correction plotted as eCDFS.
   - registration_viz
      - Plots visualizing the actual correction that was applied to the user's single molecule dataset
- bead analysis
   - displacement_tiffs
      - tiff files with same dimensionality as bead images
      - each pixel stores the displacement field function  evaluated at the center of that pixel (units = nm)
      - file names indicate relevant spatial dimension and channel
   - heatmaps
      - each heatmap visualizes the 3D chromatic error profile modeled for the specified non-reference channel 



---

## 📜 License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## 📚 Citation

If you use BeadBuddy in your research, please cite the associated publication:

```
Author(s). "Title of Paper." *Journal Name*, Year. DOI
```

*(Insert full citation here once available.)*

---

## 💬 Questions & Contact

Feel free to [open an issue](https://github.com/yourusername/BeadBuddy/issues) or reach out directly for help.

We welcome feedback, suggestions, and contributions!
