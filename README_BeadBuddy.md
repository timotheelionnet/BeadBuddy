# BeadBuddy

**BeadBuddy** is a flexible, user-friendly MATLAB tool for correcting chromatic errors in multi-channel single molecule imaging datasets. Designed to streamline the correction process through a guided graphical user interface (GUI), BeadBuddy integrates seamlessly with typical imaging workflows and data formats.

---

## 🧰 Features

- Automated localization and analysis of bead images
- Visualization of chromatic aberrations per channel
- Nearest neighbor analysis to quantify correction effectiveness
- Saves diagnostic plots and corrected data outputs
- Intuitive, step-by-step GUI

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
   git clone https://github.com/yourusername/BeadBuddy.git
   ```
2. Add the folder to your MATLAB path:
   ```matlab
   addpath(genpath('path/to/BeadBuddy'))
   ```

---

## 🚀 Getting Started

You can practice with the demo dataset included in the repository. A full [video tutorial](#) is available to guide you through the process.

### Quickstart Workflow
1. **Thermalize the microscope** (~2.5–3 hours).
2. **Acquire bead dataset**:
   - 3×3 grid of Z-stacks
   - 2.5% overlap
   - All imaging channels used in your experiment
3. **Acquire and localize single molecule data**.
   - Your localization table must include these columns:
     - `channel`
     - `x`
     - `y`
     - `z` *(omit for 2D datasets)*
     - `FOV`
4. **Launch BeadBuddy GUI**:
   - The GUI will guide you through channel assignment and correction.

---

## 📊 Output

BeadBuddy generates the following:
- Overlay plots showing chromatic aberrations between channels
- Diagnostic figures saved to the project folder
- Tables of corrected localization coordinates
- Nearest-neighbor histograms comparing pre- and post-correction

---

## 🖼️ Example GUI

*(Insert GUI screenshot here)*

![BeadBuddy GUI](path/to/example_gui.png)

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
