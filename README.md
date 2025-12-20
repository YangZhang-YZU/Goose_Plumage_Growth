## 1. Disentangling Confounded Signals: Paired Conditional GWAS

To rigorously distinguish genetic effects driven by genuine growth potential from those confounded by plumage development, we implemented a **Paired Conditional GWAS strategy**.

**Methodology:**
We performed two parallel GWAS runs for each growth trait:
1.  **Unconditional Model:** Standard LMM ($y = W\alpha + x\beta + u + \varepsilon$).
2.  **Conditional Model:** LMM including plumage color as a fixed covariate ($c$).

**The Delta ($\Delta$) Statistic:**
We quantified the shift in association strength for each SNP using the $\Delta$ statistic:

$$\Delta = -\log_{10}(P_{cond}) - [-\log_{10}(P_{uncond})]$$

**Interpretation of Results:**
The `src/02_conditional_gwas/` pipeline categorizes variants into three distinct classes:
* **$\Delta > \epsilon$ (Growth-Specific Signals):** Variants where significance increased or remained high after conditioning, indicating robust effects on body weight independent of plumage.
* **$\Delta < -\epsilon$ (Plumage-Confounded Signals):** Variants where significance dropped drastically, suggesting the original signal was driven by the correlation with feather color.
* **Stable Signals:** Variants showing negligible change.

*Resulting files are located in `output/paired_gwas_delta_results.txt`.*

## 2. Multi-Trait Integration: The Composite E-score

To synthesize association signals across static weight measurements ($W_0-W_{10}$) and dynamic growth model parameters (Gompertz/Logistic $A, k, t_0$), we developed the **E-score**. This composite statistic highlights pleiotropic variants that systematically affect the entire growth curve.

**Algorithm:**
The E-score is a weighted summation of association signals, where weights are derived from the genomic inflation factor ($\chi^2$) of each trait. This ensures that traits with higher heritability and signal-to-noise ratios contribute more to the final score.

**Mathematical Formulation:**
For a given SNP $i$ and a set of traits $\mathcal{G}$:

$$E(i) = -\sum_{j \in \mathcal{G}} \left( \frac{\max(0, \bar{\chi}_j^2 - 1)}{\sum_{k} \max(0, \bar{\chi}_k^2 - 1)} \cdot \log_{10} P_{j}(i) \right)$$

**Key Features:**
* **Variance-Weighted:** Weights are automatically adjusted based on `Mean Chi2 - 1`, penalizing traits with no genomic inflation.
* **Vectorized Computation:** The implementation in `src/03_mtag_escore/calc_E_score.py` utilizes vectorization for high-performance processing of genome-wide data.

*The final E-score summary statistics are provided in `output/escore_results.txt`, ready for Manhattan plot visualization.*
