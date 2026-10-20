# Submission checklist — "A whole-brain *Drosophila* connectome model predicts a duration- and sign-dependent potassium trap after anoxia and a reset-based route to functional recovery" (v1.0, 2026-10-20)

Author: Srihari Muralikrishnan (sole author), Kalaris Labs. ORCID 0009-0004-0220-0361. Correspondence: hello@kalarislabs.com.

## Files

| file | purpose |
|---|---|
| `paper.pdf` | camera-ready, XeLaTeX, Libertinus Serif/Math (fonts embedded), numbered references with DOI links |
| `paper_review.pdf` | review variant: line numbers, 1.5 line spacing |
| `SUPPLEMENT.pdf` | Supplementary Information (tables generated programmatically from the result files) |
| `paper.docx` | Word version built by pandoc from the same source (for journals requiring DOCX) |
| `paper_src.md`, `SUPPLEMENT_src.md` | single source of truth |
| `build_paper.py` → `build_latex.py` | rebuild pipeline: resolve tables → typeset PDF/DOCX |
| `_latex/*.tex` | generated LaTeX sources (upload these + `figures/` for arXiv) |
| `figures/Figure1–8_*.png`, `FigS1–S3_*.png` | figures (raster PNG, ≈2,000–2,400 px wide; see open items) |
| `KEY_CLAIMS.json` | 15 headline claims with numbers, sources and confidence (for communications) |

## Required elements (all present)

- [x] Title page with affiliation, ORCID and correspondence email
- [x] Abstract (~330 words)
- [x] Significance statement (~130 words)
- [x] Keywords
- [x] Main text: Introduction, Results 2.1–2.11, Discussion 3.1–3.6, Limitations
- [x] Methods 5.1–5.13
- [x] Data availability and Code availability (separate sections)
- [x] AI-use disclosure:
  - Claude (Anthropic);
  - Kalaris Labs Research Agent Skills and Myriad;
  - Aryan Kumar Jha credited for AI work.
- [x] CRediT author contributions; Sayan Chowdhury listed as a non-author contributor (Resources: secondary compute provision; institutional support)
- [x] Competing interests (author co-founder of Kalaris Labs; contributor is its Founder & CEO)
- [x] Funding (no external funding; self-funded; secondary compute in kind)
- [x] Ethics statement (computational; no animals or human data)
- [x] Acknowledgments: FlyWire Codex, the FlyWire Consortium and the dataset teams per https://flywire.ai/principles.html; dataset papers cited; Sayan Chowdhury
- [x] Numbered references in order of first citation; DOIs hyperlinked; preprints marked
- [x] Supplementary Information
- [x] Zero `[[PENDING]]` markers

## bioRxiv

- Subject area: **Neuroscience**. Article type: New Results.
- Upload `paper.pdf` as the main file, `SUPPLEMENT.pdf` as a supplementary file, and figures embedded in the PDF.
- Licence suggestion: CC BY 4.0.
- Declare competing interests and funding exactly as in the manuscript.
- bioRxiv screens for clinical claims. The abstract and Section 3.4 state explicitly that the results are model predictions with no clinical implication. Keep this wording.
- AI use: bioRxiv requires disclosure in Methods/Acknowledgments and does not allow AI authorship. Compliant.

## arXiv (q-bio.NC; cross-list physics.bio-ph or nlin.AO optional)

- Upload the LaTeX source: `_latex/paper.tex` plus `figures/`.
  - Change `\graphicspath{{../}}` to `{{./}}`, or place `figures/` beside the `.tex` file.
  - The build uses XeLaTeX with `fontspec` and Libertinus OTF loaded via `Path=`. arXiv supports XeLaTeX only through a `00README` (`compiler: xelatex`).
  - Either bundle the Libertinus `.otf` files and drop the `Path=` option, or switch to `\usepackage{libertinus}` (TeX Live) on arXiv.
  - Alternatively, upload the PDF only. arXiv accepts PDF-only submissions when they were not produced from TeX, but because this one was produced from TeX, arXiv will ask for source. Prefer source.
- References are inline `thebibliography`, so no `.bbl`/`.bib` is needed.
- Abstract field: plain text, ≤1,920 characters. Use the abstract without the bold run-in headings.
- Comments field: "27 pages, 8 figures, supplementary information; model predictions only".

## eLife (Reviewed Preprint route)

- Submit the bioRxiv preprint DOI. eLife reviews the preprint as posted.
- eLife wants an eLife digest (plain language, optional at submission). The significance statement can be adapted.
- Required: "Data availability", "Competing interests", "Funding" and "Ethics" statements (present). eLife requires a public code/data repository by acceptance. The GitHub repository is marked planned and the Zenodo DOI is not yet minted.
- Figures: eLife prefers figures uploaded separately (TIFF/EPS/PDF, ≥300 dpi at final size). Figure supplements can carry S1–S3.
- Reporting: complete the eLife transparent reporting form (computational study: seeds, CIs, preregistration date 2026-10-16, deviations in Supplement S9).

## PLOS Computational Biology

- Article type: Research Article. Required sections: Author summary (use the significance statement, 150–200 words), Abstract (≤300 words; **the current abstract is ~330 words — trim ~30 words**), Introduction, Results, Discussion, Methods (Methods may follow Discussion), Supporting information captions.
- Vancouver-style numbered references in order of citation: already satisfied.
- Figures: TIFF or EPS, 300–600 dpi, width ≤ 19.05 cm. Convert the PNG composites with PACE (https://pacev2.apexcovantage.com). Supporting tables and figures go as S1_Table etc.
- Data Availability Statement: PLOS requires data to be publicly available at acceptance. Mint the Zenodo DOI and make the repository public before submission.
- Financial disclosure: "The author received no specific funding for this work." Competing interests: as in the manuscript.
- Line numbers and double spacing: use `paper_review.pdf` (1.5 spacing). Set `\doublespacing` if strictly required.

## Open items before external submission

1. **Figures are raster PNG composites, not vector PDF.** Panels were produced by the track analysis scripts as PNG and stacked. Vector re-export means re-running each track's plotting script with `savefig('*.pdf')`; this was not done in v1.0.
2. Make https://github.com/KalarisLabs/fly-brain-resuscitation public and mint the Zenodo DOI; update the Data and Code availability sections.
3. PLOS CB: shorten the abstract to ≤300 words and provide a 150–200-word author summary.
4. arXiv: adapt font loading (see above).
5. Optional before journal submission:
   - analyse the 146 BANC replication files;
   - run sign variants B/C and the FAFB variant-A baseline, reset and 180-s arms;
   - complete the cross-connectome recoverability maps.
