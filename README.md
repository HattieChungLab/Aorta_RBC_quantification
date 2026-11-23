# Red-cell quantification in haemorrhagic lesions of the aortic media

Results for all 25 variants are in `RESULTS.md`; all settings are in `config/base.yaml`.

## Target

The number and area of red blood cells inside outlined haemorrhagic lesions of the aortic
media, per lesion and per slide.

- **Data:** 22 H&E whole-slide images (40×, 0.2524 µm/px). Each slide is annotated in QuPath with the
  media, the false lumen and the haemorrhagic lesions. Slides are referred to by code (S01–S22).
- **Counting:** a red cell is counted in a lesion if its centre lies inside the lesion outline. Red-cell
  area is the summed area of the counted cells' outlines.

## Methods

No model was trained or tuned on these slides.

- **HistoPLUS** (CellViT, 40× model, 15 cell types including red blood cell) labels every pixel with a
  cell type and also produces its own cell outlines.
- **Six general cell segmenters** outline cells but do not classify them: Cellpose-SAM v2,
  Cellpose-SAM, Cellpose-DINO (ViT-L and ViT-B backbones), InstanSeg and micro-SAM (`vit_b_lm`).
- **Two image inputs** for the segmenters: the colour H&E image, and a red-excess image,
  R − (G+B)/2, in which red cells are bright.
- **+ HP:** a segmenter outline is kept as a red cell if the majority HistoPLUS pixel label inside it is
  "red blood cell".
- **+ HP + fill:** as + HP, plus HistoPLUS's own red-cell outlines that touch no kept outline.

This gives **25 variants**: 6 segmenters × 2 inputs × (+ HP, + HP + fill), plus HistoPLUS alone.

## Validation

Two sets of annotations were used:

- **Traced cells:** 126 red cells outlined individually by the pathologist, 43 on S11 (deep-red stain)
  and 83 on S20 (pale stain). They are a sample of the cells present, so only misses can be measured.
  Each traced cell is scored in a 113 µm square around it.
  - **Found:** a detected red cell overlaps the traced outline.
  - **Counted once:** exactly one counted cell centre lies inside the traced outline (missed = none,
    split = two or more).
- **Fully counted fields:** three 100 µm squares of dense clot (on S11, S17 and S02) in which every red
  cell was marked (351 cells). Detected cells touching the square are included.
  - **TP%:** the share of marked cells that fall inside a detected red-cell outline.
  - **FP%:** the share of detected red cells with no mark inside ("–" = no detections).

**Table 1. Traced cells** (top performers; full table in `RESULTS.md`). Best value per column in bold.

| Method | Input | S11 (deep red, 43): found | S11: counted once (missed / split) | S20 (pale, 83): found | S20: counted once (missed / split) |
|---|---|---|---|---|---|
| HistoPLUS alone (reference) | colour | 93% | 81% (16% / 2%) | 57% | 49% (49% / 1%) |
| Cellpose-SAM v2 + HP + fill | red-excess | **98%** | **95%** (5% / 0%) | 83% | 72% (28% / 0%) |
| Cellpose-DINO (ViT-L) + HP + fill | red-excess | 95% | **95%** (5% / 0%) | 87% | 75% (25% / 0%) |
| Cellpose-DINO (ViT-B) + HP + fill | red-excess | 95% | 88% (12% / 0%) | 88% | 76% (23% / 1%) |
| micro-SAM + HP + fill | red-excess | **98%** | 93% (7% / 0%) | **93%** | **78%** (20% / 1%) |

**Table 2. Fully counted fields: TP% / FP%** (top performers; full table in `RESULTS.md`). Best value
per column in bold (highest TP%, lowest FP%).

| Method | Input | S11 (193 cells) | S17 (64 cells) | S02 (94 cells) | All three (351 cells) |
|---|---|---|---|---|---|
| HistoPLUS alone (reference) | colour | 81% / 15% | 52% / 24% | 85% / 31% | 77% / 22% |
| Cellpose-SAM v2 + HP + fill | red-excess | 84% / 14% | **61%** / **18%** | **88%** / **26%** | 81% / **18%** |
| Cellpose-DINO (ViT-L) + HP + fill | red-excess | 89% / 14% | 59% / 20% | 85% / 29% | 83% / 19% |
| Cellpose-DINO (ViT-B) + HP + fill | red-excess | 84% / **13%** | 56% / **18%** | 80% / 29% | 78% / **18%** |
| micro-SAM + HP + fill | red-excess | **95%** / 16% | 58% / 32% | 84% / 30% | **85%** / 22% |

Red-cell quantification uses **micro-SAM + HP + fill** on the red-excess image.

## Outputs

- `outputs/all_methods/<slide>/lesions.csv`: red-cell count and area for every lesion and every variant.
- `outputs/all_methods_report/`: the same per slide (`counts_per_slide.csv`, `area_per_slide.csv`),
  lesion and medial areas per slide (`slide_info.csv`), and a summary (`report.md`).
- `outputs/validation/`: the validation tables (`traced_cells`, `complete_fields`).
- `outputs/figures/`: pictures of every variant on the fully counted fields and on individual lesions.
