# Full validation results (all 25 variants)

Companion to `README.md`, which shows the top performers. Metrics are defined in the README
(Validation).

**Table 1. Traced cells: found and counted once (missed / split).**

| Method | Input | S11 (deep red, 43): found | S11: counted once (missed / split) | S20 (pale, 83): found | S20: counted once (missed / split) |
|---|---|---|---|---|---|
| HistoPLUS alone | colour | 93% | 81% (16% / 2%) | 57% | 49% (49% / 1%) |
| Cellpose-SAM v2 + HP | red-excess | 84% | 81% (19% / 0%) | 63% | 57% (43% / 0%) |
| Cellpose-SAM v2 + HP + fill | red-excess | 98% | 95% (5% / 0%) | 83% | 72% (28% / 0%) |
| Cellpose-SAM v2 + HP | colour | 0% | 0% (100% / 0%) | 0% | 0% (100% / 0%) |
| Cellpose-SAM v2 + HP + fill | colour | 93% | 81% (16% / 2%) | 57% | 49% (49% / 1%) |
| Cellpose-SAM + HP | red-excess | 67% | 65% (35% / 0%) | 55% | 51% (49% / 0%) |
| Cellpose-SAM + HP + fill | red-excess | 93% | 91% (9% / 0%) | 78% | 69% (31% / 0%) |
| Cellpose-SAM + HP | colour | 0% | 0% (100% / 0%) | 0% | 0% (100% / 0%) |
| Cellpose-SAM + HP + fill | colour | 93% | 81% (16% / 2%) | 57% | 49% (49% / 1%) |
| Cellpose-DINO (ViT-L) + HP | red-excess | 93% | 93% (7% / 0%) | 58% | 49% (51% / 0%) |
| Cellpose-DINO (ViT-L) + HP + fill | red-excess | 95% | 95% (5% / 0%) | 87% | 75% (25% / 0%) |
| Cellpose-DINO (ViT-L) + HP | colour | 0% | 0% (100% / 0%) | 0% | 0% (100% / 0%) |
| Cellpose-DINO (ViT-L) + HP + fill | colour | 93% | 81% (16% / 2%) | 57% | 49% (49% / 1%) |
| Cellpose-DINO (ViT-B) + HP | red-excess | 91% | 81% (19% / 0%) | 66% | 59% (40% / 1%) |
| Cellpose-DINO (ViT-B) + HP + fill | red-excess | 95% | 88% (12% / 0%) | 88% | 76% (23% / 1%) |
| Cellpose-DINO (ViT-B) + HP | colour | 0% | 0% (100% / 0%) | 0% | 0% (100% / 0%) |
| Cellpose-DINO (ViT-B) + HP + fill | colour | 93% | 81% (16% / 2%) | 57% | 49% (49% / 1%) |
| InstanSeg + HP | red-excess | 56% | 30% (70% / 0%) | 43% | 29% (71% / 0%) |
| InstanSeg + HP + fill | red-excess | 95% | 72% (28% / 0%) | 72% | 53% (46% / 1%) |
| InstanSeg + HP | colour | 0% | 0% (100% / 0%) | 0% | 0% (100% / 0%) |
| InstanSeg + HP + fill | colour | 93% | 81% (16% / 2%) | 57% | 49% (49% / 1%) |
| micro-SAM + HP | red-excess | 95% | 91% (9% / 0%) | 92% | 77% (22% / 1%) |
| micro-SAM + HP + fill | red-excess | 98% | 93% (7% / 0%) | 93% | 78% (20% / 1%) |
| micro-SAM + HP | colour | 91% | 91% (9% / 0%) | 28% | 24% (75% / 1%) |
| micro-SAM + HP + fill | colour | 95% | 95% (5% / 0%) | 70% | 59% (39% / 2%) |

**Table 2. Fully counted fields: TP% / FP%.**

| Method | Input | S11 (193 cells) | S17 (64 cells) | S02 (94 cells) | All three (351 cells) |
|---|---|---|---|---|---|
| HistoPLUS alone | colour | 81% / 15% | 52% / 24% | 85% / 31% | 77% / 22% |
| Cellpose-SAM v2 + HP | red-excess | 67% / 13% | 28% / 16% | 54% / 13% | 56% / 14% |
| Cellpose-SAM v2 + HP + fill | red-excess | 84% / 14% | 61% / 18% | 88% / 26% | 81% / 18% |
| Cellpose-SAM v2 + HP | colour | 0% / 100% | 0% / – | 0% / – | 0% / 100% |
| Cellpose-SAM v2 + HP + fill | colour | 79% / 14% | 52% / 22% | 83% / 31% | 75% / 21% |
| Cellpose-SAM + HP | red-excess | 44% / 13% | 17% / 17% | 20% / 24% | 33% / 15% |
| Cellpose-SAM + HP + fill | red-excess | 82% / 13% | 56% / 21% | 81% / 32% | 77% / 20% |
| Cellpose-SAM + HP | colour | 0% / 100% | 0% / – | 0% / – | 0% / 100% |
| Cellpose-SAM + HP + fill | colour | 79% / 14% | 52% / 22% | 83% / 31% | 75% / 21% |
| Cellpose-DINO (ViT-L) + HP | red-excess | 87% / 13% | 30% / 23% | 45% / 21% | 65% / 15% |
| Cellpose-DINO (ViT-L) + HP + fill | red-excess | 89% / 14% | 59% / 20% | 85% / 29% | 83% / 19% |
| Cellpose-DINO (ViT-L) + HP | colour | 0% / 100% | 0% / – | 0% / – | 0% / 100% |
| Cellpose-DINO (ViT-L) + HP + fill | colour | 79% / 14% | 52% / 22% | 83% / 31% | 75% / 21% |
| Cellpose-DINO (ViT-B) + HP | red-excess | 80% / 12% | 14% / 0% | 51% / 17% | 60% / 12% |
| Cellpose-DINO (ViT-B) + HP + fill | red-excess | 84% / 13% | 56% / 18% | 80% / 29% | 78% / 18% |
| Cellpose-DINO (ViT-B) + HP | colour | 0% / 100% | 0% / – | 0% / – | 0% / 100% |
| Cellpose-DINO (ViT-B) + HP + fill | colour | 79% / 14% | 52% / 22% | 83% / 31% | 75% / 21% |
| InstanSeg + HP | red-excess | 36% / 7% | 34% / 42% | 36% / 7% | 36% / 11% |
| InstanSeg + HP + fill | red-excess | 78% / 13% | 61% / 23% | 79% / 31% | 75% / 20% |
| InstanSeg + HP | colour | 0% / – | 0% / – | 0% / – | 0% / – |
| InstanSeg + HP + fill | colour | 79% / 14% | 52% / 22% | 83% / 31% | 75% / 21% |
| micro-SAM + HP | red-excess | 94% / 16% | 25% / 37% | 80% / 26% | 78% / 20% |
| micro-SAM + HP + fill | red-excess | 95% / 16% | 58% / 32% | 84% / 30% | 85% / 22% |
| micro-SAM + HP | colour | 90% / 20% | 0% / 100% | 52% / 41% | 63% / 25% |
| micro-SAM + HP + fill | colour | 93% / 20% | 50% / 25% | 79% / 39% | 81% / 26% |
