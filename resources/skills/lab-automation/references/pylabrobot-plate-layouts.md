# PyLabRobot `LiquidHandler` & Plate-Layout Templates

Load this reference when building hardware-agnostic **PyLabRobot** protocols or designing **96/384-well experimental plate layouts** (serial dilution, drug synergy checkerboard, qPCR, or immunoassay).

## 1. PyLabRobot `LiquidHandler` Pattern (Dry-Run First)

PyLabRobot separates protocol logic from hardware backends (`ChatterboxBackend` for zero-hardware dry-run validation, `OT2Backend` for Opentrons OT-2, `STARBackend` for Hamilton STAR):

```python
import asyncio
from pylabrobot.liquid_handling import LiquidHandler
from pylabrobot.liquid_handling.backends.simulation.chatterbox_backend import ChatterboxBackend
from pylabrobot.resources import Deck, TipRack, Plate, set_tip_tracking, set_volume_tracking
from pylabrobot.resources.corning_costar import Cos_96_DW_1mL, Cos_96_Vb
from pylabrobot.resources.opentrons import OT_2_Deck, opentrons_96_tiprack_300ul


async def run_serial_dilution() -> None:
    set_tip_tracking(True)
    set_volume_tracking(True)  # Tracks per-well remaining volume automatically

    deck = OT_2_Deck()
    lh = LiquidHandler(backend=ChatterboxBackend(), deck=deck)
    await lh.setup()

    tips = opentrons_96_tiprack_300ul(name="tips_300")
    reagent_trough = Cos_96_DW_1mL(name="diluent_dw")
    assay_plate = Cos_96_Vb(name="assay_96_v")

    deck.assign_child_resource(tips, slot=1)
    deck.assign_child_resource(reagent_trough, slot=2)
    deck.assign_child_resource(assay_plate, slot=3)

    # Pre-populate tracked volumes (must include >= 15 uL V-bottom dead volume)
    reagent_trough["A1"].tracker.set_volume(800)  # 600 uL needed + 200 uL overage

    # Step 1: Dispense 100 uL diluent to cols 2..7 (>= 10% of 300 uL tip = 30 uL)
    await lh.pick_up_tips(tips["A1"])
    for col in range(2, 8):
        await lh.aspirate(reagent_trough["A1"], vols=[100])
        await lh.dispense(assay_plate[f"A{col}"], vols=[100])
    await lh.drop_tips(tips["A1"])

    # Step 2: Serial 1:3 transfer (50 uL into 100 uL) with fresh tip per step or high->low mix
    for col in range(1, 7):
        await lh.pick_up_tips(tips[f"B{col}"])
        await lh.aspirate(assay_plate[f"A{col}"], vols=[50], mix=[{"volume": 75, "repetitions": 3}])
        await lh.dispense(assay_plate[f"A{col+1}"], vols=[50])
        await lh.drop_tips(tips[f"B{col}"])

    await lh.stop()


if __name__ == "__main__":
    asyncio.run(run_serial_dilution())
```

> Enable `set_tip_tracking(True)` and `set_volume_tracking(True)` during simulation so PyLabRobot catches empty wells, over-filled wells, and tip collisions before hardware execution.

## 2. Plate Geometry & Dead-Volume Reference

| Labware Format | Working Well Volume | Mandatory Dead-Volume Floor | Notes |
| --- | --- | --- | --- |
| **96-well V-bottom** | `50–200 uL` | `10–15 uL` per well | Best recovery for small-volume compounds/primers |
| **96-well flat / U-bottom** | `100–300 uL` | `15–25 uL` per well | Higher residual meniscus; never aspirate last `15 uL` |
| **96-well PCR (skirted)** | `10–100 uL` | `5–10 uL` per well | Use `0.5–1.0 mm` bottom clearance + spin-down |
| **384-well low-volume** | `10–40 uL` | `3–5 uL` per well | Strictly forbid tips `> 50 uL` nominal max |
| **12-channel / 1-well trough** | `5–15 mL` / `100–250 mL` | `+10%–15%` total overage | Account for angled dead space at trough bottom |

## 3. Production Plate-Layout Templates

### Template A: Anti-Edge-Bias 96-Well Dose-Response (10-Point + Stratified Controls)

Avoid edge evaporation bias by filling **perimeter wells (`Row A`, `Row H`, `Col 1`, `Col 12`)** with `200 uL` sterile PBS/water as a thermal/evaporation moat during multi-day cell assays, leaving **60 active interior wells (`B2–G11`)**:

- **Rows B–G, Cols 2–11 (`60 wells`)**:
  - **2 drugs/conditions x 10-point 1:3 serial dilution in duplicate** (`40 wells`), OR **1 drug x 10-point in quadruplicate** (`40 wells`).
  - **Vehicle (DMSO-matched) & Positive Controls (`12–20 wells`)**: Interleaved via **Latin-square / diagonal column offsets** across interior columns (`Cols 2, 5, 8, 11`), NEVER clustered in a single edge column.
  - **Solvent Normalization**: Match final DMSO concentration (`<= 0.5% v/v`) across every active well including vehicle controls.
- **Downstream Phi Handoff (`pgx_hts_qc`)**: Export a companion plate-map CSV with columns `plate_id, well, row, col, role, compound, concentration_uM, batch` (`role` in `{sample, pos_ctrl, neg_ctrl, blank}`), allowing `pgx_hts_qc` to compute $Z'$-factor, SSMD, row/column ANOVA, edge-effect diagnostics, and Tukey `B-score`.

### Template B: 6x10 / 8x8 Checkerboard Drug Synergy Matrix

- **96-well interior (`B2–G11`, 60 wells = 6 x 10)**:
  - **Drug A** diluted vertically (`Rows B -> G`: `C_max -> C_max/16 -> 0` in Row G = Drug B monotherapy).
  - **Drug B** diluted horizontally (`Cols 2 -> 10`: `C_max -> C_max/128 -> 0` in Col 10 = Drug A monotherapy).
  - **Col 11 (`B11–G11`)**: `3` DMSO vehicle wells + `3` positive-kill control wells (interleaved `neg, pos, neg, pos, neg, pos`).
- **384-well scaling**: Map 4 independent `6x10` quadrants (e.g., 2 biological replicates x 2 cell lines) into interior wells `B2–O23` with perimeter moat `Row A, Row P, Col 1, Col 24`.
- **Transfer rule**: Pre-dispense normalized DMSO backfill FIRST so total DMSO is identical in every well before adding Drug A and Drug B serial stocks (`new_tip="always"` per concentration step).

### Template C: MIQE qPCR & Immunoassay Calibration Plate Map

1. **RT-qPCR (`wa_qpcr_quant` compatible)**:
   - **Master-Mix First**: `distribute` gene-specific master mix (primer + probe/SYBR) into empty plate wells.
   - **Controls**: Include **NTC** (no-template control, `>= 2` wells per primer pair), **`-RT` control**, and **5-point 1:10 pooled cDNA standard curve** (efficiency `90%–110%`, `R^2 >= 0.99`) inside interior wells.
   - **Biological vs Technical Replicates**: Place technical triplicates adjacent or block-randomized per biological sample (`biological_unit` column preserved in plate map) and spread biological groups evenly across columns/plates to prevent run confounding.
   - **Downstream Export**: Record `well, sample_id, group, biological_unit, target_gene, dilution_factor, is_ntc` so `wa_qpcr_quant` validates NTC contamination, technical replicate `SD <= 0.5`, geNorm/Andersen stability, and biological-replicate stats.

2. **ELISA / Luminex / BCA (`wa_immunoassay_curve` compatible)**:
   - **Standard Curve**: `7–8` serial dilutions (typically 1:2 or 1:3) + **Blank (`0`)** in duplicate/triplicate, split across two non-adjacent plate zones (e.g., `B2–B9` and `G4–G11`) to detect intra-plate drift.
   - **QC Anchors & Unknowns**: Include Low/Mid/High QC standards (`80%–120%` recovery check) and randomized unknown sample dilutions (`dilution_factor` logged per well) for weighted 4PL/5PL `AICc` fitting and `LOD / LLOQ / ULOQ` censoring in `wa_immunoassay_curve`.
