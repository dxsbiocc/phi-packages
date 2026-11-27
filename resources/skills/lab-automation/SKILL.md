---
name: lab-automation
description: "Design, audit, and simulate liquid-handling protocols and anti-edge-bias plate layouts for Opentrons Python Protocol API v2 (OT-2 and Opentrons Flex, API 2.16+ to 2.20+) and PyLabRobot LiquidHandler workflows. Enforces physical dead-volume floors, >= 10% pipette dynamic range, cross-contamination tip policies, randomized or Latin-square control placement compatible with pgx_hts_qc, wa_qpcr_quant, and wa_immunoassay_curve, and static opentrons_simulate dry-run verification. Use for 实验自动化, 移液工作站, Opentrons, PyLabRobot, 板布局设计, liquid handling, plate layout, and serial dilution. Not for downstream qPCR/ELISA curve fitting (use wetlab-assay) or HTS hit scoring (use pharmacogenomics)."
license: MIT
metadata: {"version": "1.0.0", "category": "科学分析", "license": "MIT"}
---

# Liquid-Handling Automation & Plate-Layout Engineering

Generate, audit, and dry-run benchtop liquid-handling scripts and bias-resistant plate maps for **Opentrons Python Protocol API v2** (`OT-2` & `Opentrons Flex`) and hardware-agnostic **PyLabRobot** `LiquidHandler` setups.

## 1. Scope & Boundary

- **Use for**:
  - **Protocol Generation & Auditing**: Writing or reviewing Opentrons API v2 (`2.16+`–`2.20+`) and PyLabRobot scripts with labware definitions, Gripper/Thermocycler/Heater-Shaker/Magnetic Block sequences, tip tracking, and viscous-liquid handling.
  - **Physical & Volumetric Verification**: Enforcing dead-volume reserves, pipette dynamic-range bounds, tip-change hygiene, and static simulation checklists (`opentrons_simulate` / `ChatterboxBackend`).
  - **Plate-Layout Design**: Designing 96/384-well serial dilutions, checkerboard synergy matrices, qPCR master-mix/template layouts, and immunoassay standard curves that export clean plate maps for downstream Phi analysis tools.
- **Not for (Handoff)**:
  - High-level experimental unit, biological power, or animal/clinical study design -> `study-design`.
  - Downstream HTS plate QC ($Z'$-factor, SSMD, `B-score`) -> `pharmacogenomics` (`pgx_hts_qc`).
  - Downstream RT-qPCR quantification or ELISA/BCA 4PL/5PL curve fitting -> `wetlab-assay` (`wa_qpcr_quant`, `wa_immunoassay_curve`).

## 2. Physical & Experimental Stop Gates (Non-Negotiable)

Never emit or approve an executable liquid-handling script that violates any of these five gates:

1. **Dead-Volume & Overage Floor**:
   - Every source well/tube MUST retain physical dead volume (`>= 10–15 uL` in 96-well V-bottom plates, `>= 15–25 uL` in 96-well flat/U-bottom plates, `>= 3–5 uL` in 384-well plates, and `>= 10%–15%` overage in reagent troughs/reservoirs).
   - Required source fill = `(total aspirated volume) + dead_volume + evaporation/pipetting overage`. Halt and flag `[DEAD VOLUME VIOLATION]` if cumulative aspirations deplete a well below its floor.
2. **Pipette Dynamic Range (`>= 10%` Nominal Tip Rule)**:
   - Never aspirate below `10%` of a pipette/tip's nominal maximum volume (e.g., never aspirate `0.5 uL` or `20 uL` with a `P1000` / `Flex 1000 uL`; route `1–20 uL` to `P20` or `5–50 uL` to `Flex 50 uL`).
   - If a target transfer is `< 1 uL` (or below the smallest mounted pipette's `10%` floor), **halt** and either switch to an intermediate working dilution (`>= 5 uL` transfer) or flag `[SUB-RANGE VOLUME — DILUTION REQUIRED]`.
3. **Cross-Contamination Tip-Change Policy**:
   - Default to `new_tip="always"` for biological samples, templates, primers, antibodies, cells, and serial dilution transfers.
   - `new_tip="once"` (`distribute` / multi-dispense) is permitted ONLY when dispensing a sterile bulk reagent (buffer, media, master mix) into **empty** wells *before* any sample/template is present, using an air gap + disposal volume and zero immersion (`dispense` above liquid or `top` offset).
4. **Anti-Edge-Bias & Randomized / Latin-Square Controls**:
   - Never cluster all vehicle/positive controls in Column 1/12 or Row A/H unless an explicit edge-evaporation moat (PBS/water in perimeter wells) is configured AND controls are still stratified/interleaved across interior blocks.
   - Use randomized, balanced, or Latin-square placement for replicates and treatments so plate row/column drift does not confound biological arms—matching the spatial QC assumptions of `pgx_hts_qc`, `wa_qpcr_quant`, and `wa_immunoassay_curve`.
5. **Mandatory Dry-Run Simulation Gate**:
   - No protocol is "ready for bench" until verified via `opentrons_simulate` (API v2) or PyLabRobot `ChatterboxBackend` + deck-clearance audit. Always output the exact simulation command and pre-run checklist.

## 3. On-Demand Reference Routing

Load ONLY the reference file needed for the user's request:
- **Opentrons OT-2 / Flex API v2 (`2.16+`–`2.20+`), Modules, Gripper & Viscous Liquid Classes** -> read [references/opentrons-flex-ot2-guide.md](references/opentrons-flex-ot2-guide.md).
- **PyLabRobot `LiquidHandler` & 96/384-Well Plate-Layout Templates (Serial Dilution, Checkerboard Synergy, qPCR/Immunoassay)** -> read [references/pylabrobot-plate-layouts.md](references/pylabrobot-plate-layouts.md).

## 4. Core Workflow

1. **Lock Hardware & Assay Context**: Confirm robot (`OT-2`, `Opentrons Flex`, or PyLabRobot backend), mounted pipettes & tip racks, deck slots/modules (Thermocycler, Heater-Shaker, Magnetic Block, Gripper), labware load names, and liquid classes (aqueous, DMSO, glycerol/master-mix, volatile ethanol, cell suspension).
2. **Compute Volume Ledger & Check Redlines**: Sum per-well demand across all steps, add plate dead volumes + `10%–15%` trough overage, verify every transfer falls within `[10%, 100%]` of the chosen pipette's range, and verify tip-rack capacity (total tips needed vs loaded boxes).
3. **Design Anti-Edge-Bias Plate Map**: Assign perimeter moat (if incubated), distribute controls via Latin-square or stratified blocks, equalize final solvent `%` (e.g., `<= 0.5%` DMSO) across dose-response/synergy wells, and emit a tidy CSV/JSON well-metadata table ready for `pgx_hts_qc`, `wa_qpcr_quant`, or `wa_immunoassay_curve`.
4. **Emit Protocol & Simulation Verification**: Write clean, modular Python (`run(protocol)` for Opentrons or `async def main()` for PyLabRobot) with explicit deck-setup instructions, pre-run reagent fill table, and `opentrons_simulate` / dry-run checklist.

## 5. Output Contract

Every `lab-automation` deliverable must include:
1. **Deck & Reagent Fill Table**: Exact slot/coordinate mapping, labware load names, tip-rack count, and source volumes including dead-volume + overage margins.
2. **Audited Protocol Script** (or **Audit Findings Table** when reviewing an existing script, classifying blockers vs warnings).
3. **Plate Map Table / Schema** (when designing assays): Well-to-condition mapping compatible with downstream Phi QC/quantification tools.
4. **Dry-Run Verification Commands & Bench Checklist**: Exact `opentrons_simulate` command, physical clearance checks (Heater-Shaker latch, Thermocycler lid, Gripper collision zones), and post-run visual tip/volume checks.
