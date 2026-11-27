# Opentrons Flex & OT-2 API v2 Guide (`2.16+`–`2.20+`)

Load this reference when generating or auditing **Opentrons Python Protocol API v2** scripts for `OT-2` or `Opentrons Flex`.

## 1. Platform & API Header Differences

| Feature | Opentrons OT-2 (`apiLevel >= 2.16`) | Opentrons Flex (`apiLevel >= 2.16`–`2.20+`) |
| --- | --- | --- |
| **Header declaration** | `metadata` + `requirements = {"robotType": "OT-2", "apiLevel": "2.18"}` | `requirements = {"robotType": "Flex", "apiLevel": "2.20"}` |
| **Deck coordinates** | Numeric strings `"1"`–`"11"` + fixed `"12"` trash | Alphanumeric `"A1"`–`"D3"` (active) + `"A4"`–`"D4"` staging slots |
| **Trash / waste** | Built-in fixed trash in slot `"12"` | Explicit `protocol.load_trash_bin("A3")` or `protocol.load_waste_chute()` |
| **Pipette naming** | `p20_single_gen2`, `p300_multi_gen2`, `p1000_single_gen2` | `flex_1channel_50`, `flex_8channel_50`, `flex_1channel_1000`, `flex_8channel_1000`, `flex_96channel_1000` |
| **Dynamic range (`>= 10%`)** | P20: `2–20 uL` (min `1 uL` only w/ verified calibration); P300: `30–300 uL`; P1000: `100–1000 uL` | Flex 50: `5–50 uL`; Flex 1000: `100–1000 uL` (or `50 uL` tip w/ Flex 1000 only down to `5 uL` if 50 uL tips loaded) |
| **Plate transport** | Manual pause (`protocol.pause(...)`) | Automated Gripper (`protocol.move_labware(plate, dest, use_gripper=True)`) |

> **Never mix OT-2 and Flex tip racks or deck coordinates.** On Flex, always match `opentrons_flex_96_tiprack_50ul` / `200ul` / `1000ul` (or filter equivalents) to the pipette's active `tip_racks`.

## 2. Viscous, Volatile & Cell-Suspension Liquid Classes

Default Opentrons speeds (`100–150 uL/s`) cause under-aspiration in glycerol/master mix, dripping in ethanol/DMSO, and shearing in live cells. Override flow rates and post-aspiration behavior explicitly:

| Liquid Class | Examples | Aspirate / Dispense Rate | Post-Aspirate & Dispense Rules |
| --- | --- | --- | --- |
| **Aqueous / Buffer** | PBS, water, TE | `1.0x` default (`~94–150 uL/s`) | `air_gap(5)` for multi-dispense; `blow_out` optional in waste/dest |
| **Viscous** | qPCR 2x Master Mix, 50% glycerol, PEG, enzyme stock | `0.25x–0.4x` (`25–40%` default) | `protocol.delay(seconds=2)` after aspirate & dispense; `aspirate` `1–2 mm` above bottom; `touch_tip(radius=0.75, v_offset=-2)` or slow tip withdrawal; `blow_out(dest.top(-2))` |
| **Volatile** | 70–80% EtOH, ACN, MeOH | `0.6x–0.8x` | Pre-wet tip (`2x` aspirate/dispense in source); mandatory `air_gap(10–20)` before travel; dispense fast into dest, no `blow_out` splash |
| **Live Cells / Beads** | Cell suspension, SPRI magnetic beads | `0.3x–0.5x` | `mix(3, 0.7 * well_vol)` gently before each aspiration using wide-bore or `1000 uL` tips; avoid `blow_out` bubbles in cell wells |

## 3. Safe Module & Gripper Sequencing Rules

1. **Heater-Shaker (`heaterShakerModuleV1`)**:
   - **Always** `hs_mod.close_labware_latch()` before `set_and_wait_for_shake_speed(...)`.
   - **Always** `hs_mod.deactivate_shaker()` and `hs_mod.open_labware_latch()` before `move_labware(..., use_gripper=True)` or manual removal.
   - Adjacent-slot clearance: Never place tall tip racks (`1000 uL` / `300 uL`) or tall tube racks immediately left/right of the Heater-Shaker when using `8-channel` / `96-channel` pipettes unless verified clear in `opentrons_simulate`.
2. **Thermocycler (`thermocyclerModuleV2`)**:
   - Occupies slots `7, 8, 10, 11` (OT-2) or `A1, B1` (Flex).
   - **Always** `tc_mod.open_lid()` before pipetting or Gripper `move_labware`.
   - When heating `>= 65 C` in sealed plates, use `tc_mod.set_lid_temperature(105)` and `tc_mod.close_lid()` before thermal cycling.
3. **Magnetic Block (`magneticBlockV1` on Flex) / Magnetic Module (`magnetic module gen2` on OT-2)**:
   - On Flex, `magneticBlockV1` is passive—use `protocol.move_labware(plate, mag_block, use_gripper=True)` to engage beads and move off-magnet for resuspension.
   - Aspirate supernatant slowly (`0.2x–0.3x`) at an offset away from the bead pellet (`well.bottom(z=0.5).move(types.Point(x=1.5, y=0, z=0))`) to avoid bead carryover.

## 4. Annotated Opentrons Flex API v2.20 Pattern

```python
from opentrons import protocol_api

requirements = {
    "robotType": "Flex",
    "apiLevel": "2.20",
}

metadata = {
    "protocolName": "qPCR Reaction Setup with Viscous Master Mix",
    "author": "Lab Automation Skill",
    "description": "Master-mix-first distribution + per-sample template transfer",
}


def assert_volume_budget(required_ul: float, loaded_ul: float, dead_ul: float, label: str) -> None:
    """Fail fast during simulation if source fill violates dead-volume floor."""
    if loaded_ul - required_ul < dead_ul:
        raise ValueError(
            f"[DEAD VOLUME VIOLATION] {label}: need {required_ul:.1f} uL + {dead_ul:.1f} uL dead vol, "
            f"but only {loaded_ul:.1f} uL loaded."
        )


def run(protocol: protocol_api.ProtocolContext) -> None:
    # Deck setup
    trash = protocol.load_trash_bin("A3")
    tips_50 = [protocol.load_labware("opentrons_flex_96_filtertiprack_50ul", slot) for slot in ["C1", "C2"]]
    tips_1000 = [protocol.load_labware("opentrons_flex_96_filtertiprack_1000ul", "B1")]

    reagent_plate = protocol.load_labware("nest_96_wellplate_2ml_deep", "D1")
    sample_plate = protocol.load_labware("nest_96_wellplate_100ul_pcr_full_skirt", "D2")

    hs_mod = protocol.load_module("heaterShakerModuleV1", "D3")
    hs_adapter = hs_mod.load_adapter("opentrons_96_pcr_adapter")
    dest_plate = hs_adapter.load_labware("nest_96_wellplate_100ul_pcr_full_skirt")

    p50 = protocol.load_instrument("flex_1channel_50", "left", tip_racks=tips_50)
    p1000 = protocol.load_instrument("flex_1channel_1000", "right", tip_racks=tips_1000)

    # Pre-flight volume ledger check (24 reactions x 15 uL master mix = 360 uL)
    mm_per_rxn = 15.0
    n_rxns = 24
    assert_volume_budget(mm_per_rxn * n_rxns, loaded_ul=450.0, dead_ul=50.0, label="2x Master Mix (A1)")

    # 1. Master mix distribution (viscous liquid class, clean empty dest wells)
    hs_mod.open_labware_latch()  # Ensure accessible; close before shaking
    p50.flow_rate.aspirate = 12  # ~0.35x slow aspirate for viscous 2x MM
    p50.flow_rate.dispense = 15
    p50.pick_up_tip()
    for well in dest_plate.wells()[:n_rxns]:
        p50.aspirate(mm_per_rxn, reagent_plate["A1"].bottom(z=1.5))
        protocol.delay(seconds=1.5)
        p50.air_gap(5.0)
        p50.dispense(mm_per_rxn + 5.0, well.bottom(z=1.0))
        p50.blow_out(well.top(z=-2.0))
        p50.touch_tip(well, radius=0.75, v_offset=-2.0)
    p50.drop_tip()

    # 2. Sample template transfer (5.0 uL >= 10% of 50 uL tip, strict new_tip="always")
    p50.flow_rate.aspirate = 35
    p50.flow_rate.dispense = 35
    for src, dst in zip(sample_plate.wells()[:n_rxns], dest_plate.wells()[:n_rxns]):
        p50.pick_up_tip()
        p50.aspirate(5.0, src.bottom(z=1.0))
        p50.dispense(5.0, dst.bottom(z=1.5))
        p50.mix(3, 12.0, dst.bottom(z=1.5))
        p50.blow_out(dst.top(z=-2.0))
        p50.drop_tip()

    # 3. Safe sealed mixing on Heater-Shaker (after manual or gripper seal/transfer)
    hs_mod.close_labware_latch()
    hs_mod.set_and_wait_for_shake_speed(1200)
    protocol.delay(seconds=20)
    hs_mod.deactivate_shaker()
    hs_mod.open_labware_latch()
```

## 5. Static Simulation & Pre-Bench Audit Checklist

Run Opentrons' static simulator locally before touching the deck:

```bash
opentrons_simulate protocol.py -e
```

Verify every item below passes without warnings:
- [ ] **API & Robot Match**: `requirements["robotType"]` (`OT-2` vs `Flex`) and `apiLevel` (`>= 2.16`) match the target robot software.
- [ ] **No Sub-10% Aspirations**: Every `aspirate(vol)` is `>= 0.10 * tip_max_vol` (check intermediate dilution steps if `< 1 uL` is needed).
- [ ] **Dead-Volume Buffer Verified**: Source wells retain `>= 10–15 uL` (96-well V-bottom), `>= 3–5 uL` (384-well), or `>= 10%–15%` trough overage after the final transfer.
- [ ] **Tip Inventory Sufficient**: Total `pick_up_tip()` calls `<=` loaded tip count (`96 * n_racks`), with `new_tip="always"` on all sample/template transfers.
- [ ] **Module Latch & Lid Guards**: Heater-Shaker latch closes before shaking and opens before Gripper moves; Thermocycler lid opens before pipetting/Gripper access.
- [ ] **Deck & Height Clearance**: Tall labware/modules do not sit between an 8-/96-channel gantry move and short plates.
