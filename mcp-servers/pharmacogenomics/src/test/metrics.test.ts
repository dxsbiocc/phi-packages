import assert from "node:assert/strict";
import test from "node:test";
import { metricInfo, moreSensitive } from "../metrics.js";

test("potency, LFC and viability metrics are lower-is-stronger in every dataset", () => {
  for (const metric of ["LN_IC50", "IC50", "EC50", "AC50", "GI50", "TGI", "LC50", "LFC", "VIABILITY"]) {
    assert.equal(metricInfo("anything", metric).direction, "lower_is_stronger", metric);
  }
  assert.equal(metricInfo("gdsc1", "ln_ic50").direction, "lower_is_stronger", "the metric name is case-insensitive");
});

test("AAC is higher-is-stronger", () => {
  assert.equal(metricInfo("ctrpv2", "AAC").direction, "higher_is_stronger");
});

test("AUC direction is only asserted for datasets that document it", () => {
  for (const dataset of ["gdsc1", "gdsc2", "prism_secondary_19q4"]) {
    assert.equal(metricInfo(dataset, "AUC").direction, "lower_is_stronger", dataset);
  }
  for (const dataset of ["ctrpv2", "gcsi", "ccle", "fimm", null]) {
    const info = metricInfo(dataset, "AUC");
    assert.equal(info.direction, "unknown", String(dataset));
    assert.match(info.note ?? "", /not asserted/);
  }
});

test("unrecognised and combination metrics are unknown", () => {
  assert.equal(metricInfo("gdsc_combination_original", "SYNERGY_DELTA_AUC").direction, "unknown");
  assert.equal(metricInfo("x", "PTC").direction, "unknown");
  assert.equal(metricInfo(null, null).direction, "unknown");
});

test("moreSensitive reads the direction against which side has the larger value", () => {
  assert.equal(moreSensitive("lower_is_stronger", true), "second", "the larger-valued group is less sensitive");
  assert.equal(moreSensitive("lower_is_stronger", false), "first");
  assert.equal(moreSensitive("higher_is_stronger", true), "first");
  assert.equal(moreSensitive("higher_is_stronger", false), "second");
  assert.equal(moreSensitive("unknown", true), null);
});
