import assert from "node:assert/strict";
import { mkdtemp, rm } from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import test from "node:test";
import { ToolResultSchema, type ToolResult } from "phi-omics-schemas";
import { UpstreamClient, UpstreamError, loadConfig } from "phi-omics-runtime";
import { tools } from "../server.js";

const dir = await mkdtemp(path.join(os.tmpdir(), "phi-literature-live-"));
process.env.PHI_OMICS_CACHE_DIR = dir;
process.env.PHI_OMICS_RATE_LIMIT_RPS = "2";
process.env.PHI_OMICS_TIMEOUT_MS = "30000";
process.env.PHI_LITERATURE_EMAIL ??= "phi-literature@users.noreply.github.com";
const client = new UpstreamClient(loadConfig());

async function live(t: test.TestContext, name: string, args: Record<string, unknown>, check: (result: ToolResult) => void) {
  if (process.env.PHI_MCP_SKIP_LIVE === "1") {
    t.skip("PHI_MCP_SKIP_LIVE=1");
    return;
  }
  const found = tools.find((item) => item.name === name);
  if (!found) throw new Error(name);
  try {
    const result = ToolResultSchema.parse(await found.call(client, args));
    if (result.entity_type === "error" && (result.code === "unavailable" || result.retryable)) {
      t.skip(result.message);
      return;
    }
    check(result);
  } catch (error) {
    if (error instanceof UpstreamError && error.status == null) {
      t.skip(error.message);
      return;
    }
    throw error;
  }
}

test("live DOI-only online-first article", async (t) => {
  await live(t, "get_publication", { id: "10.1038/s43588-026-01045-2" }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type !== "page" || result.items[0]?.entity_type !== "publication") return;
    assert.equal(result.items[0].id, "doi:10.1038/s43588-026-01045-2");
    assert.equal(result.items[0].xrefs.some((id) => id.startsWith("pmid:")), false);
    assert.ok(result.items[0].title);
  });
});

test("live PMID article uses the DOI as the primary key", async (t) => {
  await live(t, "get_publication", { id: "pmid:22745249" }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type !== "page" || result.items[0]?.entity_type !== "publication") return;
    assert.equal(result.items[0].id, "doi:10.1126/science.1225829");
    assert.ok(result.items[0].xrefs.includes("pmid:22745249"));
    assert.ok(result.items[0].abstract);
  });
});

test("live bioRxiv preprint links its published version", async (t) => {
  await live(t, "get_publication", { id: "10.1101/2024.05.28.596311" }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type !== "page" || result.items[0]?.entity_type !== "publication") return;
    assert.equal(result.items[0].type, "preprint");
    assert.equal(result.items[0].published_as, "doi:10.1038/s41593-025-01960-z");
  });
});

test("live PubMed search is recent-first and truncates abstracts", async (t) => {
  await live(t, "search_publications", { query: "CRISPR[ti] AND 2024:2026[pdat]", limit: 2 }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type !== "page") return;
    assert.equal(result.item_type, "publication");
    assert.ok(result.items.length > 0);
    const first = result.items[0];
    if (first?.entity_type === "publication" && first.abstract) assert.ok(first.abstract.length <= 400);
  });
});

test("live citation expansion and CC-BY full text", async (t) => {
  await live(t, "related_publications", { id: "doi:10.1126/science.1225829", relation: "citations", limit: 2 }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") {
      assert.ok(result.items.length > 0);
      assert.equal(result.items[0]?.entity_type, "publication");
    }
  });
  await live(t, "get_fulltext", { id: "10.1111/scs.70325", section: "methods", max_chars: 1500 }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type !== "page") return;
    assert.equal(result.source_extra.fulltext, "sections");
    const sections = result.source_extra.sections as Array<{ text: string }>;
    assert.ok(sections.some((section) => section.text.length > 40));
  });
});

test("live reference list resolves a DOI and a PMID", async (t) => {
  await live(t, "resolve_references", { references: ["10.1038/s43588-026-01045-2", "PMID:22745249"] }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type !== "page") return;
    assert.equal(result.items.length, 2);
    assert.equal(result.items[0]?.id, "doi:10.1038/s43588-026-01045-2");
    assert.equal(result.items[1]?.id, "doi:10.1126/science.1225829");
  });
});

test.after(async () => rm(dir, { recursive: true, force: true }));
