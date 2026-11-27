import assert from "node:assert/strict";
import test from "node:test";
import { ToolResultSchema } from "phi-omics-schemas";
import type { JsonClient } from "phi-omics-runtime";
import { tools } from "../server.js";

const mockClient: JsonClient = {
  async getJson(url: string) {
    if (url.includes("europepmc")) {
      return {
        hitCount: 1,
        resultList: {
          result: [
            {
              id: "PPR800001",
              source: "PPR",
              doi: "10.1101/2026.01.15.600001",
              title: "Spatial multi-omics atlas of liver regeneration",
              authorList: { author: [{ fullName: "Chen Wei" }] },
              journalInfo: { journal: { title: "bioRxiv" } },
              firstPublicationDate: "2026-01-16",
              abstractText: "High-resolution spatial transcriptomics of regenerating liver.",
              isOpenAccess: "Y",
              license: "cc by",
            },
          ],
        },
      };
    }
    if (url.includes("api.biorxiv.org/details/biorxiv/10.1101/2026.01.15.600001")) {
      return {
        collection: [
          {
            doi: "10.1101/2026.01.15.600001",
            title: "Spatial multi-omics atlas of liver regeneration",
            authors: "Chen, W.; Zhang, L.",
            author_corresponding: "Wei Chen",
            author_corresponding_institution: "PKU",
            date: "2026-01-16",
            version: "1",
            type: "new results",
            license: "cc_by",
            category: "genomics",
            jatsxml: "https://www.biorxiv.org/content/early/2026/01/16/2026.01.15.600001.source.xml",
            abstract: "Version 1 abstract.",
            published: "10.1038/s41588-026-00999-1",
          },
          {
            doi: "10.1101/2026.01.15.600001",
            title: "Spatial multi-omics atlas of liver regeneration",
            authors: "Chen, W.; Zhang, L.",
            author_corresponding: "Wei Chen",
            author_corresponding_institution: "PKU",
            date: "2026-03-01",
            version: "2",
            type: "new results",
            license: "cc_by",
            category: "genomics",
            jatsxml: "https://www.biorxiv.org/content/early/2026/03/01/2026.01.15.600001.source.xml",
            abstract: "Revised version 2 abstract.",
            published: "10.1038/s41588-026-00999-1",
          },
        ],
      };
    }
    if (url.includes("api.biorxiv.org/details/biorxiv/")) {
      return {
        messages: [{ total: 1, count: 1 }],
        collection: [
          {
            doi: "10.1101/2026.01.15.600001",
            title: "Spatial multi-omics atlas of liver regeneration",
            authors: "Chen, W.",
            date: "2026-03-01",
            version: "2",
            license: "cc_by",
            category: "genomics",
            abstract: "Revised version 2 abstract.",
          },
        ],
      };
    }
    if (url.includes("api.biorxiv.org/pubs/")) {
      return {
        messages: [{ total: 1, count: 1 }],
        collection: [
          {
            biorxiv_doi: "10.1101/2026.01.15.600001",
            published_doi: "10.1038/s41588-026-00999-1",
            preprint_title: "Spatial multi-omics atlas of liver regeneration",
            preprint_authors: "Chen, W.; Zhang, L.",
            preprint_category: "genomics",
            preprint_date: "2026-01-16",
            published_date: "2026-06-10",
            published_journal: "Nature Genetics",
          },
        ],
      };
    }
    return {};
  },
  async postJson() {
    return {};
  },
};

test("every biorxiv tool returns a valid ToolResult schema", async () => {
  assert.equal(tools.length, 5);
  const calls: Record<string, Record<string, unknown>> = {
    usage_terms: {},
    biorxiv_search_preprints: { query: "spatial transcriptomics liver", server: "biorxiv" },
    biorxiv_browse_recent: { server: "biorxiv", date_from: "2026-03-01", date_to: "2026-03-07" },
    biorxiv_get_preprint: { doi: "10.1101/2026.01.15.600001" },
    biorxiv_get_published: { doi: "10.1101/2026.01.15.600001" },
  };
  for (const tool of tools) {
    const res = ToolResultSchema.parse(await tool.call(mockClient, calls[tool.name] ?? {}));
    assert.equal(res.entity_type, "page", tool.name);
    assert.equal(res.source.id, "biorxiv");
  }

  const detail = ToolResultSchema.parse(
    await tools.find((t) => t.name === "biorxiv_get_preprint")!.call(mockClient, { doi: "10.1101/2026.01.15.600001" }),
  );
  if (detail.entity_type === "page") {
    const pub = detail.items[0];
    assert.equal(pub?.entity_type, "publication");
    if (pub?.entity_type === "publication") {
      assert.equal(pub.published_as, "doi:10.1038/s41588-026-00999-1");
      assert.equal(pub.source_extra.version_count, 2);
    }
  }
});
