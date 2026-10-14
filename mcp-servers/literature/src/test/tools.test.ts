import assert from "node:assert/strict";
import test from "node:test";
import { ToolResultSchema } from "phi-omics-schemas";
import type { JsonClient } from "phi-omics-runtime";
import { tools } from "../server.js";

const longAbstract = "Clustered regularly interspaced short palindromic repeats provide adaptive immunity in bacteria and archaea. ".repeat(12);

function articleXml(pmid: string, doi: string | null, abstract = longAbstract): string {
  return `<PubmedArticle>
    <MedlineCitation><PMID Version="1">${pmid}</PMID>
      <Article>
        <Journal><Title>Science</Title><JournalIssue><PubDate><Year>2012</Year><Month>Aug</Month><Day>17</Day></PubDate></JournalIssue></Journal>
        <ArticleTitle>A programmable dual-RNA-guided DNA endonuclease.</ArticleTitle>
        <Abstract><AbstractText Label="BACKGROUND">${abstract}</AbstractText></Abstract>
        <AuthorList><Author><LastName>Jinek</LastName><ForeName>Martin</ForeName><Identifier Source="ORCID">0000-0002-7601-210X</Identifier></Author></AuthorList>
        <PublicationTypeList><PublicationType>Journal Article</PublicationType></PublicationTypeList>
      </Article>
    </MedlineCitation>
    <PubmedData><ArticleIdList>
      <ArticleId IdType="pubmed">${pmid}</ArticleId>
      ${doi ? `<ArticleId IdType="doi">${doi}</ArticleId>` : ""}
      <ArticleId IdType="pmc">PMC6286148</ArticleId>
    </ArticleIdList></PubmedData>
  </PubmedArticle>`;
}

function work(doi: string, options: { pmcid?: string; openAlex?: string; license?: string | null; title?: string } = {}) {
  const openAlex = options.openAlex ?? "W100";
  return {
    id: `https://openalex.org/${openAlex}`,
    doi: `https://doi.org/${doi}`,
    display_name: options.title ?? "Less can be more for predicting properties",
    publication_date: "2026-10-12",
    type: "article",
    ids: {
      openalex: `https://openalex.org/${openAlex}`,
      doi: `https://doi.org/${doi}`,
      ...(options.pmcid ? { pmcid: options.pmcid } : {}),
    },
    authorships: [{ author: { display_name: "Ada Lovelace", orcid: "https://orcid.org/0000-0001-2345-6789" } }],
    primary_location: { source: { display_name: "Nature Computational Science" }, landing_page_url: `https://doi.org/${doi}` },
    best_oa_location: options.license === null ? null : {
      license: options.license ?? "cc-by",
      url_for_pdf: "https://example.test/paper.pdf",
      url: "https://example.test/paper",
      version: "published",
    },
    open_access: { is_oa: options.license !== null, oa_url: options.license === null ? null : "https://example.test/paper" },
    abstract_inverted_index: { Adaptive: [0], immunity: [1], uses: [2], CRISPR: [3] },
    referenced_works: ["https://openalex.org/W101"],
    related_works: ["https://openalex.org/W102"],
    cited_by_count: 4,
  };
}

const jats = `<?xml version="1.0"?>
<article>
  <front><article-meta>
    <abstract><p>Visitors collaborate with local partners.</p></abstract>
    <license><ali:license_ref>https://creativecommons.org/licenses/by/4.0/</ali:license_ref></license>
  </article-meta></front>
  <body>
    <sec><title>Method</title><p>We visited homes and recorded collaboration partners in interviews.</p>
      <sec><title>Design</title><p>Qualitative study of preventive home visits.</p></sec>
    </sec>
    <sec><title>Results</title><p>Partners included nurses and family members.</p></sec>
  </body>
</article>`;

function doiIn(url: string): string {
  const decoded = decodeURIComponent(url);
  const match = decoded.match(/10\.\d{4,9}\/[^\s&?]+/i);
  return match?.[0]?.replace(/[.,;]+$/, "") ?? "";
}

function client(): JsonClient {
  return {
    async getJson(url: string) {
      if (url.includes("esearch.fcgi")) {
        return { esearchresult: { count: "1", idlist: ["22745249"], querytranslation: "CRISPR[Title]" } };
      }
      if (url.includes("elink.fcgi")) {
        return { linksets: [{ linksetdbs: [{ linkname: "pubmed_pubmed", links: [{ id: "27096362" }] }] }] };
      }
      if (url.includes("europepmc")) {
        if (url.includes("scs.70325")) {
          return {
            hitCount: 1,
            resultList: {
              result: [{
                id: "42828510",
                pmid: "42828510",
                doi: "10.1111/scs.70325",
                pmcid: "PMC13633696",
                source: "MED",
                title: "Collaboration for healthy ageing",
                firstPublicationDate: "2026-01-01",
                isOpenAccess: "Y",
                license: "cc by",
                abstractText: "Visitors collaborate with local partners.",
              }],
            },
          };
        }
        const doiOnly = url.includes("EXT_ID%3A999") || url.includes("EXT_ID:999");
        const withDoi = url.includes("22745249");
        if (!doiOnly && !withDoi) return { resultList: { result: [] }, hitCount: 0 };
        return {
          hitCount: 1,
          resultList: {
            result: [doiOnly ? {
              id: "999",
              pmid: "999",
              source: "MED",
              title: "Online first note without a DOI",
              firstPublicationDate: "2026-10-01",
              isOpenAccess: "N",
              abstractText: "No DOI has been issued.",
              authorList: { author: [{ fullName: "Ada Lovelace" }] },
            } : {
              id: "22745249",
              pmid: "22745249",
              doi: "10.1126/science.1225829",
              pmcid: "PMC6286148",
              source: "MED",
              title: "A programmable dual-RNA-guided DNA endonuclease in adaptive bacterial immunity.",
              firstPublicationDate: "2012-06-28",
              isOpenAccess: "Y",
              license: "cc by",
              abstractText: longAbstract,
              journalInfo: { journal: { title: "Science" } },
              authorList: { author: [{ fullName: "Martin Jinek", authorId: { type: "ORCID", value: "0000-0002-7601-210X" } }] },
            }],
          },
        };
      }
      if (url.includes("api.biorxiv.org/details/biorxiv/10.1101/2024.05.28.596311")) {
        return { messages: [{ status: "ok" }], collection: [{ doi: "10.1101/2024.05.28.596311", title: "Bilateral integration in somatosensory cortex", authors: "Park, H.", date: "2024-06-01", category: "neuroscience", version: "1", license: "cc_by_nc_nd", published: "10.1038/s41593-025-01960-z", abstract: "Sensory perception naturally requires both sides.", type: "new results" }] };
      }
      if (url.includes("api.biorxiv.org/details/medrxiv/")) return { messages: [{ status: "ok" }], collection: [] };
      if (url.includes("api.biorxiv.org/details/biorxiv/")) {
        return { messages: [{ status: "ok", total: 1, count: 1 }], collection: [{ doi: "10.1101/2024.05.28.596311", title: "Bilateral integration", authors: "Park, H.", date: "2024-06-01", category: "neuroscience", version: "1", license: "cc_by", published: "NA", abstract: "Sensory perception.", type: "new results" }] };
      }
      if (url.includes("api.openalex.org/works")) {
        if (url.includes("search=")) return { meta: { count: 1 }, results: [work("10.1126/science.1225829")] };
        if (url.includes("filter=")) {
          if (url.includes("cites%3A") || url.includes("cites:")) return { meta: { count: 2 }, results: [work("10.1126/science.1231143", { openAlex: "W200", title: "Multiplex genome engineering" })] };
          if (url.includes("openalex%3A") || url.includes("openalex:")) return { meta: { count: 1 }, results: [work("10.1000/w101", { openAlex: "W101", title: "Referenced paper" })] };
          return { meta: { count: 1 }, results: [work(doiIn(url) || "10.1126/science.1225829")] };
        }
        const doi = doiIn(url) || "10.1038/s43588-026-01045-2";
        if (doi.includes("10.9999/closed")) return work(doi, { license: null, title: "Closed article" });
        return work(doi, doi.includes("scs.70325") ? { title: "Collaboration for healthy ageing" } : {});
      }
      throw new Error(`unmocked ${url}`);
    },
    async postJson() { return {}; },
    async getText(url: string) {
      if (url.includes("efetch.fcgi")) {
        const ids = new URL(url).searchParams.get("id")?.split(",") ?? [];
        const articles = ids.map((id) => articleXml(id, id === "22745249" ? "10.1126/science.1225829" : null));
        return `<PubmedArticleSet>${articles.join("")}</PubmedArticleSet>`;
      }
      if (url.includes("fullTextXML")) return jats;
      if (url.includes("export.arxiv.org")) {
        return `<feed><entry>
          <id>http://arxiv.org/abs/2201.00001v2</id>
          <title>A q-bio preprint</title>
          <summary>Short summary of a quantitative biology preprint.</summary>
          <published>2022-01-01T00:00:00Z</published>
          <updated>2022-01-02T00:00:00Z</updated>
          <author><name>A Author</name></author>
          <link title="pdf" href="https://arxiv.org/pdf/2201.00001"/>
        </entry></feed>`;
      }
      throw new Error(`unmocked text ${url}`);
    },
  };
}

const args: Record<string, Record<string, unknown>> = {
  usage_terms: {},
  search_publications: { query: "CRISPR[ti]", limit: 2 },
  get_publication: { id: "10.1038/s43588-026-01045-2" },
  get_fulltext: { id: "10.1111/scs.70325", section: "methods", max_chars: 800 },
  related_publications: { id: "10.1126/science.1225829", relation: "citations", limit: 2 },
  resolve_references: { references: ["10.1126/science.1225829", "PMID:22745249"] },
};

test("every literature tool returns the shared schema", async () => {
  const http = client();
  assert.equal(tools.length, 6);
  for (const tool of tools) {
    const result = ToolResultSchema.parse(await tool.call(http, args[tool.name] ?? {}));
    assert.equal(result.entity_type, "page", `${tool.name} ${JSON.stringify(result).slice(0, 400)}`);
    if (result.entity_type === "page" && tool.name !== "usage_terms") {
      assert.equal(result.item_type, "publication");
      assert.ok(result.items.length > 0, tool.name);
    }
  }
});

test("search keeps the DOI as the id and truncates the abstract", async () => {
  const result = ToolResultSchema.parse(await tools.find((tool) => tool.name === "search_publications")!.call(client(), { query: "CRISPR[ti]", limit: 2 }));
  assert.equal(result.entity_type, "page");
  if (result.entity_type !== "page" || result.items[0]?.entity_type !== "publication") return;
  assert.equal(result.items[0].id, "doi:10.1126/science.1225829");
  assert.equal(result.items[0].doi, "10.1126/science.1225829");
  assert.ok(result.items[0].xrefs.includes("pmid:22745249"));
  assert.ok((result.items[0].abstract?.length ?? 0) <= 400);
  assert.equal(result.items[0].source_extra.abstract_truncated, true);
  assert.equal(result.items[0].authors[0]?.orcid, "https://orcid.org/0000-0002-7601-210X");
  assert.equal(result.source_extra.query_translation, "CRISPR[Title]");
});

test("a DOI-only article uses doi: and does not require a PMID", async () => {
  const result = ToolResultSchema.parse(await tools.find((tool) => tool.name === "get_publication")!.call(client(), { id: "HTTPS://DOI.ORG/10.1038/S43588-026-01045-2" }));
  assert.equal(result.entity_type, "page");
  if (result.entity_type !== "page" || result.items[0]?.entity_type !== "publication") return;
  assert.equal(result.items[0].id, "doi:10.1038/s43588-026-01045-2");
  assert.equal(result.items[0].xrefs.some((id) => id.startsWith("pmid:")), false);
  assert.equal(result.items[0].authors[0]?.orcid, "https://orcid.org/0000-0001-2345-6789");
});

test("a PMID article keeps the DOI as the primary key", async () => {
  const result = ToolResultSchema.parse(await tools.find((tool) => tool.name === "get_publication")!.call(client(), { id: "pmid:22745249" }));
  assert.equal(result.entity_type, "page");
  if (result.entity_type !== "page" || result.items[0]?.entity_type !== "publication") return;
  assert.equal(result.items[0].id, "doi:10.1126/science.1225829");
  assert.ok(result.items[0].xrefs.includes("pmid:22745249"));
  assert.ok((result.items[0].abstract?.length ?? 0) > 400);
});

test("a record with no DOI falls back to pmid:", async () => {
  const result = ToolResultSchema.parse(await tools.find((tool) => tool.name === "get_publication")!.call(client(), { id: "pmid:999" }));
  assert.equal(result.entity_type, "page");
  if (result.entity_type !== "page" || result.items[0]?.entity_type !== "publication") return;
  assert.equal(result.items[0].id, "pmid:999");
  assert.equal(result.items[0].doi, null);
});

test("a bioRxiv preprint links the published DOI without merging the records", async () => {
  const result = ToolResultSchema.parse(await tools.find((tool) => tool.name === "get_publication")!.call(client(), { id: "10.1101/2024.05.28.596311" }));
  assert.equal(result.entity_type, "page");
  if (result.entity_type !== "page" || result.items[0]?.entity_type !== "publication") return;
  assert.equal(result.items[0].id, "doi:10.1101/2024.05.28.596311");
  assert.equal(result.items[0].type, "preprint");
  assert.equal(result.items[0].published_as, "doi:10.1038/s41593-025-01960-z");
  assert.equal(result.items[0].license, "cc-by-nc-nd");
  assert.equal(result.items[0].preprint_of, null);
});

test("an arXiv id is the primary key when the preprint has no DOI", async () => {
  const result = ToolResultSchema.parse(await tools.find((tool) => tool.name === "get_publication")!.call(client(), { id: "arxiv:2201.00001v2" }));
  assert.equal(result.entity_type, "page");
  if (result.entity_type !== "page" || result.items[0]?.entity_type !== "publication") return;
  assert.equal(result.items[0].id, "arxiv:2201.00001");
  assert.equal(result.items[0].type, "preprint");
});

test("open full text inlines CC-BY methods and withholds text when the license is unknown", async () => {
  const open = ToolResultSchema.parse(await tools.find((tool) => tool.name === "get_fulltext")!.call(client(), { id: "10.1111/scs.70325", section: "methods", max_chars: 800 }));
  assert.equal(open.entity_type, "page");
  if (open.entity_type === "page") {
    assert.equal(open.source_extra.fulltext, "sections");
    assert.equal(open.items[0]?.entity_type === "publication" && open.items[0].xrefs.includes("pmcid:PMC13633696"), true);
    const sections = open.source_extra.sections as Array<{ heading: string; text: string }>;
    assert.ok(sections.some((section) => /method/i.test(section.heading)));
    assert.ok(sections.some((section) => section.text.includes("interviews") || section.text.includes("home")));
  }
  const closed = ToolResultSchema.parse(await tools.find((tool) => tool.name === "get_fulltext")!.call(client(), { id: "10.9999/closed", section: "methods" }));
  assert.equal(closed.entity_type, "page");
  if (closed.entity_type === "page") {
    assert.equal(closed.source_extra.fulltext, "links");
    assert.deepEqual(closed.source_extra.sections, []);
  }
});

test("requested sources merge duplicate DOIs", async () => {
  const result = ToolResultSchema.parse(await tools.find((tool) => tool.name === "search_publications")!.call(client(), {
    query: "CRISPR",
    sources: ["pubmed", "openalex"],
    limit: 5,
  }));
  assert.equal(result.entity_type, "page");
  if (result.entity_type !== "page") return;
  assert.equal(result.items.length, 1);
  if (result.items[0]?.entity_type === "publication") assert.equal(result.items[0].id, "doi:10.1126/science.1225829");
});

test("bioRxiv date browse and invalid routing return schema objects", async () => {
  const browse = ToolResultSchema.parse(await tools.find((tool) => tool.name === "search_publications")!.call(client(), {
    type: "preprint",
    server: "biorxiv",
    date_from: "2024-06-01",
    date_to: "2024-06-02",
    limit: 5,
  }));
  assert.equal(browse.entity_type, "page");
  if (browse.entity_type === "page" && browse.items[0]?.entity_type === "publication") {
    assert.equal(browse.items[0].type, "preprint");
    assert.equal(browse.items[0].published_as, null);
  }
  const missing = ToolResultSchema.parse(await tools.find((tool) => tool.name === "search_publications")!.call(client(), {}));
  assert.equal(missing.entity_type, "error");
  if (missing.entity_type === "error") assert.equal(missing.code, "invalid_input");
  const keyword = ToolResultSchema.parse(await tools.find((tool) => tool.name === "search_publications")!.call(client(), { query: "alzheimers", sources: ["biorxiv"] }));
  assert.equal(keyword.entity_type, "error");
  if (keyword.entity_type === "error") assert.equal(keyword.code, "invalid_input");
});
