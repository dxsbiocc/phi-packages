import assert from "node:assert/strict";
import test from "node:test";
import { ToolResultSchema } from "phi-omics-schemas";
import type { JsonClient } from "phi-omics-runtime";
import { tools } from "../server.js";

const mockClient: JsonClient = {
  async getJson(url: string) {
    if (url.includes("esearch.fcgi")) {
      return {
        esearchresult: {
          count: "1",
          idlist: ["38000001"],
          querytranslation: "TP53[All Fields]",
        },
      };
    }
    if (url.includes("elink.fcgi")) {
      return { linksets: [{ linksetdbs: [{ linkname: "pubmed_pubmed", links: [{ id: "38000001" }] }] }] };
    }
    if (url.includes("pubtator3-api/publications/export/biocjson")) {
      return {
        PubTator3: [
          {
            id: "38000001",
            pmid: 38000001,
            passages: [
              {
                text: "TP53 guards genomic stability in hepatocytes.",
                annotations: [{ text: "TP53", infons: { type: "Gene", identifier: "7157", name: "TP53" } }],
              },
            ],
          },
        ],
      };
    }
    if (url.includes("europepmc")) {
      return {
        hitCount: 1,
        resultList: {
          result: [
            {
              id: "38000001",
              source: "MED",
              pmid: "38000001",
              pmcid: "PMC9999999",
              doi: "10.1038/s41586-026-00001-1",
              title: "TP53 regulation in human liver cancer.",
              isOpenAccess: "Y",
              license: "cc by",
            },
          ],
        },
      };
    }
    if (url.includes("openalex.org")) {
      return {
        id: "https://openalex.org/W38000001",
        doi: "https://doi.org/10.1038/s41586-026-00001-1",
        display_name: "TP53 regulation in human liver cancer.",
        publication_date: "2026-01-10",
        type: "article",
        ids: {
          openalex: "https://openalex.org/W38000001",
          doi: "https://doi.org/10.1038/s41586-026-00001-1",
          pmid: "https://pubmed.ncbi.nlm.nih.gov/38000001",
          pmcid: "PMC9999999",
        },
        authorships: [{ author: { display_name: "Mei Lin" } }],
        primary_location: { source: { display_name: "Nature" }, landing_page_url: "https://doi.org/10.1038/s41586-026-00001-1" },
        best_oa_location: { license: "cc-by", url: "https://pmc.ncbi.nlm.nih.gov/articles/PMC9999999/", version: "published" },
        open_access: { is_oa: true, oa_url: "https://pmc.ncbi.nlm.nih.gov/articles/PMC9999999/" },
        referenced_works: [],
        related_works: [],
      };
    }
    return {};
  },
  async postJson() {
    return {};
  },
  async getText(url: string) {
    if (url.includes("efetch.fcgi")) {
      return `<?xml version="1.0"?>
        <PubmedArticleSet>
          <PubmedArticle>
            <MedlineCitation>
              <PMID Version="1">38000001</PMID>
              <Article>
                <Journal><Title>Nature</Title><JournalIssue><PubDate><Year>2026</Year></PubDate></JournalIssue></Journal>
                <ArticleTitle>TP53 regulation in human liver cancer.</ArticleTitle>
                <Abstract><AbstractText>TP53 guards genomic stability in hepatocytes.</AbstractText></Abstract>
                <AuthorList><Author><LastName>Lin</LastName><ForeName>Mei</ForeName></Author></AuthorList>
                <PublicationTypeList><PublicationType>Journal Article</PublicationType></PublicationTypeList>
              </Article>
            </MedlineCitation>
            <PubmedData>
              <ArticleIdList>
                <ArticleId IdType="pubmed">38000001</ArticleId>
                <ArticleId IdType="pmc">PMC9999999</ArticleId>
                <ArticleId IdType="doi">10.1038/s41586-026-00001-1</ArticleId>
              </ArticleIdList>
            </PubmedData>
          </PubmedArticle>
        </PubmedArticleSet>`;
    }
    if (url.includes("fullTextXML")) {
      return `<?xml version="1.0"?><article><front><article-meta><license><ali:license_ref>https://creativecommons.org/licenses/by/4.0/</ali:license_ref></license></article-meta></front><body><sec><title>Methods</title><p>Single-cell RNA sequencing was performed.</p></sec></body></article>`;
    }
    return "";
  },
};

test("every pubmed tool returns a valid ToolResult schema", async () => {
  assert.equal(tools.length, 7);
  const sampleCalls: Record<string, Record<string, unknown>> = {
    usage_terms: {},
    pubmed_search_articles: { query: "TP53 AND liver cancer" },
    pubmed_get_article: { id: "pmid:38000001" },
    pubmed_get_fulltext: { id: "pmid:38000001", section: "methods" },
    pubmed_related_articles: { id: "pmid:38000001", relation: "similar" },
    pubmed_lookup_citations: { references: ["PMID:38000001"] },
    pubmed_entity_annotations: { id: "pmid:38000001" },
  };
  for (const tool of tools) {
    const res = ToolResultSchema.parse(await tool.call(mockClient, sampleCalls[tool.name] ?? {}));
    assert.equal(res.entity_type, "page", `${tool.name} expected page`);
    assert.equal(res.source.id, "pubmed");
  }
});
