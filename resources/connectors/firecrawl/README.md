Firecrawl is a web scraping and crawling engine built for AI workflows, converting single pages or entire websites into clean Markdown and schema-aligned JSON while handling dynamic rendering and link discovery.

### Capabilities

- Search the web and return page titles, URLs, and snippets to pinpoint primary sources.
- Scrape individual URLs into clean Markdown, stripping boilerplate navigation while preserving headings, tables, metadata, and links.
- Map site URL structures and crawl specific path prefixes or domains with configurable depth and page limits.
- Extract structured JSON fields from unstructured web pages using explicit field descriptions or schemas.

### Use Cases

- Index documentation sites, database manuals, and software changelogs to support reproducible analysis workflows.
- Harvest publicly available dataset landing pages, supplementary material indexes, and consortium announcements.
- Track conference programs, funding calls, and consortium updates across institutional websites.
- Extract standardized comparison tables of tools, reagents, cohorts, or benchmarks from multiple web sources.

### Usage Notes

- Requires a Firecrawl API key configured in the connector header settings.
- Page scrapes, crawls, and structured extractions consume credits based on your Firecrawl plan; scope path filters and page limits carefully before crawling large sites.
- Homepage: [Firecrawl MCP Server](https://docs.firecrawl.dev/mcp-server). Respect target websites' `robots.txt` directives, terms of service, and copyright terms.
