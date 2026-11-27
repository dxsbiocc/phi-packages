Browser Use runs a cloud-hosted headless browser that lets agents navigate live websites, click links, fill out forms, paginate through listings, and extract structured content from dynamic pages.

### Capabilities

- Execute multi-step browser tasks from natural-language instructions, including navigation, search, clicking, and pagination.
- Fill inputs, select dropdown options, submit queries, and capture results from interactive web applications.
- Extract tables, lists, and body text from JavaScript-rendered portals that do not expose a public REST API.
- Inspect execution steps and refine instructions when page layouts or form fields change.
- Consolidate multi-column web listings into structured tables for downstream filtering and verification.

### Use Cases

- Query biobanks, reagent catalogs, conference schedules, and funding agency portals without repetitive manual clicking.
- Collect publicly posted supplementary notes, dataset documentation, and protocol pages across publisher sites.
- Monitor calls for papers, grant deadlines, and job postings scattered across institutional websites.
- Run recurring page checks or multi-step form workflows and flag layout changes or broken links for review.

### Usage Notes

- Requires a Browser Use API key configured in the connector header settings.
- Cloud browser sessions consume account credits according to your Browser Use plan; review pricing and concurrency limits before running large batches.
- Homepage: [Browser Use MCP Server](https://docs.browser-use.com/cloud/guides/mcp-server). Respect target websites' terms of service and rate limits when automating access.
