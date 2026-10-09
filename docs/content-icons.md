# Content icons

Content icons travel with their package source and installed files. Add one optional
`icon.svg`, `icon.png`, `icon.webp`, `icon.jpg`, or `icon.jpeg` beside
`phi-package.yaml`, `SKILL.md`, or `wrapper.yaml`. The lookup order is SVG, PNG,
WebP, JPG, then JPEG. Keep the file non-empty and at most 256 KiB. A missing or
unusable icon uses Phi's existing type fallback.

For a generated wrapper family package, place its shared icon at the family root.
The registry builder checks that root first, then its `wrapper/` adapter directory;
it does not promote images from unrelated nested components. An individual wrapper
surface may use the icon beside its own `wrapper.yaml`. Plugin and plugin-component
icons belong to their respective source roots.

The builder includes the source icon in `files.json` and the archive, then copies its
exact bytes to `icons/<type>-<id>-<version>.<extension>` in the registry. The optional
`iconAsset: { path, sha256, size }` index field describes that copy. Its SHA-256 and
size are covered by an index signature when the registry is signed. The local reader
accepts only safe relative image paths and bounded sizes; the app verifies sidecar
bytes before showing them. Installed icons continue to use the package's verified
file allowlist. No manifest field is needed, and registries without `iconAsset`
remain valid. See [package contract](contracts/package.md) § 4.2.

Brand assets are displayed as inert images. The renderer receives an opaque icon
reference rather than a resource filesystem path. SVGs are checked for executable
or external content before display; invalid images retain the type fallback.

## Existing connector attribution

The following 18 icons were copied byte-for-byte on 2026-10-08 from Phi's former
`src/renderer/src/features/mcp/assets/` catalog assets into
`resources/connectors/<id>/icon.<extension>`. Each affected connector's package
version changed from 1.0.0 to 1.0.1 because its payload changed. BioMCP was added
separately from its official website, as recorded below. These provider assets identify services;
this attribution does not change their owners' licensing or trademark terms.

| Icon                  | Listing                                                                                                 |
| --------------------- | ------------------------------------------------------------------------------------------------------- |
| Google Drive          | https://claude.com/marketplace/connectors/google-drive                                                  |
| Gmail                 | https://claude.com/marketplace/connectors/gmail                                                         |
| Notion                | https://claude.com/marketplace/connectors/notion                                                        |
| Composio              | https://logos.composio.dev/api/composio (official logo CDN, retrieved 2026-09-29)                       |
| Tavily                | https://www.tavily.com/brand (official brand mark)                                                      |
| SerpApi               | https://us-west.serpapi.com/media (official logo SVG)                                                   |
| Firecrawl             | https://github.com/firecrawl/firecrawl/blob/main/img/firecrawl_logo.png                                 |
| Browser Use           | https://browser-use.com/logo.png                                                                        |
| Linear                | https://claude.com/marketplace/connectors/linear                                                        |
| Slack                 | https://claude.com/marketplace/connectors/slack                                                         |
| Figma                 | https://claude.com/marketplace/connectors/figma                                                         |
| Canva                 | https://claude.com/marketplace/connectors/canva                                                         |
| BioRender             | https://claude.com/marketplace/connectors/biorender                                                     |
| PubMed                | https://claude.com/marketplace/connectors/pubmed                                                        |
| bioRxiv               | https://claude.com/marketplace/connectors/biorxiv                                                       |
| Clinical Trials       | https://claude.com/marketplace/connectors/clinical-trials                                               |
| cBioPortal            | https://www.cbioportal.org/images/cbioportal_icon.png (official website favicon, retrieved 2026-10-08)  |
| Open Targets Platform | https://opentargets.org/branding (official helix SVG); MCP: https://github.com/opentargets/platform-mcp |

## BioMCP website icon

BioMCP's package-local `resources/connectors/biomcp/icon.png` was copied
byte-for-byte on 2026-10-08 from [its official website logo](https://biomcp.org/assets/icon.png),
which [the website](https://biomcp.org/) references as `assets/icon.png`. The PNG
is 389 × 380 pixels, 90,550 bytes, and has SHA-256
`0e448634f051ccdd5645189482251d948d0a9e4b69b5b5463b7bd05e8afff06c`.
The connector package changed from 1.0.0 to 1.0.1 for this payload addition;
its upstream BioMCP pin stays at 0.9.1. The current 1.1.0 package retains these
exact icon bytes while moving installation into its managed native environment.
This asset is separate from the 18-icon migration and retains its owner's licensing
and trademark terms.

The signed official catalog must be rebuilt from the updated tracked source
and publish the 1.1.0 archive, connector manifest sidecar, icon sidecar,
`index.json`, and matching `index.sig.json` together. Updating this source alone
does not change the published catalog or an installed connector. See
[publishing](publishing.md).

## Migration byte record

These hashes preserve the original image bytes without a frontend asset dependency.

| Package file                                    | SHA-256                                                            | Bytes |
| ----------------------------------------------- | ------------------------------------------------------------------ | ----: |
| `resources/connectors/biorender/icon.jpg`       | `5058279cd72d517919bce9623c803c0d334a918740b2785ac1b47ea87c33b311` |  3566 |
| `resources/connectors/biorxiv/icon.png`         | `ce412025016e5535bbb5aa1d8a3b7e84c67f22dfe2c638d7eb2447188df80d52` | 59601 |
| `resources/connectors/browser-use/icon.png`     | `c0da08801542a6c68f9837530c93820b80da0dcdf086c009866e119f9e13c109` | 43719 |
| `resources/connectors/canva/icon.jpg`           | `ecfbe1a4a05cd342b14b9d7b2323a5637a72a3469b5ce6462a75c37b33af4012` |  2512 |
| `resources/connectors/cbioportal/icon.png`      | `d19e84d8f61376ffec1408a882ed38447a0ba168ddf75a947b2b0c9a0d7522f0` |  2587 |
| `resources/connectors/clinical-trials/icon.png` | `6136f7db817a75fd75a5ef4dfeb548733bfa7d875c2d7e917a2d1485943e0420` | 24196 |
| `resources/connectors/composio/icon.svg`        | `8da32432515f89f3f81992340ae2feda38c6e2f2b807a1a509922a8ca8d1cd53` |  7479 |
| `resources/connectors/figma/icon.jpg`           | `bdece31fd7626aaad1f75c9b96026a867d0d34bebb0ba85cff280f573b0cd604` |  1834 |
| `resources/connectors/firecrawl/icon.png`       | `37e85b43b4b9cdbd689e9e710410caa4c8a240b42dd21ef8edd5edde9a6116df` |  8927 |
| `resources/connectors/gmail/icon.svg`           | `8ae53aecbb00adfecbefb557e10234e44183bd2d826a756f10393e29c0cb2691` |  2418 |
| `resources/connectors/google-drive/icon.svg`    | `373f7efb45157b05017318550023f7e807519ae12e767a40f107326d85187e8b` |  1521 |
| `resources/connectors/linear/icon.svg`          | `a213893f4d76eba22c3bebe05e9627704a3469e0f81f8a732bdf7ebbebf6f4d7` |   499 |
| `resources/connectors/notion/icon.svg`          | `f34ddba7d44c7c8701b73c43dd3762355e460c54f6eb0f0f628e01af4629aeb8` |  1472 |
| `resources/connectors/open-targets/icon.svg`    | `319081e74bcb3d6c2dd2d472e0699cff39488b93713576441c54285ed4c2bb94` |   745 |
| `resources/connectors/pubmed/icon.svg`          | `a677ecb9865f7a0dbfff9863672a05acdb642af1774da954145c83140ddc51cf` |   886 |
| `resources/connectors/serpapi/icon.svg`         | `53656f8526033e4323a0d6c502cf1caf3f4ab17be72c2f83be37bd0cec9163e8` |  1848 |
| `resources/connectors/slack/icon.svg`           | `fd9e7be2ea0056bcbab54e03b027c5be9b1d1c2d157821c06de2d563982c05fb` |  1738 |
| `resources/connectors/tavily/icon.svg`          | `320ef37b5a13331bee4d756c6aed45422e61e73a02816ac3a4768815f6ed5153` |  2282 |
