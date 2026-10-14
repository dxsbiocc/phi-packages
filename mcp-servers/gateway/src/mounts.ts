import type { HttpMount } from "phi-omics-runtime";
import { serverModule as cancer } from "phi-omics-cancer";
import { serverModule as epigenomics } from "phi-omics-epigenomics";
import { serverModule as metabolomics } from "phi-omics-metabolomics";
import { serverModule as proteomics } from "phi-omics-proteomics";
import { serverModule as transcriptomics } from "phi-omics-transcriptomics";

export const gatewayMounts: HttpMount[] = [
  { path: "/transcriptomics/mcp", module: transcriptomics },
  { path: "/cancer/mcp", module: cancer },
  { path: "/epigenomics/mcp", module: epigenomics },
  { path: "/proteomics/mcp", module: proteomics },
  { path: "/metabolomics/mcp", module: metabolomics },
];
