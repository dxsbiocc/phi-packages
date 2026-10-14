#!/usr/bin/env node
import { runModule } from "phi-omics-runtime";
import { serverModule } from "./server.js";

runModule(serverModule, "phi-omics-transcriptomics").catch((error: unknown) => {
  process.stderr.write(`${error instanceof Error ? error.message : String(error)}\n`);
  process.exitCode = 1;
});
