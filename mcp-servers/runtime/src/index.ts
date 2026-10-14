export { cacheDirFor, loadConfig, type RuntimeConfig } from "./config.js";
export { InputError, UpstreamClient, UpstreamError, type JsonClient, type UpstreamRequest } from "./http.js";
export {
  asArray,
  asRecord,
  cellType,
  compound,
  dataset,
  drugResponse,
  donor,
  expression,
  failure,
  guard,
  file,
  gene,
  isRecord,
  labels,
  locusOf,
  makePage,
  metabolite,
  numberOrNull,
  pageInfo,
  paging,
  protein,
  publication,
  sample,
  sliceList,
  sourceOf,
  text,
  tissue,
  usagePage,
  variant,
} from "./page.js";
export {
  createServerFromTools,
  helpText,
  readListenOptions,
  startHttp,
  startStdio,
  toolResult,
  type HttpMount,
  type OmicsTool,
  type ServerModule,
} from "./mcp.js";
import { loadConfig } from "./config.js";
import { UpstreamClient } from "./http.js";
import { helpText, readListenOptions, startHttp, startStdio, type ServerModule } from "./mcp.js";

export async function runModule(mod: ServerModule, bin: string, argv: string[] = process.argv.slice(2)): Promise<void> {
  const options = readListenOptions(argv);
  if (options.help) {
    process.stdout.write(helpText(bin, `${mod.name} ${mod.version}`));
    return;
  }
  const client = new UpstreamClient(loadConfig());
  if (!options.http) {
    startStdio(mod, client);
    return;
  }
  const listener = await startHttp([{ path: mod.mcpPath, module: mod }], options, client);
  process.stderr.write(`${mod.name} listening on http://${options.host}:${listener.port}${options.basePath}/health and ${mod.mcpPath}\n`);
}
