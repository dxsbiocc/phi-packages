import { AsyncLocalStorage } from "node:async_hooks";

/**
 * What one tool call learned while it ran. The upstream client reports when each response was
 * really fetched (a cached response is older than the call), and row mappers report rows they
 * had to skip. The MCP wrapper turns both into the result's provenance.
 */
export interface CallContext {
  fetchedAt: number[];
  skipped: number;
}

const storage = new AsyncLocalStorage<CallContext>();

export async function withCallContext<T>(fn: () => Promise<T>): Promise<{ value: T; context: CallContext }> {
  const context: CallContext = { fetchedAt: [], skipped: 0 };
  const value = await storage.run(context, fn);
  return { value, context };
}

/** Record that an upstream response fetched at `timestampMs` fed the current call. */
export function noteFetch(timestampMs: number): void {
  storage.getStore()?.fetchedAt.push(timestampMs);
}

/** Record rows dropped because they did not fit the shared schema. */
export function noteSkipped(count: number): void {
  const context = storage.getStore();
  if (context) context.skipped += count;
}

/** The oldest upstream fetch behind this call, or null when nothing was fetched upstream. */
export function oldestFetch(context: CallContext): number | null {
  return context.fetchedAt.length === 0 ? null : Math.min(...context.fetchedAt);
}
