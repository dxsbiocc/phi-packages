import { z } from "zod";
import { EntitySchema, EntityTypeSchema } from "./entities.js";

export const PageInfoSchema = z.object({
  offset: z.number().int().nonnegative(),
  limit: z.number().int().positive(),
  returned: z.number().int().nonnegative(),
  total: z.number().int().nonnegative().nullable(),
  next_offset: z.number().int().nonnegative().nullable(),
});

export const SourceSchema = z.object({
  id: z.string().min(1),
  name: z.string().min(1),
  version: z.string().nullable(),
  homepage: z.string().min(1),
  retrieved_at: z.string().min(1),
});

export const UsageSchema = z.object({
  summary: z.string(),
  url: z.string().nullable(),
});

export const ErrorCodeSchema = z.enum([
  "upstream_http",
  "not_found",
  "invalid_input",
  "unavailable",
]);

export const PageSchema = z.object({
  entity_type: z.literal("page"),
  item_type: EntityTypeSchema,
  items: z.array(EntitySchema),
  page: PageInfoSchema,
  source: SourceSchema,
  usage: UsageSchema,
  source_extra: z.record(z.string(), z.unknown()),
}).superRefine((page, ctx) => {
  if (page.items.length !== page.page.returned) {
    ctx.addIssue({ code: "custom", message: "returned does not match items", path: ["page", "returned"] });
  }
  page.items.forEach((item, index) => {
    if (item.entity_type !== page.item_type) {
      ctx.addIssue({
        code: "custom",
        message: `item entity_type must be ${page.item_type}`,
        path: ["items", index, "entity_type"],
      });
    }
  });
  const { offset, limit, returned, total, next_offset } = page.page;
  if (returned > limit) {
    ctx.addIssue({ code: "custom", message: "returned exceeds limit", path: ["page", "returned"] });
  }
  if (total != null && offset + returned < total && next_offset == null) {
    ctx.addIssue({ code: "custom", message: "next_offset is missing", path: ["page", "next_offset"] });
  }
  if (total != null && offset + returned >= total && next_offset != null) {
    ctx.addIssue({ code: "custom", message: "next_offset must be null on the last page", path: ["page", "next_offset"] });
  }
});

export const ErrorObjectSchema = z.object({
  entity_type: z.literal("error"),
  code: ErrorCodeSchema,
  message: z.string().min(1),
  source: SourceSchema,
  retryable: z.boolean(),
  source_extra: z.record(z.string(), z.unknown()),
});

export const ToolResultSchema = z.discriminatedUnion("entity_type", [PageSchema, ErrorObjectSchema]);
