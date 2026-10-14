export {
  CURIE_PATTERN,
  CurieSchema,
  curie,
  ensemblGeneCurie,
  humanChrom,
  ontologyCurie,
  sexOf,
  strandOf,
  taxonCurie,
} from "./ids.js";
export {
  ENTITY_TYPES,
  AccessSchema,
  CellTypeSchema,
  DatasetSchema,
  DonorSchema,
  EntitySchema,
  EntityTypeSchema,
  ExpressionSchema,
  ExpressionUnitSchema,
  FileSchema,
  GeneSchema,
  LocusSchema,
  MetaboliteSchema,
  ProteinSchema,
  SampleSchema,
  SexSchema,
  SourceExtraSchema,
  StrandSchema,
  TissueSchema,
  VariantSchema,
} from "./entities.js";
export {
  ErrorCodeSchema,
  ErrorObjectSchema,
  PageInfoSchema,
  PageSchema,
  SourceSchema,
  ToolResultSchema,
  UsageSchema,
} from "./envelope.js";

import type { z } from "zod";
import {
  CellTypeSchema,
  DatasetSchema,
  DonorSchema,
  ExpressionSchema,
  FileSchema,
  GeneSchema,
  LocusSchema,
  MetaboliteSchema,
  ProteinSchema,
  SampleSchema,
  TissueSchema,
  VariantSchema,
} from "./entities.js";
import { ErrorObjectSchema, PageSchema, SourceSchema, ToolResultSchema, UsageSchema } from "./envelope.js";

export type Gene = z.infer<typeof GeneSchema>;
export type Sample = z.infer<typeof SampleSchema>;
export type Donor = z.infer<typeof DonorSchema>;
export type Dataset = z.infer<typeof DatasetSchema>;
export type OmicsFile = z.infer<typeof FileSchema>;
export type CellType = z.infer<typeof CellTypeSchema>;
export type Tissue = z.infer<typeof TissueSchema>;
export type Protein = z.infer<typeof ProteinSchema>;
export type Metabolite = z.infer<typeof MetaboliteSchema>;
export type Expression = z.infer<typeof ExpressionSchema>;
export type Variant = z.infer<typeof VariantSchema>;
export type Entity = Gene | Sample | Donor | Dataset | OmicsFile | CellType | Tissue | Protein | Metabolite | Expression | Variant;
export type Locus = z.infer<typeof LocusSchema>;
export type Page = z.infer<typeof PageSchema>;
export type ErrorObject = z.infer<typeof ErrorObjectSchema>;
export type ToolResult = z.infer<typeof ToolResultSchema>;
export type Source = z.infer<typeof SourceSchema>;
export type Usage = z.infer<typeof UsageSchema>;
export type EntityType = Entity["entity_type"];
export type Access = "open" | "controlled" | "unknown";
export type ExpressionUnit = "TPM" | "log2FC" | "count" | "UMI" | "protein_level" | "z_score";
