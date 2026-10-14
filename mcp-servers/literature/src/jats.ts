import { XMLParser } from "fast-xml-parser";
import { allowsFullText, asRecord, list, normalizeLicense, plain } from "./util.js";

const parser = new XMLParser({
  ignoreAttributes: false,
  attributeNamePrefix: "@_",
  textNodeName: "#text",
  removeNSPrefix: true,
  processEntities: true,
});

export interface TextSection {
  heading: string;
  text: string;
}

export function parseXml(xml: string): unknown {
  return parser.parse(xml.replace(/<!DOCTYPE[\s\S]*?>/i, "")) as unknown;
}

export function licenseFromXml(xml: string): string | null {
  const path = xml.match(/creativecommons\.org\/licenses\/([a-z0-9-]+)/i);
  if (path?.[1]) return normalizeLicense(path[1]);
  if (/creativecommons\.org\/publicdomain\/zero/i.test(xml)) return "cc0";
  return null;
}

function collect(node: unknown, out: TextSection[]): void {
  const record = asRecord(node);
  for (const section of list(record.sec)) {
    const heading = plain(asRecord(section).title).replace(/\s+/g, " ").trim() || "Section";
    const copy = { ...asRecord(section) };
    delete copy.title;
    delete copy.sec;
    const text = plain(copy).replace(/\s+/g, " ").trim();
    if (text) out.push({ heading, text });
    collect(section, out);
  }
}

export function parseJats(xml: string): { sections: TextSection[]; license: string | null } {
  const license = licenseFromXml(xml);
  const root = asRecord(parseXml(xml));
  const article = asRecord(root.article ?? root);
  const meta = asRecord(asRecord(article.front)["article-meta"]);
  const sections: TextSection[] = [];
  const abstract = plain(meta.abstract).replace(/\s+/g, " ").trim();
  if (abstract) sections.push({ heading: "Abstract", text: abstract });
  collect(article.body, sections);
  return { sections: sections.slice(0, 30), license };
}

const GROUP = {
  methods: /method|material/i,
  results: /result/i,
  discussion: /discussion/i,
} as const;

export function clipSections(sections: TextSection[], maxChars: number, wanted: "abstract" | "methods" | "results" | "discussion" | "all"): { sections: TextSection[]; truncated: boolean } {
  let chosen = sections;
  if (wanted === "abstract") chosen = sections.filter((section) => /abstract/i.test(section.heading)).slice(0, 3);
  else if (wanted === "all") chosen = sections.filter((section) => !/reference|bibliography|acknowledgment/i.test(section.heading));
  else {
    const order = ["methods", "results", "discussion"] as const;
    const later = order.slice(order.indexOf(wanted) + 1).map((key) => GROUP[key]);
    const matched: TextSection[] = [];
    let on = false;
    for (const section of sections) {
      if (GROUP[wanted].test(section.heading)) on = true;
      else if (on && later.some((pattern) => pattern.test(section.heading))) break;
      if (on) matched.push(section);
    }
    chosen = matched;
  }
  let remaining = maxChars;
  let truncated = chosen.length === 0 && wanted !== "all";
  const out: TextSection[] = [];
  for (const section of chosen) {
    if (remaining <= 0) {
      truncated = true;
      break;
    }
    const text = section.text.length > remaining ? `${section.text.slice(0, remaining).trimEnd()}…` : section.text;
    if (text.length < section.text.length) truncated = true;
    remaining -= Math.min(section.text.length, remaining);
    if (text.trim()) out.push({ heading: section.heading, text: text.trim() });
  }
  return { sections: out, truncated };
}

export function redistributableLicense(recordLicense: string | null, xmlLicense: string | null): string | null {
  if (allowsFullText(xmlLicense)) return xmlLicense;
  if (allowsFullText(recordLicense)) return recordLicense;
  return xmlLicense ?? recordLicense;
}
