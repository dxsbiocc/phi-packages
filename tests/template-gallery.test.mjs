import assert from 'node:assert/strict'
import { spawnSync } from 'node:child_process'
import { existsSync, mkdirSync, mkdtempSync, readdirSync, readFileSync, rmSync, statSync, symlinkSync, writeFileSync } from 'node:fs'
import { tmpdir } from 'node:os'
import { join, sep } from 'node:path'
import test from 'node:test'
import { fileURLToPath } from 'node:url'

const root = process.env.PHI_PACKAGES_TEST_ROOT ?? fileURLToPath(new URL('..', import.meta.url))
const skillRoot = join(root, 'resources/plugins/visualization/skills/omics-visualization')
const VARIANT_PREVIEW = /^scripts\/[a-z0-9_]+\/[a-z0-9_]+\/preview_[A-Za-z0-9][A-Za-z0-9_-]{0,63}\.png$/
const CATALOG_PREVIEW = /^scripts\/[a-z0-9_]+\/[a-z0-9_]+\/preview\.png$/

test('template gallery matches the omics-visualization catalogs', () => {
  const check = spawnSync('python3', ['scripts/build_template_gallery.py', '--check'], {
    cwd: skillRoot,
    encoding: 'utf8'
  })
  assert.equal(check.status, 0, `${check.stdout}\n${check.stderr}`)

  const gallery = JSON.parse(readFileSync(join(skillRoot, 'references/template_gallery.json'), 'utf8'))
  assert.equal(gallery.schemaVersion, 1)
  assert.equal(gallery.skill, 'omics-visualization')
  assert.equal(gallery.previewBackground, 'transparent')
  assert.ok(Array.isArray(gallery.families))

  const templates = gallery.families.flatMap((family) => family.templates)
  // 211 catalog templates, each with its primary preview.png.
  const primaries = templates.filter((template) => template.variant === undefined)
  assert.equal(primaries.length, 211)
  assert.equal(new Set(primaries.map((template) => template.id)).size, primaries.length)
  // A variant repeats the template id, so uniqueness is (id, variant).
  const keys = templates.map((template) => `${template.id}\0${template.variant ?? ''}`)
  assert.equal(new Set(keys).size, keys.length)
  for (const template of templates) {
    assert.match(template.preview, template.variant ? VARIANT_PREVIEW : CATALOG_PREVIEW)
    assert.equal(existsSync(join(skillRoot, template.preview)), true, template.preview)
  }

  // Every shipped preview_<name>.png (outside the layout audit set) is listed.
  const onDisk = []
  function walk(directory) {
    for (const name of readdirSync(directory)) {
      const path = join(directory, name)
      if (statSync(path).isDirectory()) walk(path)
      else if (/^preview_[A-Za-z0-9][A-Za-z0-9_-]{0,63}\.png$/.test(name)) onDisk.push(path)
    }
  }
  walk(join(skillRoot, 'scripts'))
  const listed = new Set(templates.filter((template) => template.variant).map((template) => template.preview))
  for (const path of onDisk) {
    if (path.includes(`${sep}layouts${sep}`)) continue
    const relative = path.slice(skillRoot.length + 1).split(sep).join('/')
    assert.equal(listed.has(relative), true, relative)
  }
  assert.equal(listed.size, onDisk.filter((path) => !path.includes(`${sep}layouts${sep}`)).length)

  const manifest = readFileSync(join(root, 'resources/plugins/visualization/phi-package.yaml'), 'utf8')
  assert.match(manifest, /^version: 1\.2\.3$/m)
})

test('gallery lists preview_<name>.png drawings and skips escaped or unnamed files', () => {
  const fixture = mkdtempSync(join(tmpdir(), 'phi-gallery-variants-'))
  const skill = join(fixture, 'skill')
  const basic = join(skill, 'scripts', 'scatter', 'basic')
  mkdirSync(join(skill, 'references', 'catalog'), { recursive: true })
  mkdirSync(basic, { recursive: true })
  mkdirSync(join(fixture, 'outside'), { recursive: true })
  for (const name of ['preview.png', 'preview_zeta.png', 'preview_with_trend.png', 'preview_.png', 'notes.png']) {
    writeFileSync(join(basic, name), 'png')
  }
  writeFileSync(join(fixture, 'outside', 'preview.png'), 'png')
  symlinkSync(join(fixture, 'outside', 'preview.png'), join(basic, 'preview_leaked.png'))
  writeFileSync(
    join(skill, 'references', 'plots.yaml'),
    `groups:
  scatter:
    title: "Scatter"
    catalog: "references/catalog/scatter.yaml"
`
  )
  writeFileSync(
    join(skill, 'references', 'catalog', 'scatter.yaml'),
    `templates:
  - id: scatter-basic
    title: "Basic scatter plot"
    preview: "scripts/scatter/basic/preview.png"
`
  )
  const script = `
import json, sys
sys.path.insert(0, ${JSON.stringify(join(skillRoot, 'scripts'))})
from build_template_gallery import gallery_document
document, errors = gallery_document(${JSON.stringify(skill)})
json.dump({"document": document, "errors": errors}, sys.stdout)
`
  try {
    const result = spawnSync('python3', ['-c', script], { encoding: 'utf8' })
    assert.equal(result.status, 0, result.stderr)
    const parsed = JSON.parse(result.stdout)
    assert.deepEqual(parsed.errors, [])
    assert.deepEqual(parsed.document.families[0].templates, [
      {
        id: 'scatter-basic',
        title: 'Basic scatter plot',
        preview: 'scripts/scatter/basic/preview.png'
      },
      {
        id: 'scatter-basic',
        title: 'Basic scatter plot · with trend',
        preview: 'scripts/scatter/basic/preview_with_trend.png',
        variant: 'with_trend'
      },
      {
        id: 'scatter-basic',
        title: 'Basic scatter plot · zeta',
        preview: 'scripts/scatter/basic/preview_zeta.png',
        variant: 'zeta'
      }
    ])
  } finally {
    rmSync(fixture, { recursive: true, force: true })
  }
})
