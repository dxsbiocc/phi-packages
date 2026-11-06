import assert from 'node:assert/strict'
import { spawnSync } from 'node:child_process'
import { existsSync, readFileSync } from 'node:fs'
import { join } from 'node:path'
import test from 'node:test'
import { fileURLToPath } from 'node:url'

const root = process.env.PHI_PACKAGES_TEST_ROOT ?? fileURLToPath(new URL('..', import.meta.url))
const skillRoot = join(root, 'resources/plugins/visualization/skills/omics-visualization')

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
  assert.equal(templates.length, 211)
  const ids = new Set(templates.map((template) => template.id))
  assert.equal(ids.size, templates.length)
  for (const template of templates) {
    assert.match(template.preview, /^scripts\/[a-z0-9_]+\/[a-z0-9_]+\/preview\.png$/)
    assert.equal(existsSync(join(skillRoot, template.preview)), true, template.preview)
  }

  const manifest = readFileSync(join(root, 'resources/plugins/visualization/phi-package.yaml'), 'utf8')
  assert.match(manifest, /^version: 1\.2\.0$/m)
})
