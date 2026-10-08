import assert from 'node:assert/strict'
import { createHash } from 'node:crypto'
import { existsSync, readFileSync, readdirSync } from 'node:fs'
import { dirname, join } from 'node:path'
import { fileURLToPath } from 'node:url'
import test from 'node:test'

const root = process.env.PHI_PACKAGES_TEST_ROOT ?? fileURLToPath(new URL('..', import.meta.url))
const provenance = JSON.parse(readFileSync(join(root, 'SOURCE.json'), 'utf8'))

test('connector icon migration preserves original bytes, attribution, and package versions', () => {
  const migration = provenance.iconMigration
  assert.equal(migration.assetCount, 18)
  assert.equal(migration.files.length, 18)
  assert.equal(migration.assetBytes, 167830)
  const attribution = readFileSync(join(root, 'docs/content-icons.md'), 'utf8')
  for (const file of migration.files) {
    assert.match(file.path, /^resources\/connectors\/[a-z][a-z0-9-]+\/icon\.(svg|png|webp|jpg|jpeg)$/)
    const bytes = readFileSync(join(root, file.path))
    assert.equal(bytes.length, file.size, file.path)
    assert.equal(createHash('sha256').update(bytes).digest('hex'), file.sha256, file.path)
    assert.ok(attribution.includes(file.path), `attribution lacks ${file.path}`)
    assert.ok(attribution.includes(file.sha256), `byte record lacks ${file.path}`)
    assert.match(readFileSync(join(root, dirname(file.path), 'phi-package.yaml'), 'utf8'), /^version: 1\.0\.1$/m)
  }
})

test('BioMCP keeps its 1.0.0 version and type fallback without an authored logo', () => {
  const connectorRoot = join(root, 'resources/connectors/biomcp')
  assert.match(readFileSync(join(connectorRoot, 'phi-package.yaml'), 'utf8'), /^version: 1\.0\.0$/m)
  assert.equal(readdirSync(connectorRoot).some(name => /^icon\.(svg|png|webp|jpg|jpeg)$/.test(name)), false)
})

test('icon attribution stays outside connector package roots and core resources remain absent', () => {
  const connectorRoot = join(root, 'resources/connectors')
  for (const name of readdirSync(connectorRoot).filter(name => name !== '.DS_Store')) {
    assert.ok(existsSync(join(connectorRoot, name, 'phi-package.yaml')), `${name} is not a connector package`)
  }
  for (const path of ['resources/palettes', 'resources/runtime']) assert.equal(existsSync(join(root, path)), false)
})
