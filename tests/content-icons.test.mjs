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

test('BioMCP packages its official website icon with provenance and an updated version', () => {
  const connector = provenance.addedContent.find(item => item.path === 'resources/connectors/biomcp/')
  assert.equal(connector.packageVersion, '1.1.0')
  const icon = connector.icon
  assert.equal(icon.path, 'resources/connectors/biomcp/icon.png')
  assert.equal(icon.sourceUrl, 'https://biomcp.org/assets/icon.png')
  const bytes = readFileSync(join(root, icon.path))
  assert.equal(bytes.length, icon.size)
  assert.ok(bytes.length > 0 && bytes.length <= 256 * 1024)
  assert.deepEqual(bytes.subarray(0, 8), Buffer.from([137, 80, 78, 71, 13, 10, 26, 10]))
  assert.equal(createHash('sha256').update(bytes).digest('hex'), icon.sha256)
  const attribution = readFileSync(join(root, 'docs/content-icons.md'), 'utf8')
  for (const value of [icon.path, icon.sourceUrl, icon.sha256]) assert.ok(attribution.includes(value))
  assert.match(readFileSync(join(root, 'resources/connectors/biomcp/phi-package.yaml'), 'utf8'), /^version: 1\.1\.0$/m)
})

test('icon attribution stays outside connector package roots and core resources remain absent', () => {
  const connectorRoot = join(root, 'resources/connectors')
  for (const name of readdirSync(connectorRoot).filter(name => name !== '.DS_Store')) {
    assert.ok(existsSync(join(connectorRoot, name, 'phi-package.yaml')), `${name} is not a connector package`)
  }
  for (const path of ['resources/palettes', 'resources/runtime']) assert.equal(existsSync(join(root, path)), false)
})
