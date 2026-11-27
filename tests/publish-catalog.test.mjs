import assert from 'node:assert/strict'
import { spawnSync } from 'node:child_process'
import { createHash } from 'node:crypto'
import { mkdtempSync, rmSync, writeFileSync } from 'node:fs'
import { tmpdir } from 'node:os'
import { join, resolve } from 'node:path'
import test from 'node:test'
import { fileURLToPath } from 'node:url'

import { planUpload, verifyCatalogAssets } from '../scripts/publish-catalog.mjs'

const repoRoot = resolve(fileURLToPath(new URL('..', import.meta.url)))
const scriptPath = join(repoRoot, 'scripts/publish-catalog.mjs')

function sha256(content) {
  return createHash('sha256').update(content).digest('hex')
}

function makeFixtureCatalog() {
  const dir = mkdtempSync(join(tmpdir(), 'phi-publish-test-'))
  const archiveA = Buffer.from('archive-a-bytes')
  const iconA = Buffer.from('<svg/>')
  const archiveB = Buffer.from('archive-b-bytes')

  writeFileSync(join(dir, 'plugin-a-1.0.0.tar.gz'), archiveA)
  writeFileSync(join(dir, 'plugin-a-1.0.0.svg'), iconA)
  writeFileSync(join(dir, 'plugin-b-2.0.0.tar.gz'), archiveB)

  const index = {
    packages: [
      {
        id: 'plugin-a',
        type: 'plugin',
        version: '1.0.0',
        archive: 'plugin-a-1.0.0.tar.gz',
        size: archiveA.length,
        sha256: sha256(archiveA),
        iconAsset: { path: 'plugin-a-1.0.0.svg', size: iconA.length, sha256: sha256(iconA) }
      },
      {
        id: 'plugin-b',
        type: 'plugin',
        version: '2.0.0',
        archive: 'plugin-b-2.0.0.tar.gz',
        size: archiveB.length,
        sha256: sha256(archiveB)
      }
    ]
  }
  writeFileSync(join(dir, 'index.json'), JSON.stringify(index, null, 2))
  return dir
}

test('publish-catalog verifies SHA-256/sizes and plans incremental uploads with index last', () => {
  const dir = makeFixtureCatalog()
  try {
    const verified = verifyCatalogAssets(dir)
    assert.equal(verified.packageCount, 2)
    assert.deepEqual(verified.immutableAssets, [
      'plugin-a-1.0.0.svg',
      'plugin-a-1.0.0.tar.gz',
      'plugin-b-2.0.0.tar.gz'
    ])

    const plan = planUpload(verified.immutableAssets, ['plugin-a-1.0.0.tar.gz', 'plugin-a-1.0.0.svg'])
    assert.deepEqual(plan.skipped, ['plugin-a-1.0.0.svg', 'plugin-a-1.0.0.tar.gz'])
    assert.deepEqual(plan.toUpload, ['plugin-b-2.0.0.tar.gz'])
    assert.deepEqual(plan.indexFiles, ['index.json'])

    const existingPath = join(dir, 'existing.json')
    writeFileSync(existingPath, JSON.stringify(['plugin-a-1.0.0.tar.gz']))
    const res = spawnSync(
      process.execPath,
      [scriptPath, '--assets', dir, '--existing-assets-json', existingPath, '--dry-run'],
      { encoding: 'utf8' }
    )
    assert.equal(res.status, 0, res.stderr)
    const out = JSON.parse(res.stdout)
    assert.equal(out.packagesVerified, 2)
    assert.equal(out.skippedExisting, 1)
    assert.equal(out.toUploadNew, 2)
    assert.deepEqual(out.indexFilesLast, ['index.json'])
  } finally {
    rmSync(dir, { recursive: true, force: true })
  }
})

test('publish-catalog rejects corrupted assets before any upload', () => {
  const dir = makeFixtureCatalog()
  try {
    writeFileSync(join(dir, 'plugin-b-2.0.0.tar.gz'), Buffer.from('corrupted-bytes'))
    assert.throws(() => verifyCatalogAssets(dir), /sha256 mismatch/)
  } finally {
    rmSync(dir, { recursive: true, force: true })
  }
})

test('publish-catalog detects unbumped immutable release asset digest/size conflicts', () => {
  const dir = makeFixtureCatalog()
  try {
    const verified = verifyCatalogAssets(dir)
    assert.throws(
      () =>
        planUpload(
          verified.immutableAssets,
          [{ name: 'plugin-a-1.0.0.tar.gz', size: 999999, digest: 'sha256:0000' }],
          verified.assetMeta
        ),
      /immutable release asset conflict/
    )
  } finally {
    rmSync(dir, { recursive: true, force: true })
  }
})
