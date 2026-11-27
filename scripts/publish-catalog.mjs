#!/usr/bin/env node
/**
 * Incremental, integrity-checked publisher for Phi's package release assets.
 *
 * Enforces docs/publishing.md invariants:
 *   1. Verifies every archive, iconAsset, and manifestAsset listed in index.json
 *      exists in --assets with exact byte size and SHA-256 before uploading.
 *   2. Treats existing package/icon/manifest assets on the release as immutable
 *      and uploads only new files (skipping already-published filenames).
 *   3. Uploads index.json LAST with --clobber so clients
 *      never observe a partial catalog generation.
 */

import { spawnSync } from 'node:child_process'
import { createHash } from 'node:crypto'
import { existsSync, readFileSync, statSync } from 'node:fs'
import { join, resolve } from 'node:path'

const INDEX_FILES = ['index.json']

function parseArgs(argv) {
  if (argv.includes('--help') || argv.includes('-h')) {
    console.log(
      'Usage: node scripts/publish-catalog.mjs --assets <dir> [--tag catalog-v1] [--repo owner/name] [--batch-size 25] [--dry-run]'
    )
    process.exit(0)
  }
  const get = (flag, fallback) => {
    const idx = argv.indexOf(flag)
    if (idx < 0) return fallback
    const val = argv[idx + 1]
    if (!val || val.startsWith('--')) throw new Error(`${flag} requires a value`)
    return val
  }
  const assetsDir = get('--assets')
  if (!assetsDir) throw new Error('--assets <dir> is required')
  const batchSize = Number.parseInt(get('--batch-size', '25'), 10)
  if (!Number.isInteger(batchSize) || batchSize < 1) {
    throw new Error('--batch-size must be a positive integer')
  }
  return {
    assetsDir: resolve(assetsDir),
    tag: get('--tag', 'catalog-v1'),
    repo: get('--repo', ''),
    batchSize,
    dryRun: argv.includes('--dry-run'),
    existingAssetsJson: get('--existing-assets-json', '')
  }
}

export function sha256File(filePath) {
  const hash = createHash('sha256')
  hash.update(readFileSync(filePath))
  return hash.digest('hex')
}

export function verifyCatalogAssets(assetsDir) {
  for (const name of INDEX_FILES) {
    const p = join(assetsDir, name)
    if (!existsSync(p)) throw new Error(`missing required catalog file: ${name}`)
  }
  const index = JSON.parse(readFileSync(join(assetsDir, 'index.json'), 'utf8'))
  if (!Array.isArray(index.packages)) {
    throw new Error('index.json is missing packages array')
  }

  const referenced = new Map()
  const check = (asset) => {
    const filePath = join(assetsDir, asset.path)
    if (!existsSync(filePath)) {
      throw new Error(`referenced asset missing on disk: ${asset.path}`)
    }
    const stat = statSync(filePath)
    if (stat.size !== asset.size) {
      throw new Error(`size mismatch for ${asset.path}: expected ${asset.size}, got ${stat.size}`)
    }
    const digest = sha256File(filePath)
    if (digest !== asset.sha256) {
      throw new Error(`sha256 mismatch for ${asset.path}: expected ${asset.sha256}, got ${digest}`)
    }
    referenced.set(asset.path, { size: asset.size, sha256: asset.sha256 })
  }

  for (const pkg of index.packages) {
    check({ path: pkg.archive, size: pkg.size, sha256: pkg.sha256 })
    if (pkg.iconAsset) check(pkg.iconAsset)
    if (pkg.manifestAsset) check(pkg.manifestAsset)
  }

  const immutableAssets = [...referenced.keys()].sort()
  return {
    packageCount: index.packages.length,
    immutableAssets,
    assetMeta: referenced
  }
}

export function planUpload(immutableAssets, existingAssets, localMeta = new Map()) {
  const existingMap = new Map()
  for (const item of existingAssets) {
    if (typeof item === 'string') {
      existingMap.set(item, null)
    } else if (item && typeof item === 'object' && item.name) {
      const rawDigest = typeof item.digest === 'string' && item.digest.startsWith('sha256:')
        ? item.digest.slice('sha256:'.length).toLowerCase()
        : (typeof item.sha256 === 'string' ? item.sha256.toLowerCase() : undefined)
      existingMap.set(item.name, {
        size: typeof item.size === 'number' ? item.size : undefined,
        sha256: rawDigest
      })
    }
  }

  const skipped = []
  const toUpload = []
  for (const name of immutableAssets) {
    if (existingMap.has(name)) {
      const remote = existingMap.get(name)
      const local = localMeta instanceof Map ? localMeta.get(name) : undefined
      if (remote && local) {
        if (remote.size !== undefined && remote.size !== local.size) {
          throw new Error(
            `immutable release asset conflict for ${name}: remote size ${remote.size} != local size ${local.size} (bump package version or pass --wrapper-version)`
          )
        }
        if (remote.sha256 && remote.sha256 !== local.sha256) {
          throw new Error(
            `immutable release asset conflict for ${name}: remote sha256 ${remote.sha256} != local sha256 ${local.sha256} (bump package version or pass --wrapper-version)`
          )
        }
      }
      skipped.push(name)
    } else {
      toUpload.push(name)
    }
  }
  return { skipped, toUpload, indexFiles: [...INDEX_FILES] }
}

function runGh(args) {
  const res = spawnSync('gh', args, { encoding: 'utf8' })
  if (res.error) throw res.error
  if (res.status !== 0) {
    throw new Error(`gh ${args.slice(0, 3).join(' ')} failed (${res.status}): ${(res.stderr || res.stdout).trim()}`)
  }
  return res.stdout
}

function fetchExistingAssets(tag, repo) {
  const args = ['release', 'view', tag, '--json', 'assets']
  if (repo) args.push('--repo', repo)
  const out = runGh(args)
  const data = JSON.parse(out)
  return (data.assets || []).map((a) => ({
    name: a.name,
    size: a.size,
    digest: a.digest
  }))
}

function main() {
  const opts = parseArgs(process.argv.slice(2))
  const verified = verifyCatalogAssets(opts.assetsDir)

  let existingAssets = []
  if (opts.existingAssetsJson) {
    existingAssets = JSON.parse(readFileSync(resolve(opts.existingAssetsJson), 'utf8'))
  } else {
    try {
      existingAssets = fetchExistingAssets(opts.tag, opts.repo)
    } catch (err) {
      if (opts.dryRun) {
        console.error(`[publish-catalog] warning: could not query release '${opts.tag}' (${err.message}); assuming 0 existing assets for dry-run`)
        existingAssets = []
      } else {
        throw err
      }
    }
  }

  const plan = planUpload(verified.immutableAssets, existingAssets, verified.assetMeta)
  const summary = {
    tag: opts.tag,
    dryRun: opts.dryRun,
    packagesVerified: verified.packageCount,
    immutableAssetsTotal: verified.immutableAssets.length,
    skippedExisting: plan.skipped.length,
    toUploadNew: plan.toUpload.length,
    indexFilesLast: plan.indexFiles
  }

  if (opts.dryRun) {
    console.log(JSON.stringify(summary, null, 2))
    return
  }

  // 1. Upload all new immutable package/icon/manifest assets in batches first
  for (let i = 0; i < plan.toUpload.length; i += opts.batchSize) {
    const batch = plan.toUpload.slice(i, i + opts.batchSize).map((name) => join(opts.assetsDir, name))
    const args = ['release', 'upload', opts.tag, ...batch]
    if (opts.repo) args.push('--repo', opts.repo)
    console.error(`[publish-catalog] uploading batch ${Math.floor(i / opts.batchSize) + 1} (${batch.length} assets)...`)
    runGh(args)
  }

  // 2. Upload index.json LAST with --clobber
  const indexPaths = plan.indexFiles.map((name) => join(opts.assetsDir, name))
  const finalArgs = ['release', 'upload', opts.tag, ...indexPaths, '--clobber']
  if (opts.repo) finalArgs.push('--repo', opts.repo)
  console.error('[publish-catalog] uploading index.json last...')
  runGh(finalArgs)

  console.log(JSON.stringify({ ...summary, published: true }, null, 2))
}

if (import.meta.url === `file://${process.argv[1]}`) {
  try {
    main()
  } catch (error) {
    console.error(error instanceof Error ? error.message : String(error))
    process.exitCode = 1
  }
}
