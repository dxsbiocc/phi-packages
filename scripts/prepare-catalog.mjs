import { spawnSync } from 'node:child_process'
import { dirname, resolve } from 'node:path'
import { fileURLToPath } from 'node:url'

const sourceRoot = resolve(dirname(fileURLToPath(import.meta.url)), '..')
const argv = process.argv.slice(2)
if (argv.includes('--help')) {
  console.log('node scripts/prepare-catalog.mjs --phi ../Phi --out /tmp/phi-catalog --key ~/.phi/publishing/phi-packages-ed25519.pem [--wrapper-version 0.1.1]')
  process.exit(0)
}
function argument(name, fallback) {
  const index = argv.indexOf(name)
  if (index < 0 && fallback !== undefined) return fallback
  const value = index < 0 ? undefined : argv[index + 1]
  if (!value || value.startsWith('--')) throw new Error(`${name} requires a value`)
  return value
}
try {
  const phiRoot = resolve(argument('--phi', resolve(sourceRoot, '../Phi')))
  const command = [
    '--import', resolve(phiRoot, 'scripts/test-loader.mjs'),
    resolve(phiRoot, 'scripts/packages/prepare-official-release.ts'),
    '--source', sourceRoot,
    '--out', resolve(argument('--out')),
    '--key', resolve(argument('--key')),
    '--wrapper-version', argument('--wrapper-version', '0.1.1')
  ]
  const result = spawnSync(process.execPath, command, { cwd: phiRoot, stdio: 'inherit' })
  if (result.error) throw result.error
  process.exitCode = result.status ?? 1
} catch (error) {
  console.error(error instanceof Error ? error.message : String(error))
  process.exitCode = 1
}
