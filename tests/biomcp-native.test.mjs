import assert from 'node:assert/strict'
import { execFileSync } from 'node:child_process'
import { createHash } from 'node:crypto'
import {
  cpSync,
  existsSync,
  mkdirSync,
  mkdtempSync,
  readFileSync,
  readdirSync,
  rmSync,
  statSync
} from 'node:fs'
import { tmpdir } from 'node:os'
import { basename, join, resolve } from 'node:path'
import { fileURLToPath, pathToFileURL } from 'node:url'
import test from 'node:test'

// Retired launcher coverage now belongs to Phi's envs-application-artifacts,
// envs-application-native, and envs-ensure suites: digest/size verification,
// atomic cache publication, concurrency, offline reuse, unsafe ZIP members,
// cancellation, and runtime readiness. This suite verifies the actual package
// declarations, original upstream pins, ownership, and distributable payload.
// Run with Phi's existing TypeScript loader; no test dependency is vendored here.
const root = process.env.PHI_PACKAGES_TEST_ROOT ?? fileURLToPath(new URL('..', import.meta.url))
const phiRoot = process.env.PHI_SOURCE_ROOT ?? resolve(root, '../Phi')
const phiModule = (relative) => import(pathToFileURL(join(phiRoot, relative)).href)
const { parsePackageManifestText } = await phiModule('src/main/agent/packages/manifest.ts')
const { computeEnvId, parseEnvironmentSpec, parseExplicitLock } = await phiModule(
  'src/main/agent/envs/contract.ts'
)
const { describeMcpPackageEnvironment } = await phiModule(
  'src/main/agent/mcp/package-environment.ts'
)
const { buildRegistry } = await phiModule('scripts/packages/build-registry.ts')
const { parseTarGz } = await phiModule('src/main/agent/packages/archive.ts')
const { installNativeApplication } = await phiModule('src/main/agent/envs/applications/native.ts')
const connectorRoot = join(root, 'resources/connectors/biomcp')
const platforms = ['darwin-arm64', 'darwin-x64', 'linux-x64']
const originalPins = {
  'darwin-arm64': {
    filename: 'biomcp_cli-0.9.1-py3-none-macosx_11_0_arm64.whl',
    sha256: '934befced41c457e3eea72a33d2c06d88fdbdcbb7780a4614f216db5ac9da859',
    size: 15340685
  },
  'darwin-x64': {
    filename: 'biomcp_cli-0.9.1-py3-none-macosx_10_12_x86_64.whl',
    sha256: '36a94801839a76bf5236263a56145a911b1ba7c79d3a7fed153cc337c6f4ea85',
    size: 16534235
  },
  'linux-x64': {
    filename: 'biomcp_cli-0.9.1-py3-none-manylinux_2_28_x86_64.whl',
    sha256: '83bad39c58aca37afc1aca9a2099f47447c0bd07e3464fab12465613715d8864',
    size: 16802942
  }
}

function manifest() {
  return parsePackageManifestText(readFileSync(join(connectorRoot, 'phi-package.yaml'), 'utf8'))
}

function environment() {
  const parsed = parseEnvironmentSpec(readFileSync(join(connectorRoot, 'environment.yml'), 'utf8'))
  assert.equal(parsed.ok, true, parsed.errors?.join('; '))
  return parsed.spec
}

function id(descriptor) {
  return computeEnvId({
    scope: descriptor.scope,
    owner: descriptor.owner,
    name: descriptor.spec.name,
    platform: descriptor.platform,
    lockText: descriptor.lockText,
    sourcePackages: descriptor.spec.sourcePackages,
    installation: descriptor.spec.installation
  })
}

async function withTemporaryRoot(body) {
  const directory = mkdtempSync(join(tmpdir(), 'phi-biomcp-native-test-'))
  try {
    return await body(directory)
  } finally {
    rmSync(directory, { recursive: true, force: true })
  }
}

test('BioMCP launches the managed native executable directly and gates unsupported app versions', () => {
  const parsed = manifest()
  assert.equal(parsed.id, 'biomcp')
  assert.equal(parsed.version, '1.1.0')
  assert.equal(parsed.minAppVersion, '1.0.1')
  assert.equal(parsed.connector.transport, 'stdio')
  assert.equal(parsed.connector.environment, './environment.yml')
  assert.equal(parsed.connector.command, 'biomcp')
  assert.deepEqual(parsed.connector.args, ['serve'])
})

test('BioMCP declares a native installation without Conda or host-runtime dependencies', () => {
  const spec = environment()
  assert.equal(spec.name, 'biomcp')
  assert.deepEqual(spec.channels, [])
  assert.deepEqual(spec.dependencies, [])
  assert.equal(spec.host, undefined)
  assert.equal(spec.sourcePackages, undefined)
  assert.equal(spec.installation.backend, 'native')
  assert.equal(spec.installation.executable, 'biomcp')
  assert.deepEqual(Object.keys(spec.installation.artifacts).sort(), platforms.toSorted())
})

for (const platform of platforms) {
  test(`${platform}: original official BioMCP 0.9.1 wheel and exact native member remain pinned`, () => {
    const artifact = environment().installation.artifacts[platform]
    const pin = originalPins[platform]
    const url = new URL(artifact.url)
    assert.equal(url.protocol, 'https:')
    assert.equal(url.hostname, 'files.pythonhosted.org')
    assert.equal(url.username, '')
    assert.equal(url.password, '')
    assert.equal(basename(url.pathname), pin.filename)
    assert.equal(artifact.sha256, pin.sha256)
    assert.equal(artifact.size, pin.size)
    assert.equal(artifact.format, 'zip')
    assert.equal(artifact.member, 'biomcp_cli-0.9.1.data/scripts/biomcp')
  })

  test(`${platform}: explicit lock is valid and contains no Conda packages`, () => {
    const text = readFileSync(join(connectorRoot, 'locks', `${platform}.txt`), 'utf8')
    assert.equal(text, '@EXPLICIT\n')
    const lock = parseExplicitLock(text, { allowEmpty: true })
    assert.equal(lock.ok, true, lock.errors?.join('; '))
    assert.deepEqual(lock.entries, [])
  })

  test(`${platform}: real MCP environment lookup owns and hashes the native installation`, () => {
    const record = describeMcpPackageEnvironment(manifest(), connectorRoot, { platform })
    assert.ok(record)
    assert.equal(record.kind, 'package')
    assert.equal(record.descriptor.scope, 'mcp')
    assert.equal(record.descriptor.owner, 'biomcp')
    assert.equal(record.descriptor.platform, platform)
    assert.equal(record.descriptor.spec.installation.backend, 'native')
    assert.equal(record.envId, id(record.descriptor))
    assert.match(record.envId, /^mcp-biomcp-biomcp-[a-f0-9]{12}$/)
  })
}

test('native artifact changes cannot reuse the prior environment identity', () => {
  const record = describeMcpPackageEnvironment(manifest(), connectorRoot, {
    platform: 'darwin-arm64'
  })
  const changed = structuredClone(record.descriptor)
  changed.spec.installation.artifacts['darwin-arm64'].sha256 = '0'.repeat(64)
  assert.notEqual(id(changed), record.envId)
})

test('BioMCP has no package-owned download launcher or second artifact-pin document', () => {
  for (const name of ['server.py', 'upstream.json', '__pycache__'])
    assert.equal(existsSync(join(connectorRoot, name)), false)
  assert.deepEqual(
    readdirSync(join(connectorRoot, 'locks')).sort(),
    platforms.map((platform) => `${platform}.txt`).sort()
  )
})

test('upstream provenance and the MIT notice remain with the native package', () => {
  const provenance = JSON.parse(readFileSync(join(root, 'SOURCE.json'), 'utf8'))
  const connector = provenance.addedContent.find(
    (item) => item.path === 'resources/connectors/biomcp/'
  )
  assert.equal(connector.packageVersion, '1.1.0')
  assert.equal(connector.upstream.repository, 'https://github.com/genomoncology/biomcp')
  assert.equal(connector.upstream.package, 'biomcp-cli')
  assert.equal(connector.upstream.version, '0.9.1')
  assert.equal(connector.upstream.index, 'https://pypi.org/project/biomcp-cli/0.9.1/')
  assert.equal(connector.upstream.pins, 'environment.yml')
  const license = readFileSync(join(connectorRoot, 'LICENSE.upstream'), 'utf8')
  assert.match(license, /MIT License/)
  assert.match(license, /Copyright \(c\) 2025 Ian Maurer/)
})

test('the real registry builder includes native pins, every lock, icon, and license in the verified payload', async () => {
  await withTemporaryRoot(async (directory) => {
    const source = join(directory, 'source')
    for (const name of ['connectors', 'plugins', 'skills', 'wrappers'])
      mkdirSync(join(source, 'resources', name), { recursive: true })
    cpSync(connectorRoot, join(source, 'resources/connectors/biomcp'), { recursive: true })
    execFileSync('git', ['init', '-q', source])
    execFileSync('git', ['-C', source, 'add', 'resources'])
    const output = join(directory, 'registry')
    const index = buildRegistry({ repoRoot: source, outDir: output })
    assert.equal(index.packages.length, 1)
    const entry = index.packages[0]
    assert.equal(entry.version, '1.1.0')
    assert.equal(entry.minAppVersion, '1.0.1')
    const archived = parseTarGz(readFileSync(join(output, entry.archive)))
    const files = new Map(archived.map((file) => [file.path, file.data]))
    const allowlist = JSON.parse(files.get('files.json').toString('utf8')).files
    for (const name of [
      'environment.yml',
      'icon.png',
      'LICENSE.upstream',
      ...platforms.map((platform) => `locks/${platform}.txt`)
    ]) {
      const bytes = readFileSync(join(connectorRoot, name))
      assert.deepEqual(files.get(name), bytes)
      const listed = allowlist.find((file) => file.path === name)
      assert.equal(listed.size, bytes.length)
      assert.equal(listed.sha256, createHash('sha256').update(bytes).digest('hex'))
    }
    for (const name of ['server.py', 'upstream.json']) assert.equal(files.has(name), false)
    assert.deepEqual(readFileSync(join(output, entry.iconAsset.path)), files.get('icon.png'))
    assert.deepEqual(
      readFileSync(join(output, entry.manifestAsset.path)),
      files.get('phi-package.yaml')
    )
  })
})

// Optional offline integration against official wheels independently downloaded
// and verified into a disposable directory; regular content tests never fetch.
for (const platform of platforms) {
  test(
    `${platform}: Phi native installer extracts the declared official wheel and reuses its verified artifact offline`,
    { skip: !process.env.PHI_BIOMCP_WHEEL_DIR },
    async () => {
      await withTemporaryRoot(async (directory) => {
        const spec = environment()
        const artifact = spec.installation.artifacts[platform]
        const bytes = readFileSync(
          join(process.env.PHI_BIOMCP_WHEEL_DIR, originalPins[platform].filename)
        )
        assert.equal(bytes.length, artifact.size)
        assert.equal(createHash('sha256').update(bytes).digest('hex'), artifact.sha256)
        const prefix = join(directory, 'prefix')
        const input = {
          root: directory,
          prefix,
          sourceDir: connectorRoot,
          platform,
          installation: spec.installation
        }
        const result = await installNativeApplication({
          ...input,
          fetch: async (url) => {
            assert.equal(String(url), artifact.url)
            return new Response(bytes)
          }
        })
        const binary = join(prefix, 'bin/biomcp')
        const installed = readFileSync(binary)
        assert.ok(installed.length > 0)
        assert.equal(
          installed.subarray(0, 4).toString('hex'),
          platform === 'linux-x64' ? '7f454c46' : 'cffaedfe'
        )
        assert.equal(statSync(binary).mode & 0o777, 0o700)
        assert.deepEqual(readdirSync(prefix), ['bin'])
        assert.equal(result.backend, 'native')
        assert.equal(result.artifacts[0].sha256, artifact.sha256)
        await installNativeApplication({
          ...input,
          fetch: async () => {
            throw new Error('must stay offline')
          }
        })
        assert.deepEqual(readFileSync(binary), installed)
      })
    }
  )
}
