import assert from 'node:assert/strict'
import { existsSync, readdirSync, readFileSync, statSync } from 'node:fs'
import { join } from 'node:path'
import { fileURLToPath } from 'node:url'
import test from 'node:test'

const root = process.env.PHI_PACKAGES_TEST_ROOT ?? fileURLToPath(new URL('..', import.meta.url))

function listSubdirs(dir) {
  return readdirSync(dir)
    .filter(name => name !== '.DS_Store' && statSync(join(dir, name)).isDirectory())
    .sort()
}

function countWrapperAdapters(dir) {
  let total = 0
  const stack = [dir]
  while (stack.length > 0) {
    const current = stack.pop()
    for (const entry of readdirSync(current, { withFileTypes: true })) {
      if (entry.name === '.DS_Store') continue
      const full = join(current, entry.name)
      if (entry.isDirectory()) stack.push(full)
      else if (entry.name === 'wrapper.yaml' && current.endsWith('/wrapper')) total += 1
    }
  }
  return total
}

test('READMEs and SOURCE.json counts stay synchronized with repository directories', () => {
  const provenance = JSON.parse(readFileSync(join(root, 'SOURCE.json'), 'utf8'))
  const readmeZh = readFileSync(join(root, 'README.md'), 'utf8')
  const readmeEn = readFileSync(join(root, 'README.en.md'), 'utf8')

  const skills = listSubdirs(join(root, 'resources/skills')).length
  const connectors = listSubdirs(join(root, 'resources/connectors')).length
  const plugins = listSubdirs(join(root, 'resources/plugins'))
  const pluginSkills = plugins.reduce((sum, plugin) => {
    const skillsDir = join(root, 'resources/plugins', plugin, 'skills')
    return sum + (existsSync(skillsDir) ? listSubdirs(skillsDir).filter(s => existsSync(join(skillsDir, s, 'SKILL.md'))).length : 0)
  }, 0)
  const wrappers = countWrapperAdapters(join(root, 'resources/wrappers'))

  assert.equal(provenance.counts.standaloneSkills, skills)
  assert.equal(provenance.counts.connectors, connectors)
  assert.equal(provenance.counts.plugins, plugins.length)
  assert.equal(provenance.counts.pluginSkills, pluginSkills)
  assert.equal(provenance.counts.wrapperAdapters, wrappers)

  assert.match(readmeZh, new RegExp(`\\| MCP connector 定义 \\| ${connectors} \\|`))
  assert.match(readmeZh, new RegExp(`\\| Phi 复合插件 \\| ${plugins.length} \\|`))
  assert.match(readmeZh, new RegExp(`\\| Wrapper 适配器 \\| ${wrappers} \\|`))
  assert.match(readmeEn, new RegExp(`\\| MCP connector definitions \\| ${connectors} \\|`))
  assert.match(readmeEn, new RegExp(`\\| Composite Phi plugins \\| ${plugins.length} \\|`))
  assert.match(readmeEn, new RegExp(`\\| Wrapper adapters \\| ${wrappers} \\|`))

  const addedByPath = new Map((provenance.addedContent ?? []).map(item => [item.path.replace(/\/$/, ''), item]))
  for (const plugin of plugins) {
    const pluginDir = join(root, 'resources/plugins', plugin)
    const manifestPath = join(pluginDir, 'phi-package.yaml')
    assert.ok(existsSync(manifestPath), `${plugin} missing phi-package.yaml`)
    assert.ok(existsSync(join(pluginDir, 'README.md')), `${plugin} missing README.md`)
    assert.ok(existsSync(join(pluginDir, 'CHANGELOG.md')), `${plugin} missing CHANGELOG.md`)
    assert.ok(readdirSync(pluginDir).some(f => /^icon\.(svg|png|webp|jpg|jpeg)$/.test(f)), `${plugin} missing icon`)
    if (plugin !== 'visualization') {
      const entry = addedByPath.get(`resources/plugins/${plugin}`)
      assert.ok(entry, `${plugin} missing in SOURCE.json addedContent`)
      const manifestVer = readFileSync(manifestPath, 'utf8').match(/^version:\s*(\S+)$/m)?.[1]
      assert.equal(entry.packageVersion, manifestVer, `${plugin} version mismatch between phi-package.yaml and SOURCE.json`)
    }
  }
})

test('content source omits Phi-owned shared runtime and palette trees', () => {
  for (const relative of ['resources/palettes', 'resources/runtime']) {
    assert.equal(existsSync(join(root, relative)), false, `${relative} belongs to the Phi engine`)
  }
})

test('visualization retains the palette data its installed helpers read', () => {
  const skillRoot = join(root, 'resources/plugins/visualization/skills/omics-visualization')
  assert.ok(readFileSync(join(skillRoot, 'references/palettes.yaml'), 'utf8').length > 0)
  const catalog = JSON.parse(readFileSync(join(skillRoot, 'references/palettes/colors.json'), 'utf8'))
  assert.ok(catalog && typeof catalog === 'object')
})

const family = join(root, 'resources/wrappers/modules/local/differential-expression')
const dependencies = {
  deseq2: ['bioconductor-deseq2', 'bioconductor-limma', 'bioconductor-biocparallel', 'r-ashr'],
  edger: ['bioconductor-edger', 'bioconductor-limma'],
  limma: ['bioconductor-edger', 'bioconductor-limma'],
  qc: ['bioconductor-deseq2', 'r-ggplot2', 'r-pheatmap'],
  timeseries: ['bioconductor-variancepartition', 'bioconductor-edger', 'bioconductor-limma', 'bioconductor-biocparallel', 'r-lme4'],
  'visualization/pca': ['bioconductor-deseq2', 'r-ggplot2'],
  'visualization/heatmap': ['r-pheatmap'],
  'visualization/volcano': ['r-ggplot2']
}

test('retired differential-expression image source is absent', () => {
  assert.equal(existsSync(join(root, 'resources/wrappers/images/differential-expression-r')), false)
})

for (const [name, required] of Object.entries(dependencies)) {
  test(`${name} keeps executable code and declares its required R dependencies locally`, () => {
    const moduleRoot = join(family, name)
    const source = readFileSync(join(moduleRoot, 'main.nf'), 'utf8')
    assert.ok(!/phi\/differential-expression-r|images\/differential-expression-r/.test(source), 'retired image reference remains')
    const reference = source.match(/conda\s+["']\$\{moduleDir\}\/([^"']+)["']/)
    assert.ok(reference, 'module needs a package-relative Conda environment')
    assert.equal(reference[1].includes('..'), false, 'environment must stay with its module')
    const spec = readFileSync(join(moduleRoot, reference[1]), 'utf8')
    assert.match(spec, /^dependencies:\s*$/m)
    const packages = new Set([...spec.matchAll(/^\s*-\s*(?:[\w-]+::)?([\w-]+)(?:[=<>][^\n#]+)?\s*(?:#.*)?$/gm)].map(match => match[1]))
    for (const pkg of required) assert.ok(packages.has(pkg), `missing required dependency ${pkg}`)
    assert.ok(existsSync(join(moduleRoot, 'wrapper/wrapper.yaml')), 'keep the callable wrapper')
    assert.ok(existsSync(join(moduleRoot, 'wrapper/main.nf')), 'keep the wrapper entrypoint')
  })
}
